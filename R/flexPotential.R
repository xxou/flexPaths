#' Causal Path-Specific Potential Outcome Calculation for Single Treatment
#'
#' This function calculates causal path-specific potential outcomes for a single binary treatment.
#' Users can specify different levels of the treatment for each mediator (M1, M2, ..., Mk) and directly define the outcome (Y) to represent the potential outcome for a specific path or a combination of paths.
#' The function supports various estimation methods, including EIF (Efficient Influence Function), IPW (Inverse Probability Weighting), and G-computation.
#'
#' @param data A data frame containing all the variables required for the analysis.
#' @param pathsFit An object of class \code{pathsFit}, which is the output from the \code{pathsFit} function, containing the necessary fitted models and data.
#' @param active A vector specifying the active values for each mediator (M1 to Mk) and the direct value setting for the outcome (Y) in order.
#' The last value in the vector corresponds to the treatment value for Y. The length of this vector should be equal to the number of mediators plus one (for the direct effect on Y).
#'
#' @return An object of class \code{flexPotential}, which includes the following components:
#' \describe{
#'   \item{\code{call}}{The matched call to the \code{flexPotential} function.}
#'   \item{\code{pathsFit}}{The original \code{pathsFit} object used for the analysis.}
#'   \item{\code{active}}{A vector representing the active values used for each mediator and the direct outcome.}
#'   \item{\code{results}}{A data frame containing the calculated potential outcome (value) and its standard error (SE).}
#'   \item{\code{potential_data}}{A vector containing the calculated potential outcomes for each observation in the dataset.}
#' }
#'
#' @examples
#' # Example usage
#' data("singTreat")
#'
#' fit <- pathsFit(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
#'                 M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
#'                 estimation = "EIF",
#'                 model.outcome = list(~ glm(family = gaussian())),
#'                 model.propensity = ~ bart(verbose = FALSE, ndpost = 200))
#'
#' potential_outcome <- flexPotential(
#'   pathsFit = fit,
#'   active = c(1, 0, 1, 1)  # Setting active values for M1, M2, M3, and Y
#' )
#'
#' @import stats purrr
#' @import SuperLearner dbarts
#'
#' @export
#'
#'
flexPotential <-  function(
    pathsFit,active  # M1, M2, M3,...,Direc
){
  cl <- match.call()
  list2env(pathsFit, envir = environment())

  #### 1. data preparation ---------
  # A is a binary;
  I = data[[A]]
  # counterfactual data
  A0_data <- A1_data <- data;
  A0_data[[A]] <- 0; A1_data[[A]] <- 1

  ### cov_x and mediators
  Km <- length(M.list)    # have Km mediators
  if((Km + 1)!=length(active)){
    stop("Error: The length of 'active' must be equal to the length of 'M.list' + 1.
         Please ensure both inputs are of the same length.")
  }

  # combine the same value of path
  M.index <- create_index(active)[-(Km + 1)]

  M.list_s <- combine_item(M.list, M.index)
  active_s <- rle(active)$values
  cum_mediators <- Reduce(c, c(list(cov_x), M.list_s), accumulate = TRUE)

  K = length(cum_mediators) # K-1 combined mediators and cov


  #### 2. propensity model ----------
  if(estimation != "G"){
  Pmodel = Pmodel.lists[c(1,unique(M.index[M.index!=0])+1)]

  pro.details<- map(Pmodel,extract_model_details)
  pro.fit_names <- map(pro.details, ~ paste0("fl_", .x$fit_name))
  pro.calls <- map(pro.details, "model_call")


  pro.modelfits<- pmap(list(cum_mediators, pro.fit_names, pro.calls),
                       function(cum_mediators, pro.fit_names, pro.calls){
                         model_fun <- match.fun(pro.fit_names)
                         model<-model_fun(data, X = cum_mediators,Y = A, fl_call=pro.calls)
                         return(model)
                       })

  # P(A=1|cov_X); P(A=1|M1,cov_X); P(A=1|M1,...,MK,X)
  gA1_Mk_X.list <- map2(cum_mediators, pro.modelfits, ~ {predict(.y, newdata = data[,.x, drop=FALSE])})
  # P(A=0|cov_X); P(A=0|M1,cov_X); P(A=0|M1,...,MK,X)
  gA0_Mk_X.list <- map(gA1_Mk_X.list, ~ {1-.x})
}

  #### 3. outcome model modelO.list ----------
  # check only one
  if(estimation !="IPW"){
  Omodel <-  Omodel.lists[max(c(1,unique(M.index[M.index!=0])+1))]
  out.details<-  map(Omodel,extract_model_details)
  out.fit_names <- map(out.details, ~ paste0("fl_", .x$fit_name))
  out.calls <- map(out.details, "model_call")

  Imodel <- Omodel.lists[c(1,unique(M.index[M.index!=0])+1)][-K]
  iter.details<- map(Imodel,extract_model_details)
  iter.fit_names <- map(iter.details, ~ paste0("fl_", .x$fit_name))
  iter.calls <- map(iter.details, "model_call")
  iter.calls <-map(iter.calls, ~ replace_family(.x))

  fit_names <- c(iter.fit_names, out.fit_names)
  mu.calls <- c(iter.calls, out.calls)


  # Innermost function
  mu.all <- matrix(NA,nc=K+1,nr=nrow(data)) # EE..E()|A,X)|..), ...., E(Y|mk,A,x), Y
  mu.all[,K+1] <- data[,Y]

  # iterative outcome model
  for(i in K:1){
    data$mu <- mu.all[,i+1]

    model_fun <- match.fun(fit_names[[i]])
    modelfit <- model_fun(data, X = c(A,cum_mediators[[i]]),Y = "mu", fl_call=mu.calls[[i]])

    Ai_data <- if(active_s[i]==1){ A1_data }else{ A0_data }
    mu.all[,i] <- predict(modelfit, Ai_data[,c(A,cum_mediators[[i]])])
  }
}
  #### 4. calculate phi() and Pn( ) ------
  if(estimation == "EIF"){
  # phi() is a value for each row of EIF formula
  phi.all <- matrix(NA,nc=K+1,nr=nrow(data)) # calculation from muK to mu1
  phi.all[,1] <- mu.all[,1]


  for (i in K:1) {
    if(active_s[i] == 1){
      gA_X <- gA1_Mk_X.list[[1]]
      ratio_index = -1
    }else{
      gA_X <- gA0_Mk_X.list[[1]]
      ratio_index = 1
    }

    # Bayes density ratio
    product <- 1 # Initialize the product
    # Iterate through previous indices j < i
    for (j in seq_len(i - 1)) {
      if (active_s[i] != active_s[j]) {
        product <- product *
          ( (gA1_Mk_X.list[[j+1]]/gA0_Mk_X.list[[j+1]])*(gA0_Mk_X.list[[j]] / gA1_Mk_X.list[[j]]) )^(ratio_index)
      }
    }


    # sub EIF function for each row
    phi.all[,i+1] = (I == active_s[i])/gA_X*product*(mu.all[,i+1]-mu.all[,i])
  }
   potential_data = rowSums(phi.all)
  }


  if(estimation == "G"){
    potential_data = mu.all[,1]
  }

  if(estimation =="IPW"){
    if(active_s[K] == 1){ gA_X <- gA1_Mk_X.list[[1]]; ratio_index = -1
    }else{ gA_X <- gA0_Mk_X.list[[1]]; ratio_index = 1 }

    # Bayes density ratio
    product <- 1 # Initialize the product
    # Iterate through previous indices j < i
    for (j in seq_len(K - 1)) {
      if (active_s[K] != active_s[j]) {
        product <- product *
          ( (gA1_Mk_X.list[[j+1]]/gA0_Mk_X.list[[j+1]])*(gA0_Mk_X.list[[j]] / gA1_Mk_X.list[[j]]) )^(ratio_index)
      }
    }
    potential_data = (I == active_s[K])/gA_X*product*data[,Y]
  }

  results = data.frame(active = paste0(as.character(active), collapse = ""),
                   value = mean(potential_data,na.rm=T),
                   SE = sqrt(var(potential_data,na.rm=T)/nrow(data))
  )

  output <- list(call = cl,
                 pathsFit =pathsFit,
                 active = active,
                 results = results,
                 potential_data = potential_data)
  class(output) <- 'flexPotential'
  return(output)
}
