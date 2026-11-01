#' Causal Path-Specific Potential Outcome Calculation for Multiple Treatments
#'
#' This function calculates causal path-specific potential outcomes for multiple treatments.
#' Users can specify different levels of each treatment for each mediator (M1 to Mk) and directly define the outcome (Y) to represent the potential outcome for a single path or several combined paths.
#' The function supports various estimation methods, including EIF (Efficient Influence Function), IPW (Inverse Probability Weighting), and G-computation.
#'
#' @param pathsInfo An object of class \code{pathsInfo}, returned by \code{pathsInfo()},
#' containing the analysis settings and model specifications. The analysis data
#' are taken from \code{pathsInfo$data}.
#' @param active A list of numeric vectors, with one vector per treatment in the order of \code{pathsInfo$A}. Integer vectors are also accepted.
#' Each vector must have length \code{length(pathsInfo$M.list) + 1}, with entries for mediator stages M1 to Mk followed by the outcome Y.
#' Values must be 0 or 1. Use \code{NA} for mediator stages that occur before the corresponding treatment; \code{NaN} is not allowed.
#'
#' @return An object of class \code{mflexPotential}, which includes the following components:
#' \describe{
#'   \item{\code{call}}{The matched call to the \code{mflexPotential} function.}
#'   \item{\code{pathsInfo}}{The original \code{pathsInfo} object used for the analysis.}
#'   \item{\code{active}}{A matrix representing the active values used for each mediator and the direct outcome.}
#'   \item{\code{results}}{A data frame containing the average potential outcome (value) and its standard error (SE).}
#'   \item{\code{potential_data}}{A vector containing the calculated potential outcomes for each observation in the dataset.}
#' }
#'
#'
#' @import stats purrr
#' @import SuperLearner dbarts
#'
#' @export
#'
#' @example inst/examples/mflexPotential-example.R
#'
#'
mflexPotential <-  function(
   pathsInfo,active
){
  cl <- match.call()
  # Variables populated by list2env().
  data <- A <- Y <- cov_x <- M.list <- estimation <- NULL
  Pmodel.lists <- Omodel.lists <- Imodel.lists <- Imodel.source <- NULL
  list2env(pathsInfo, envir = environment())

  # **********************************
  #    1. data preparation ---------
  # **********************************
  ## cov_x and mediators
  N = nrow(data)
  Km <- length(M.list)

  if(!is.list(active) || !is.null(dim(active))){
    stop("'active' must be a list of numeric vectors.")
  }
  if(length(active) != length(A)){
    stop("The length of 'active' must equal the number of treatments.")
  }
  for(i in seq_along(active)){
    a <- active[[i]]
    if(length(a) != Km + 1){
      stop("Each vector in 'active' must have length equal to the length of 'M.list' + 1.")
    }
    if(!is.numeric(a) || is.complex(a) || !is.null(dim(a)) ||
       !all(a %in% c(0, 1, NA))){
      stop(sprintf("'active[[%d]]' must be a numeric vector containing only 0, 1 or NA.", i))
    }
  }

  # combine mediators
  # mediators for model
  active_m <- do.call(rbind, active)
  M.index <- mcreate_index(active_m)[-(Km + 1)]
  M.list_s <- combine_item(M.list, M.index)
  cum_mediators <- cumulate_vars(c(list(cov_x),M.list_s))

  # treatment value for model
  active_s <- simplify_active(active_m)
  # treatment variables list; assigned value list
  a.list <- map(seq_len(ncol(active_s)), ~ A[which(!is.na(active_s[, .x]))])
  active_s.list<-map(1:ncol(active_s), ~ as.numeric(na.omit(active_s[, .x])))

  # cumulative expression for model
  cum_a.list <- map(a.list,~cumulate_vars(.x,first_null = T))
  cum_active_s.list<-map(active_s.list,~cumulate_vars(.x,first_null = T))

  K = length(cum_mediators)

  # **********************************
  #     2. outcome model  ----------
  # **********************************
  if(estimation != "IPW"){
  out.details <- map(Omodel.lists, extract_model_details)
  out.fit_names <- map(out.details, ~ paste0("fl_", .x$fit_name))
  out.calls <-  map(out.details, "model_call")

  iter.details <- map(Imodel.lists, extract_model_details)
  iter.fit_names <- map(iter.details, ~ paste0("fl_", .x$fit_name))
  iter.calls <- map(iter.details, "model_call")
  if(Imodel.source == "default"){iter.calls <- map(iter.calls, ~ replace_family(.x))}

  fit_names <- c(rep(iter.fit_names, K-1), out.fit_names)
  mu.calls <- c(rep(iter.calls, K-1), out.calls)

  # Innermost function
  mu.all <- matrix(NA, ncol = K+1, nrow = N) # EE..E()|A,X)|..), ...., E(Y|mk,A,x), Y
  mu.all[ , K+1] <- data[ , Y]


  # iterative outcome model
  for(i in K:1){
    data$mu <- mu.all[ , i+1]
    Ai_data = data
    # counterfactual data according to assigned treatment values
    Ai_data[ , a.list[[i]]] = suppressMessages({map_dfc(active_s.list[[i]], ~rep(.x, N))})


    model_fun <- match.fun(fit_names[[i]])
    modelfit <- model_fun(data, X = c(a.list[[i]], cum_mediators[[i]]), Y = "mu", fl_call = mu.calls[[i]])

    # cat("num:", i, "\n"); print(a.list[[i]]); print(active_s.list[[i]]); print(modelfit$model)
    mu.all[ , i] <- predict(modelfit, Ai_data[ , c(a.list[[i]], cum_mediators[[i]]), drop = FALSE])
  }

  }

  # **********************************
  #    3. propensity model ----------
  # **********************************
  if(estimation !="G"){
  pro.details <- extract_model_details(Pmodel.lists[[1]])
  pro.fit_names <-paste0("fl_", pro.details$fit_name)
  pro.calls <- pro.details$model_call


  # Nested list: propensity for p(t1|...);p(t2|t1, ...);p(t3|...)

  pro.modelfits<-pmap(list(cum_mediators, cum_a.list, a.list),
                      function(mediator, a_predictor, a_out){

                        model_fun <- match.fun(pro.fit_names)
                        modelfit <- map2(a_predictor, a_out,
                                         ~{model_fun(data, X = c(mediator, .x), Y = .y, fl_call = pro.calls)})

                        # print(mediator); cat("____\n");
                        # print(a_predictor); cat("____\n");
                        # print(a_out); cat("____\n");
                        # print(map(modelfit, "model"))
                        # cat("\nNext----------------------------------------------\n")

                        return(modelfit)
                      })



  # propensity functions in the denominators of Indicator
  I.model.index<- map_int(unique(flatten(cum_a.list) ), ~ which(map_lgl(flatten(cum_a.list), identical, .x))[1])
  model.I <- flatten(pro.modelfits)[I.model.index]  # print(map(model.I, "model"))
  model.I.list <-map(seq_len(ncol(active_s)), ~ model.I[which(!is.na(active_s[, .x]))])
  }

  # print(map(model.I.list[[5]], "model"))


  # **********************************
  # 4. potential outcome data --------
  # **********************************
  if(estimation == "EIF"){
  phi.all <- matrix(NA, ncol = K+1, nrow = N) # calculation from muK to mu1
  phi.all[ ,1] <- mu.all[ ,1]

  for (i in K:1) {
    active_i <- active_s.list[[i]]

    ## indicator I(t1 = , t2 =, ...)
    Indicator_expr <- paste(map2(a.list[[i]], active_i, ~ paste0(.x, "==", .y)),
                            collapse = " & ")
    I = eval(parse(text = Indicator_expr), envir = data)


    ## product.I:
    # cat(i, ": product.I \n")
    product.I = prod_assign_modellist(model_list = model.I.list[[i]],
                                      a_vars = a.list[[i]],
                                      a_values = active_s.list[[i]] )



    ## Bayes density ratio
    # cat(i, ": Bayes density ratio \n")

    product.Bayes <- 1 # Initialize the product
    # Iterate through previous indices j < i
    for (j in seq_len(i - 1)) {
      if (!all(active_s[,i] == active_s[,j], na.rm = T) ) {
        product.Bayes <- product.Bayes *
          (  prod_assign_modellist(pro.modelfits[[j+1]], a.list[[j]], active_s.list[[j]]) /  # use the value of itself to predict
               prod_assign_modellist(pro.modelfits[[j+1]], a.list[[j]], active_s.list[[i]])    # use the value of i_th to predict
          )*
          (prod_assign_modellist(pro.modelfits[[j]],a.list[[j]],active_s.list[[i]])/
             prod_assign_modellist(pro.modelfits[[j]],a.list[[j]],active_s.list[[j]]))
      }
    }


    # sub EIF function for each row
    phi.all[ ,i+1] = I/product.I * product.Bayes * ( mu.all[ , i+1] - mu.all[ , i] )

    # cat("-------------------------------")
    }

    potential_data = rowSums(phi.all)
  }

  if(estimation=="G"){
    potential_data = mu.all[,1]
  }

  if(estimation == "IPW"){
      active_K <- active_s.list[[K]]

      ## indicator I(t1 = , t2 =, ...)
      Indicator_expr <- paste(map2(a.list[[K]], active_K, ~ paste0(.x, "==", .y)),
                              collapse = " & ")
      I = eval(parse(text = Indicator_expr),envir = data)


      ## product.I:
      product.I = prod_assign_modellist(model_list = model.I.list[[K]],
                                        a_vars = a.list[[K]],
                                        a_values = active_s.list[[K]] )

      ## Bayes density ratio
      product.Bayes <- 1 # Initialize the product
      # Iterate through previous indices j < i
      for (j in seq_len(K - 1)) {
        if (!all(active_s[,K] == active_s[,j],na.rm = T) ) {
          product.Bayes <- product.Bayes *
            (  prod_assign_modellist(pro.modelfits[[j+1]],a.list[[j]],active_s.list[[j]]) /  # use the value of itself to predict
                 prod_assign_modellist(pro.modelfits[[j+1]],a.list[[j]],active_s.list[[K]])    # use the value of i to predict
            )*
            (prod_assign_modellist(pro.modelfits[[j]],a.list[[j]],active_s.list[[K]])/
               prod_assign_modellist(pro.modelfits[[j]],a.list[[j]],active_s.list[[j]]))
        }
      }


      # sub EIF function for each row
      potential_data = I/product.I*product.Bayes*data[,Y]
  }

  # **********************************
  #       5. results --------------
  # **********************************
  out = data.frame(active = paste0(map_chr(active, ~ paste0(na.omit(.x), collapse = "")),
                                   collapse = ";"),
                   value = mean(potential_data,na.rm=T),
                   SE = sqrt(var(potential_data,na.rm=T)/nrow(data))
                   )

  output <- list(call = cl,
                 pathsInfo = pathsInfo,
                 active = active_m,
                 results = out,
                 potential_data = potential_data)
  class(output) <- 'mflexPotential'
  return(output)

}



