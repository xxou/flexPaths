#' Path-Specific Mean Estimation for a Single Binary Treatment
#'
#' Estimates a path-specific counterfactual mean for a single binary treatment
#' and causally ordered mediators. The \code{assign} vector specifies treatment
#' values on the edges into each mediator stage and the outcome, encoding the
#' pathways of interest. Under the identification assumptions, the target mean
#' is identified by the edge g-formula. Available estimators are inverse
#' probability weighting (IPW), sequential regression-based g-computation (G),
#' and efficient influence function-based estimation (EIF).
#'
#' @param pathsInfo An object of class \code{pathsInfo} containing the analysis settings
#' and model specifications provided by the \code{pathsInfo} function.
#' The analysis data are taken from \code{pathsInfo$data}.
#' @param assign A numeric vector containing only 0 and 1, without missing values, specifying the treatment values for each mediator stage (M1 to Mk) and the outcome (Y) in order. Integer vectors are also accepted.
#' The last value in the vector corresponds to the treatment value for Y. The length of this vector should be equal to the number of mediators plus one (for the direct effect on Y).
#' @param nboot An optional integer of at least 2 specifying the number of bootstrap
#' samples used to estimate the standard error. If \code{NULL}, no bootstrapping is performed.
#'
#' @return An object of class \code{flexMean}, which includes the following components:
#' \describe{
#'   \item{\code{call}}{The matched call to the \code{flexMean} function.}
#'   \item{\code{pathsInfo}}{The original \code{pathsInfo} object used for the analysis.}
#'   \item{\code{assign}}{A vector representing the assigned treatment values used for each mediator and the direct outcome.}
#'   \item{\code{results}}{A data frame containing the calculated potential mean (value).
#'   The standard error (SE) is included only for EIF estimation; G and IPW omit this column.}
#'   \item{\code{potential_data}}{A vector containing the calculated potential outcomes for each observation in the dataset.}
#'   \item{\code{boot_results}}{If bootstrapping is performed, a data frame containing
#'   the treatment assignments, the original potential mean (value), its bootstrap standard
#'   error (boot.SE), and the number of bootstrap samples (nboot).}
#' }
#'
#' @details Bootstrapping resamples individuals with replacement and refits all
#' required models in each sample. It is available for EIF, G, and IPW estimation.
#' The bootstrap standard error is the sample standard deviation of the bootstrap means.
#'
#' @import stats purrr
#' @import SuperLearner dbarts
#'
#' @export
#'
#' @example inst/examples/flexMean-example.R
#'
flexMean <- function(
  pathsInfo, assign, nboot = NULL  # M1, M2, M3,...,Direc
) {
  cl <- match.call()
  check_nboot(nboot)
  # Variables populated by list2env().
  data <- A <- Y <- cov_x <- M.list <- estimation <- NULL
  Pmodel.lists <- Omodel.lists <- Imodel.lists <- Imodel.source <- NULL
  list2env(pathsInfo, envir = environment())

  #### 1. data preparation ---------
  # A is a binary;
  I = data[[A]]
  # counterfactual data
  A0_data <- A1_data <- data;
  A0_data[[A]] <- 0; A1_data[[A]] <- 1

  ### cov_x and mediators
  Km <- length(M.list)    # have Km mediators
  if ((Km + 1) != length(assign)) {
    stop("Error: The length of 'assign' must be equal to the length of 'M.list' + 1.
         Please ensure both inputs are of the same length.")
  }
  if (!is.numeric(assign) || is.complex(assign) || !is.null(dim(assign)) ||
    anyNA(assign) || !all(assign %in% c(0, 1))) {
    stop("'assign' must be a numeric vector containing only 0 and 1, without missing values.")
  }

  # combine the same value of path
  M.index <- create_index(assign)[-(Km + 1)]

  M.list_s <- combine_item(M.list, M.index)
  assign_s <- rle(assign)$values
  cum_mediators <- Reduce(c, c(list(cov_x), M.list_s), accumulate = TRUE)

  K = length(cum_mediators) # K-1 combined mediators and cov


  #### 2. propensity model ----------
  if (estimation != "G") {
    Pmodel = Pmodel.lists[c(1, unique(M.index[M.index != 0]) + 1)]

    pro.details <- map(Pmodel, extract_model_details)
    pro.fit_names <- map(pro.details, ~ paste0("fl_", .x$fit_name))
    pro.calls <- map(pro.details, "model_call")


    pro.modelfits <- pmap(list(cum_mediators, pro.fit_names, pro.calls),
      function(cum_mediators, pro.fit_names, pro.calls) {
        model_fun <- match.fun(pro.fit_names)
        model <- model_fun(data, X = cum_mediators, Y = A, fl_call = pro.calls)
        return(model)
      })

  # P(A=1|cov_X); P(A=1|M1,cov_X); P(A=1|M1,...,MK,X)
    gA1_Mk_X.list <- map2(cum_mediators, pro.modelfits, ~ {predict(.y, newdata = data[, .x, drop = FALSE])})
  # P(A=0|cov_X); P(A=0|M1,cov_X); P(A=0|M1,...,MK,X)
    gA0_Mk_X.list <- map(gA1_Mk_X.list, ~ {1 - .x})
  }

  #### 3. outcome model modelO.list ----------
  # check only one
  if (estimation != "IPW") {
    Omodel <- Omodel.lists[max(c(1, unique(M.index[M.index != 0]) + 1))]
    out.details <- map(Omodel, extract_model_details)
    out.fit_names <- map(out.details, ~ paste0("fl_", .x$fit_name))
    out.calls <- map(out.details, "model_call")

    Imodel <- Imodel.lists[c(1, unique(M.index[M.index != 0])[-(K - 1)] + 1)]
    iter.details <- map(Imodel, extract_model_details)
    iter.fit_names <- map(iter.details, ~ paste0("fl_", .x$fit_name))
    iter.calls <- map(iter.details, "model_call")
    if (Imodel.source == "default") {iter.calls <- map(iter.calls, ~ replace_family(.x))}

    fit_names <- c(iter.fit_names, out.fit_names)
    mu.calls <- c(iter.calls, out.calls)


  # Innermost function
    mu.all <- matrix(NA, ncol = K + 1, nrow = nrow(data)) # EE..E()|A,X)|..), ...., E(Y|mk,A,x), Y
    mu.all[, K + 1] <- data[, Y]

  # iterative outcome model
    for (i in K:1) {
      data$mu <- mu.all[, i + 1]

      model_fun <- match.fun(fit_names[[i]])
      modelfit <- model_fun(data, X = c(A, cum_mediators[[i]]), Y = "mu", fl_call = mu.calls[[i]])

      Ai_data <- if (assign_s[i] == 1) { A1_data } else { A0_data }
      mu.all[, i] <- predict(modelfit, Ai_data[, c(A, cum_mediators[[i]]), drop = FALSE])
    }
  }
  #### 4. calculate phi() and Pn( ) ------
  if (estimation == "EIF") {
  # phi() is a value for each row of EIF formula
    phi.all <- matrix(NA, ncol = K + 1, nrow = nrow(data)) # calculation from muK to mu1
    phi.all[, 1] <- mu.all[, 1]


    for (i in K:1) {
      if (assign_s[i] == 1) {
        gA_X <- gA1_Mk_X.list[[1]]
        ratio_index = -1
      } else {
        gA_X <- gA0_Mk_X.list[[1]]
        ratio_index = 1
      }

    # Bayes density ratio
      product <- 1 # Initialize the product
    # Iterate through previous indices j < i
      for (j in seq_len(i - 1)) {
        if (assign_s[i] != assign_s[j]) {
          product <- product *
            ((gA1_Mk_X.list[[j + 1]] / gA0_Mk_X.list[[j + 1]]) * (gA0_Mk_X.list[[j]] / gA1_Mk_X.list[[j]])) ^ (ratio_index)
        }
      }


    # sub EIF function for each row
      phi.all[, i + 1] = (I == assign_s[i]) / gA_X * product * (mu.all[, i + 1] - mu.all[, i])
    }
    potential_data = rowSums(phi.all)
  }


  if (estimation == "G") {
    potential_data = mu.all[, 1]
  }

  if (estimation == "IPW") {
    if (assign_s[K] == 1) { gA_X <- gA1_Mk_X.list[[1]]; ratio_index = -1
    } else { gA_X <- gA0_Mk_X.list[[1]]; ratio_index = 1 }

    # Bayes density ratio
    product <- 1 # Initialize the product
    # Iterate through previous indices j < i
    for (j in seq_len(K - 1)) {
      if (assign_s[K] != assign_s[j]) {
        product <- product *
          ((gA1_Mk_X.list[[j + 1]] / gA0_Mk_X.list[[j + 1]]) * (gA0_Mk_X.list[[j]] / gA1_Mk_X.list[[j]])) ^ (ratio_index)
      }
    }
    potential_data = (I == assign_s[K]) / gA_X * product * data[, Y]
  }

  results = data.frame(assign = paste0(as.character(assign), collapse = ""),
    value = mean(potential_data, na.rm = T)
  )
  if (estimation == "EIF") {
    results$SE = sqrt(var(potential_data, na.rm = T) / nrow(data))
  }

  output <- list(call = cl,
    pathsInfo = pathsInfo,
    assign = assign,
    results = results,
    potential_data = potential_data)
  if (!is.null(nboot)) {
    boot.values <- map_dbl(seq_len(nboot),
      ~ one_boot.potential(pathsInfo, assign, flexMean))
    output$boot_results <- data.frame(assign = results$assign, value = results$value,
      boot.SE = sd(boot.values), nboot = nboot)
  }
  class(output) <- 'flexMean'
  return(output)
}
