#' Causal Path-Specific Analysis Setup
#'
#' This function sets up the necessary components for conducting a causal path-specific analysis.
#' It prepares the data, analysis settings, and model specifications required
#' for further mean and effect calculations.
#'
#' @param A A character string or a vector of character strings representing the name(s) of one or more treatments.
#' Each treatment must be a numeric vector containing only 0 and 1 (apart from missing values),
#' or a factor with exactly two non-missing levels, including an ordered factor.
#' Factors are encoded using their existing level order: \code{levels(x)[1]} becomes 0
#' and \code{levels(x)[2]} becomes 1. Set \code{levels} or use \code{relevel} before calling
#' this function to choose the reference level. All treatment columns in the returned
#' \code{data} are converted to numeric 0/1; subsequent \code{assign} values use this encoding.
#' @param Y A character string indicating the name of the outcome variable. Y should be numeric variable.
#' @param cov_x A vector of character strings indicating the name(s) of covariates. At least one covariate is required.
#' @param M.list A list of \eqn{K} mediators, where \eqn{K} is the number of mediators. Each mediator can be multivariate.
#' Use the form \code{list(M1 = c("mediator1_1", "mediator1_2"), M2 = c("mediator2"), ... )}, with \code{M_i} to index the order of mediators.
#' Alternatively, you can list them in order using \code{list(c("mediator1_1", "mediator1_2"), c("mediator2"), ... )}.
#' @param estimation A character string indicating the type of estimation to be used. Choose from "EIF" (Efficient Influence Function),
#' "IPW" (Inverse Probability Weighting), or "G" (G-computation).
#' Direct standard errors are provided only for EIF potential means and EIF effects
#' on the difference scale. Use bootstrapping for uncertainty estimation otherwise.
#' @param model.treatment Not used for "G" estimation. A model specification for
#' the conditional treatment probabilities \eqn{P(A \mid \ldots)} required by
#' the estimator. Supply one formula or a list containing one formula to apply
#' the same fitting specification to all required treatment regressions.
#' This applies to both single and repeated treatments.
#' \itemize{
#'   \item For a single treatment only, a list of \eqn{K+1} formulas may instead
#'   specify \eqn{P(A \mid X)}, \eqn{P(A \mid M_1, X)}, ...,
#'   \eqn{P(A \mid M_1, \ldots, M_K, X)}, in that order. The list may be named
#'   \code{cov_x}, \code{M1}, ..., \code{MK}. For repeated treatments, use one
#'   common formula. A common specification is recommended; see Details.
#'   \item Formulas may use \code{glm}, \code{bart}, or \code{SuperLearner}.
#'   Supply the fitting arguments; the package supplies the response and predictors.
#'   To create a custom model, copy the template located by:
#'   \preformatted{
#'   system.file("templates", "fl_model_template.R",
#'               package = "flexPaths")
#'   }
#'   Adapt the fitting and prediction functions, and source it into the global environment.
#' }
#' @param model.outcome Not used for "IPW" estimation. Supply one formula or a
#' list containing one formula for outcome regression. For repeated treatments,
#' this common fitting specification is used with the required predictor sets.
#' For a single treatment only, a list of \eqn{K+1} formulas may instead specify
#' \eqn{E(Y \mid A, X)}, \eqn{E(Y \mid A, M_1, X)}, ...,
#' \eqn{E(Y \mid A, M_1, \ldots, M_K, X)}, in that order, with the same naming
#' convention as \code{model.treatment}.
#' @param model.iter Not used for "IPW" estimation. Specifies the sequential
#' regressions, which regress fitted values from the next stage on the current
#' treatment and covariate history.
#' \itemize{
#'   \item If \code{NULL}, the corresponding \code{model.outcome} specifications
#'   are reused with the Gaussian family.
#'   \item Supply one formula or a list containing one formula for all sequential
#'   regressions. For repeated treatments, use this common specification.
#'   We recommend this common specification to reduce errors in specifying
#'   and ordering separate regression models.
#'   \item For a single treatment only, a list of \eqn{K} formulas may instead
#'   specify regressions with predictors \eqn{(A, X)}, \eqn{(A, M_1, X)}, ...,
#'   \eqn{(A, M_1, \ldots, M_{K-1}, X)}, in that order.
#' }
#' @param data A wide-format data frame with one row per individual. Treatments and mediators measured at different times must be stored in separate columns. Treatment types, levels, and non-missing values are checked before rows with missing values are dropped. A summary of dropped rows and final sample size will be provided. The input data frame is not modified.
#'
#' @details We recommend a common \code{SuperLearner} specification within each
#' of \code{model.treatment} and \code{model.outcome}, with an appropriate
#' \code{family} and a library of candidate learners, to improve prediction and
#' reduce reliance on a single parametric model. Using one specification for
#' each argument is simpler than specifying every regression separately.
#'
#' @return If \code{cov_x} is missing or empty, a warning is issued and \code{NULL} is returned.
#' Otherwise, an object of class \code{pathsInfo}, which is a list containing the following elements:
#' \describe{
#'   \item{\code{data}}{The data frame used for the analysis, with missing rows removed and all treatment columns encoded as numeric 0/1. Other columns retain their original types and values for the retained rows.}
#'   \item{\code{A}}{A vector or list representing the treatment(s) in the causal analysis.}
#'   \item{\code{Y}}{A character representing the outcome variable.}
#'   \item{\code{cov_x}}{A vector of covariates to be adjusted for in the models.}
#'   \item{\code{M.list}}{A list of mediators in the causal pathway, reordered if necessary.}
#'   \item{\code{estimation}}{A character string indicating the type of estimation to be used, such as "EIF" (Efficient Influence Function), "IPW" (Inverse Probability Weighting), or "G" (G-computation).}
#'   \item{\code{Pmodel.lists}}{A list of formulas for the propensity score models corresponding to each mediator in the causal pathway.}
#'   \item{\code{Omodel.lists}}{A list of formulas for the outcome models corresponding to each mediator in the causal pathway.}
#'   \item{\code{Imodel.lists}}{A list of formulas for potential iterative models used in estimation procedures.}
#'   \item{\code{Imodel.source}}{A character for source of iterative model}

#' }
#'
#' @import stats purrr
#'
#' @export
#'
#' @example inst/examples/pathsInfo-example.R
#'
pathsInfo <- function(
  A, Y, cov_x,
  data,
  M.list = list(),
  estimation = "EIF",       # EIF IPW G
  model.treatment = NULL,  # P(A|X), P(A|M1,X), P(A|M1, M2,X), P(A|M1, M2, M3,X)
  model.outcome = NULL,     # E(Y|A,X);  E(Y|M1, A,X); E(Y|M1, M2, A,X)
  model.iter = NULL         # one model, or K-1 models; E(E(Y|M1, A,X)|A,X); E(E(Y|M1, M2, A,X)|M1,A,X)
) {

  # 1. check inputs correct, missing input; correct input type -----------------
  if (!is.data.frame(data)) {stop("'data' must be a data.frame.")}
  if (missing(cov_x) || length(cov_x) == 0) {
    warning("Please provide at least one covariate in 'cov_x'.", call. = FALSE)
    return(NULL)
  }

  if (!is.character(A) || length(A) == 0L) {stop("'A' must be a character string or a non-empty vector of character strings.")}
  if (!all(A %in% names(data))) {stop("All elements of 'A' must exist in the data.")}

  # Validate the original treatment columns, including rows with missing covariates.
  # Keep declared factor levels so the encoding does not change after row removal.
  for (a in A) {
    treatment <- data[[a]]
    if (is.factor(treatment) && is.null(dim(treatment))) {
      if (nlevels(treatment) != 2L || anyNA(levels(treatment))) {
        stop(sprintf("Treatment '%s' must be a factor with exactly two non-missing levels.", a))
      }
    } else if (is.numeric(treatment) && !is.complex(treatment) && is.null(dim(treatment))) {
      if (any(!is.na(treatment) & !(treatment %in% c(0, 1)))) {
        stop(sprintf("Numeric treatment '%s' must contain only 0 and 1 (apart from missing values).", a))
      }
    } else {
      stop(sprintf("Treatment '%s' must be numeric 0/1 or a factor with exactly two non-missing levels.", a))
    }
  }

  na_rows <- which(rowSums(is.na(data)) > 0)
  # Report and drop rows with NA
  if (length(na_rows) > 0) {
    message(sprintf("Dropping %d rows due to missing values. Final sample size: %d.",
      length(na_rows), nrow(data) - length(na_rows)))
    data <- data[-na_rows, , drop = FALSE]
  } else {
    message(sprintf("No missing values found. Final sample size: %d.", nrow(data)))
  }

  if (!is.character(Y) || length(Y) != 1) {stop("'Y' must be a single character string.")}
  if (!(Y %in% names(data))) { stop("'Y' must exist in the data.")}
  if (!(is.numeric(data[[Y]]))) {
    stop("'Y' must be numeric.")
  }


  if (!is.character(cov_x)) {stop("'cov_x' must be a character string or a vector of character strings.")}
  if (!all(cov_x %in% names(data))) {stop("All elements of 'cov_x' must exist in the data.")}

  if (!is.list(M.list)) {stop("'M.list' must be a list.")}
  if (!all(unlist(M.list) %in% names(data))) { stop("All mediators in 'M.list' must exist in the data.")}

  if (!estimation %in% c("EIF", "IPW", "G")) {stop("'estimation' must be one of 'EIF', 'IPW', or 'G'.")}

  for (a in A) {
    treatment <- data[[a]]
    data[[a]] <- if (is.factor(treatment)) as.numeric(treatment) - 1 else as.numeric(treatment)
  }


  # 2. reorder M.list, model.list --------------------------------
  M.list <- reorder_list(M.list)
  Nk <- length(M.list) + 1


  # 3. check model list input
  # ensure all input model as form of list
  model.treatment = if (purrr::is_formula(model.treatment)) list(model.treatment) else model.treatment
  model.outcome = if (purrr::is_formula(model.outcome)) list(model.outcome) else model.outcome
  model.iter = if (purrr::is_formula(model.iter)) list(model.iter) else model.iter
  # a list or NULL

  # models for further effect calculation
  Pmodel.lists = Omodel.lists = Imodel.lists = NULL
   ## (1) single treatment scenario ------
  if (length(A) == 1) {
    type = "single treatment"
    # number of model.treatment
    if (length(model.treatment) == 1) {
      # a list of one formula --> get list of Nk models
      Pmodel.lists = rep(model.treatment, Nk)
    } else if (length(model.treatment) == Nk) {
      # a list of Nk model
      Pmodel.lists = reorder_list(model.treatment)
    }

    # number of model.outcome
    if (length(model.outcome) == 1) {
      # a list of one formula --> get list of Nk models
      Omodel.lists = rep(model.outcome, Nk)
    } else if (length(model.outcome) == Nk) {
      # a list of Nk model
      Omodel.lists = reorder_list(model.outcome)
    }

    # number of model.iter
    Imodel.source = "user input"
    if (is.null(model.iter)) {
      Imodel.source = "default"
      Imodel.lists = Omodel.lists[-Nk]
      # a list of (Nk - 1) model from Omodel.lists or NULL
    } else if (length(model.iter) == 1) {
      Imodel.lists = rep(model.iter, Nk - 1)
    } else if (length(model.iter) == Nk - 1) {
      Imodel.lists = reorder_list(model.iter)
    }

  }

  ## (2) multiple treatments scenario ------
  if (length(A) > 1) {
    type = "multiple treatments"

    if (length(model.treatment) == 1) {
      Pmodel.lists = model.treatment
    }

    if (length(model.outcome) == 1) {
      Omodel.lists = model.outcome
    }
    #
    Imodel.source = "user input"
    if (is.null(model.iter)) {
      Imodel.source = "default"
      Imodel.lists = Omodel.lists
    } else if (length(model.iter) == 1) {
      Imodel.lists = model.iter
    }
  }
  # if there is no input or the length is not correct, Pmodel.lists, Omodel.lists, Imodel.lists will be NULL


  if (estimation == "G") {
    if (!is.null(model.treatment)) {
      message("model.treatment will not be used in the plug-in G computation estimation")}
    if (is.null(Omodel.lists) | is.null(Imodel.lists)) {
      stop("Please input 1 or (k+1) model.outcome, k is the number of mediators \n Please input NULL or 1 or k model.iter")}}

  if (estimation == "IPW") {
    if (!is.null(model.outcome)) {
      message("model.outcome will not be used in the IPW estimation")}
    if (is.null(Pmodel.lists)) {
      stop("Please input 1 or (k+1) model.treatment, k is the number of mediators ")}}

  if (estimation == "EIF") {
    if (is.null(Pmodel.lists) | is.null(Omodel.lists) | is.null(Imodel.lists)) {
      if (length(A) == 1) {
        stop("Please input 1 or (k+1) model.treatment, k is the number of mediators \n Please input 1 or (k+1) model.outcome, k is the number of mediators \n Please input NULL or 1 or k model.iter")
      } else {
        stop("Please input 1 model.treatment \n Please input 1 model.outcome \n Please input NULL or 1 model.iter")
      }
    }
  }

 # output for next calculation for paths effect
  output <- list(
    data = data,
    type = type,
    A = A,
    Y = Y,
    cov_x = cov_x,
    M.list = M.list,
    estimation = estimation,
    Pmodel.lists = Pmodel.lists,  # list
    Omodel.lists = Omodel.lists,  # list
    Imodel.lists = Imodel.lists,  # list
    Imodel.source = Imodel.source
  )

  # Report successful validation.
  message("Input checks passed successfully.")

  if (estimation != "EIF") {
    message("Please use bootstrapping to estimate standard errors in the mean or effect estimation functions.")
  }


  class(output) <- "pathsInfo"
  return(output)

}
