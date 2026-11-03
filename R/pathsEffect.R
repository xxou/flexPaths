#' Causal Path-Specific Effect Estimation for Single Treatment
#'
#' Estimates causal path-specific effects (PSEs) for a single binary treatment
#' with K (>= 1) causally ordered mediators. The function estimates the effect
#' along each mediator-specific pathway, A -> M_k -> ... -> Y, for k = 1, ..., K,
#' as well as the direct effect along A -> Y.
#' The estimation is based on analysis settings and model specifications
#' provided by the \code{pathsInfo} function.
#'
#' @details **Decomposition methods:** Both options compare means defined by
#' different treatment assignments, as specified by \code{assign} in \code{flexMean}.
#' The entries correspond to the mediator stages followed by the outcome.
#' \itemize{
#'   \item \code{"refer0"} gives natural path-specific effects (NPSEs).
#'   Each pathway-specific assignment is compared with the same all-zero
#'   reference. These effects need not sum to the total effect.
#'   \item \code{"sequential"} gives cumulative path-specific effects (CPSEs).
#'   Adjacent assignments are compared as treatment values are changed from
#'   1 to 0, one mediator stage at a time, followed by the outcome.
#'   On the difference scale, these effects sum to the total effect.
#' }
#' The examples below show the assignment contrasts for three mediator stages
#' and illustrate how to express them using \code{flexMean} and \code{flexEffect}.
#'
#' @param pathsInfo An object of class \code{pathsInfo} containing the analysis settings
#' and model specifications provided by the \code{pathsInfo} function.
#' @param decomposition A character string specifying the decomposition method:
#' \code{"refer0"} for natural path-specific effects (NPSEs), or
#' \code{"sequential"} for cumulative path-specific effects (CPSEs).
#' The default is \code{"refer0"}.
#' @param scale A character string specifying the scale of the effect. Options are
#' \code{"diff"}, \code{"risk"}, and \code{"oddsratio"}. Direct standard errors,
#' confidence intervals, and p-values are returned only for EIF estimation with
#' \code{scale = "diff"}. Otherwise, use \code{nboot} for bootstrap standard
#' errors and confidence intervals.
#' @param CI_level A numeric value indicating the confidence level for confidence intervals, typically set between 0 and 1. The default is 0.95, which corresponds to 95% confidence intervals.
#' @param nboot An optional numeric value specifying the number of bootstrap samples to draw for estimating the standard errors and confidence intervals of the effects. If NULL, no bootstrapping is performed.
#' @param m.cores An optional numeric value indicating the number of cores to use for parallel processing during bootstrapping. If NULL, the function will use a single core.
#'
#' @return An object of class \code{pathsEffect}, containing the following components:
#' \describe{
#'   \item{\code{call}}{The matched call to the \code{pathsEffect} function.}
#'   \item{\code{pathsInfo}}{The original \code{pathsInfo} object that was used as input.}
#'   \item{\code{decomposition}}{A character string indicating the method of decomposition used ("refer0" or "sequential decomposition").}
#'   \item{\code{estimation}}{The type of estimation method used, such as "EIF", "IPW", or "G-computation".}
#'   \item{\code{scale}}{The scale used for the effect estimation, such as "diff", "risk", or "oddsratio".}
#'   \item{\code{CI_level}}{The confidence level used for the confidence intervals.}
#'   \item{\code{data}}{The original data set that was used for the analysis.}
#'   \item{\code{potential_data}}{A data frame of observation-level contributions
#'   to the mean estimators, with one column for each treatment assignment used
#'   in the effect contrasts.}
#'   \item{\code{results}}{A data frame containing the estimated effects for each path.
#'   SE, CI.lower, CI.upper, and P.value are included only for EIF estimation
#'   with \code{scale = "diff"}; otherwise, these columns are omitted.}
#'   \item{\code{boot_results}}{If bootstrapping is performed, a data frame containing
#'   the effect estimates, bootstrap standard errors and percentile confidence intervals.
#'   P-values based on a normal approximation using the bootstrap standard errors
#'   are included only for \code{scale = "diff"}.}
#' }
#'
#' @import stats purrr
#' @import SuperLearner dbarts
#' @importFrom dplyr mutate group_by summarise relocate
#' @importFrom parallel mclapply detectCores
#'
#' @export
#'
#' @example inst/examples/pathsEffect-example.R
#'
#'
pathsEffect <- function(
  pathsInfo, decomposition = "refer0", scale = "diff",
  CI_level = 0.95, nboot = NULL, m.cores = NULL) {
  # evironment
  Call <- match.call()
  # Variables populated by list2env().
  data <- type <- A <- Y <- cov_x <- M.list <- estimation <- NULL
  Pmodel.lists <- Omodel.lists <- Imodel.lists <- Imodel.source <- NULL
  list2env(pathsInfo, envir = environment())
  # Column names used in dplyr expressions.
  Path <- Effect <- NULL

  if (type != "single treatment") {
    stop("This function is only used for single treatment")
  }
  #### PSE decomposition -----
  if (decomposition == "refer0") {
    out <- refer0_effect(data = data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
      CI_level = CI_level, scale = scale, estimation = estimation,
      Pmodel.lists = Pmodel.lists,
      Omodel.lists = Omodel.lists,
      Imodel.lists = Imodel.lists,
      Imodel.source = Imodel.source)
  }
  if (decomposition == "sequential") {
    out <- sequential_effect(data = data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
      CI_level = CI_level, scale = scale, estimation = estimation,
      Pmodel.lists = Pmodel.lists,
      Omodel.lists = Omodel.lists,
      Imodel.lists = Imodel.lists,
      Imodel.source = Imodel.source)
  }


  #### bootstrap ------
  if (is.numeric(nboot)) {
     # check_SuperLearner return warning

    if (is.numeric(m.cores)) {
      m.cores.detected <- parallel::detectCores()
      if (m.cores > m.cores.detected) {m.cores <- m.cores.detected}

      # windows, m.cores = 1
      os <- tolower(Sys.info()[['sysname']])
      if (os == "windows") {
        m.cores <- 1
        message("Running on Windows. Setting cores to 1 since mclapply is not available.")
      }

      boot.list <- parallel::mclapply(seq_len(nboot), function(i) {
        one_boot.pathsEffect(out)}, mc.cores = m.cores)
    } else {
      boot.list = purrr::map(seq_len(nboot), ~ one_boot.pathsEffect(out))
    }

    results.combine <- purrr::map_dfr(boot.list, ~.x, .id = "source")
    output <- results.combine %>%
      group_by(Path) %>%
      summarise(
        boot.SE = sd(Effect, na.rm = T),
        boot.CI.lower = quantile(Effect, probs = (1 - CI_level) / 2, na.rm = T),
        boot.CI.upper = quantile(Effect, probs = 1 - (1 - CI_level) / 2, na.rm = T)) %>%
      mutate(Effect = out$results$Effect[match(Path, out$results$Path)],
             # p value for risk and ratio???
        nboot = nboot) %>%
      relocate(Effect, .after = Path) %>% as.data.frame()

    if (scale == "diff") {output$boot.P.value = round(2 * (1 - pnorm(abs(output$Effect) / output$boot.SE)), 4)}
    out$boot_results <- output
  }

  out$call <- Call; out$pathsInfo = pathsInfo
  class(out) <- "pathsEffect"
  return(out)
}
