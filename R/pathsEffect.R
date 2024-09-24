#' Causal Path-Specific Effect Estimation for Single Treatment
#'
#' This function estimates causal path-specific effects (PSE) for a single treatment,
#' focusing on specific pathways through which the treatment influences the outcome.
#' The estimation is based on the fitted models and data provided by the \code{pathsFit} function.
#'
#' **Decomposition Methods**:
#' Path-specific effects can be decomposed in various ways to reveal how a treatment's effect is mediated through different pathways.
#' This function supports two primary methods of decomposition:
#'
#' 1. **Cumulative Path-Specific Effects (CPSE)**:
#'    - This method captures the effect of a treatment by sequentially activating all pathways that involve a particular mediator while blocking others.
#'    - These CPSEs can be summed to yield the total effect of the treatment, allowing the decomposition of the total effect into path-specific components. This method is particularly useful for calculating the proportion of mediation.
#'
#' 2. **Natural Path-Specific Effects (NPSE)**:
#'    - NPSEs focus on the specific pathways of interest by treating the treatment variable as set to an 'active' value (A = 1) in some pathways while keeping it at a 'baseline' value (A = 0) in others. This approach helps in isolating the main effect of each mediator.
#'    - For example, NPSEs can be interpreted as the average change in the outcome Y if a marginalized group had the same mediator values as a privileged group.
#'    - Unlike CPSEs, NPSEs cannot be simply summed to represent the total effect in the presence of interactions between the treatment and mediators.
#'
#' @param pathsFit An object of class \code{pathsFit}, which is the output from the \code{pathsFit} function, containing the necessary fitted models and data.
#' @param decomposition A character string specifying the method of decomposition to use for effect calculation. Options are "refer0" (aligned with CPSE) and "sequential" (aligned with NPSE). The default is "refer0".
#' @param scale A character string specifying the scale of the effect. Options are \code{"diff"}, \code{"risk"}, and \code{"oddsratio"}. Only the "diff" scale will provide standard errors and p-values directly from the estimated data. To obtain these for the other two scales, bootstrapping should be used.
#' @param CI_level A numeric value indicating the confidence level for confidence intervals, typically set between 0 and 1. The default is 0.95, which corresponds to 95% confidence intervals.
#' @param nboot An optional numeric value specifying the number of bootstrap samples to draw for estimating the standard errors and confidence intervals of the effects. If NULL, no bootstrapping is performed.
#' @param m.cores An optional numeric value indicating the number of cores to use for parallel processing during bootstrapping. If NULL, the function will use a single core.
#'
#' @return An object of class \code{pathsEffect}, containing the following components:
#' \describe{
#'   \item{\code{call}}{The matched call to the \code{pathsEffect} function.}
#'   \item{\code{pathsFit}}{The original \code{pathsFit} object that was used as input.}
#'   \item{\code{decomposition}}{A character string indicating the method of decomposition used ("refer0" or "sequential decomposition").}
#'   \item{\code{estimation}}{The type of estimation method used, such as "EIF", "IPW", or "G-computation".}
#'   \item{\code{scale}}{The scale used for the effect estimation, such as "diff", "risk", or "oddsratio".}
#'   \item{\code{CI_level}}{The confidence level used for the confidence intervals.}
#'   \item{\code{data}}{The original data set that was used for the analysis.}
#'   \item{\code{potential_data}}{The potential outcome data generated during the analysis.}
#'   \item{\code{results}}{A data frame containing the calculated effects for each path, including standard errors, confidence intervals, and p-values for the "diff" scale.}
#'   \item{\code{boot_results}}{If bootstrapping is performed, a data frame containing the bootstrap estimates of the effects, including standard errors, confidence intervals, and p-values.}
#' }
#'
#' @import stats purrr
#' @import SuperLearner dbarts
#' @importFrom dplyr mutate group_by summarise relocate
#' @importFrom parallel mclapply detectCores
#'
#' @export
#'
#' @examples
#' data("singTreat")
#' EIF_fit <- pathsFit(
#'   data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
#'   M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
#'   estimation = "EIF",
#'   model.outcome = list(~ glm(family = gaussian())),
#'   model.propensity = ~ bart(verbose = FALSE, ndpost = 200)
#' )
#' effect_results1 <- pathsEffect(
#'   pathsFit = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95
#' )
#' effect_results2 <- pathsEffect(
#'   pathsFit = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95,
#'   nboot = 100, m.cores = 2
#' )
#' effect_results3 <- pathsEffect(
#'   pathsFit = EIF_fit, decomposition = "sequential", scale = "diff", CI_level = 0.95,
#'   nboot = 100, m.cores = 2
#' )
#'
#'
pathsEffect<-function(
    pathsFit , decomposition = "refer0", scale = "diff",
    CI_level = 0.95 , nboot = NULL, m.cores = NULL){
  # evironment
  Call <- match.call()
  list2env(pathsFit, envir = environment())

  #### PSE decomposition -----
  if(decomposition == "refer0"){
    out<- refer0_effect(data = data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
                        CI_level = CI_level, scale = scale, estimation =estimation,
                        Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists)
  }
  if(decomposition == "sequential"){
    out<- sequential_effect(data = data, A = A, Y = Y, cov_x = cov_x, M.list = M.list,
                            CI_level = CI_level, scale = scale, estimation =estimation,
                            Pmodel.lists = Pmodel.lists, Omodel.lists = Omodel.lists, Imodel.lists = Imodel.lists)
  }


  #### bootstrap ------
  if(is.numeric(nboot)){
     # check_SuperLearner return warning

    if (is.numeric(m.cores)) {
      m.cores.detected <- parallel::detectCores()
      if(m.cores > m.cores.detected){m.cores <- m.cores.detected}

      # windows, m.cores = 1
      os <- tolower(Sys.info()[['sysname']])
      if (os == "windows") {
        m.cores <- 1
        message("Running on Windows. Setting cores to 1 since mclapply is not available.")
      }

      boot.list <- parallel::mclapply(seq_len(nboot), function(i) {
        one_boot.pathsEffect(out)}, mc.cores = m.cores)
    }else{
      boot.list = purrr::map(seq_len(nboot), ~one_boot.pathsEffect(out))
    }

    results.combine <- purrr::map_dfr(boot.list, ~ .x , .id = "source")
    output<- results.combine %>%
      group_by(Path) %>%
      summarise(
        boot.SE = sd(Effect,na.rm = T),
        boot.CI.lower = quantile(Effect, probs = (1-CI_level)/2,na.rm = T),
        boot.CI.upper = quantile(Effect, probs = 1-(1-CI_level)/2,na.rm = T)) %>%
      mutate(Effect = out$results$Effect,
             # p value for risk and ratio???
             nboot = nboot) %>%
      relocate(Effect, .after = Path) %>% as.data.frame()

    if(scale =="diff"){output$boot.P.value = round(2*(1-pnorm(abs(output$Effect)/output$boot.SE)),4)}
    out$boot_results <- output
  }

  out$call <- Call; out$pathsFit = pathsFit
  class(out) <- "pathsEffect"
  return(out)
}
