#' Simulated Data with a Single Binary Treatment
#'
#' Data simulated by Xiaxian Ou for the single-treatment examples in flexPaths.
#' Each row represents one simulated individual, with baseline covariates,
#' three ordered mediator stages, and continuous and binary outcomes.
#'
#' @format A data frame with 1000 rows and 9 variables:
#' \describe{
#'   \item{Treat}{Binary treatment, stored as integer 0 or 1.}
#'   \item{X1, X2}{Continuous baseline covariates.}
#'   \item{Med1}{Continuous mediator at the first stage.}
#'   \item{Med2_1, Med2_2}{Continuous mediators jointly forming the second stage.}
#'   \item{Med3}{Continuous mediator at the third stage.}
#'   \item{Outcome1}{Continuous outcome.}
#'   \item{Outcome2}{Binary outcome, stored as integer 0 or 1.}
#' }
#' @details
#' The data use wide format and contain no missing values. In the examples,
#' the mediator stages are specified as
#' \code{list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"), M3 = "Med3")}.
#' @source Simulated by Xiaxian Ou for flexPaths.
#' @seealso [pathsInfo()], [pathsEffect()], [flexPotential()]
#' @examples
#' data("singTreat")
#' head(singTreat)
#' @keywords datasets
"singTreat"


#' Simulated Data with Multiple Binary Treatments
#'
#' Data simulated by Xiaxian Ou for the multiple-treatment examples in flexPaths.
#' Each row represents one simulated individual, with three binary treatments,
#' six ordered mediator stages, a baseline covariate, and a continuous outcome.
#'
#' @format A data frame with 1000 rows and 11 variables:
#' \describe{
#'   \item{X}{Continuous baseline covariate.}
#'   \item{T1, T2, T3}{Binary treatment variables, stored as integer 0 or 1.}
#'   \item{M1, M2, M3, M4, M5, M6}{Continuous mediators, one per ordered stage.}
#'   \item{Y}{Continuous outcome.}
#' }
#' @details
#' The data use wide format and contain no missing values. Treatments and
#' mediators at different stages are stored in separate columns. The
#' \code{active} argument of \code{mflexPotential()} specifies the treatment
#' values for each mediator stage and the outcome.
#' @source Simulated by Xiaxian Ou for flexPaths.
#' @seealso [pathsInfo()], [mflexPotential()], [flexEffect()]
#' @examples
#' data("multiTreat")
#' head(multiTreat)
#' @keywords datasets
"multiTreat"


#' Medical Expenditure Panel Survey 2004 Sample
#'
#' A person-level sample from the 2004 Medical Expenditure Panel Survey (MEPS),
#' containing demographic information, health measures, insurance coverage,
#' healthcare spending, and service use.
#'
#' @format A data frame with 19386 rows and 32 variables:
#' \describe{
#'   \item{duid}{Identifier for the dwelling unit.}
#'   \item{dupersid}{Unique identifier for the person.}
#'   \item{hieuidx}{Identifier for the health-insurance eligibility unit.}
#'   \item{wtdper}{Person-level sampling weight.}
#'   \item{age}{Age in years.}
#'   \item{female}{Indicator for female sex.}
#'   \item{race_bl}{Indicator for Black race.}
#'   \item{race_oth}{Indicator for race other than White or Black.}
#'   \item{eth_hisp}{Indicator for Hispanic ethnicity.}
#'   \item{famsize}{Annualized size of the responding family.}
#'   \item{ed_hs}{Indicator for high-school education.}
#'   \item{ed_hsplus}{Indicator for some college education.}
#'   \item{ed_col}{Indicator for college education.}
#'   \item{ed_colplus}{Indicator for education beyond college.}
#'   \item{lninc}{Natural logarithm of family income.}
#'   \item{reg_midw}{Indicator for residence in the Midwest.}
#'   \item{reg_south}{Indicator for residence in the South.}
#'   \item{reg_west}{Indicator for residence in the West.}
#'   \item{anylim}{Indicator for any disability.}
#'   \item{mcs12}{SF-12 mental health score.}
#'   \item{pcs12}{SF-12 physical health score.}
#'   \item{ins_mcare}{Indicator for Medicare coverage.}
#'   \item{ins_mcaid}{Indicator for Medicaid coverage.}
#'   \item{ins_unins}{Indicator for having no insurance.}
#'   \item{ins_dent}{Prorated dental-insurance coverage.}
#'   \item{exp_tot}{Total expenditure on medical care.}
#'   \item{exp_dent}{Expenditure on dental care.}
#'   \item{exp_self}{Expenditure paid by the person or family.}
#'   \item{use_disch}{Number of hospital discharges.}
#'   \item{use_los}{Number of nights spent in hospital.}
#'   \item{use_dent}{Number of visits for dental care.}
#'   \item{use_rx}{Number of prescriptions, including refills.}
#' }
#' @details
#' This is the unchanged \code{meps} data object distributed in version 0.1.0
#' of the \pkg{twopartm} package, including its column types and attributes.
#' It has one row per person and no missing values.
#'
#' The columns \code{female}, \code{race_bl}, \code{race_oth}, \code{eth_hisp},
#' \code{ed_hs}, \code{ed_hsplus}, \code{ed_col}, \code{ed_colplus},
#' \code{reg_south}, \code{reg_west}, \code{anylim}, \code{ins_mcare},
#' \code{ins_mcaid}, and \code{ins_unins} are factors with levels \code{"0"}
#' and \code{"1"}. The columns \code{reg_midw} and \code{ins_dent} retain
#' the upstream \code{haven_labelled} numeric class. The column \code{hieuidx}
#' is character; the remaining columns are numeric.
#' @source
#' The \code{meps} dataset in \pkg{twopartm} 0.1.0, by Yajie Duan, Birol Emir,
#' Griffith Bell, and Javier Cabrera:
#' \url{https://CRAN.R-project.org/package=twopartm}.
#' Its documentation identifies the original source as the datasets accompanying
#' \emph{Health Econometrics Using Stata}:
#' \url{https://www.stata-press.com/data/heus.html}.
#' @references
#' Partha Deb, Edward C. Norton, and Willard G. Manning (2017).
#' \emph{Health Econometrics Using Stata}. College Station, TX: Stata Press.
#' @examples
#' data("meps")
#' head(meps[, c("age", "female", "exp_tot")])
#' @keywords datasets
"meps"
