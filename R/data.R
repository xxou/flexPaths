#' Simulated Data with a Single Binary Treatment
#'
#' Data simulated for the single-treatment examples in flexPaths.
#' Each row represents one simulated individual, with baseline covariates,
#' three ordered mediator stages, and continuous and binary outcomes.
#'
#' @format A data frame with 1000 rows and 9 variables:
#' \describe{
#'   \item{Treat}{Binary treatment, stored as integer 0 or 1.}
#'   \item{X1, X2}{Continuous baseline covariates.}
#'   \item{Med1}{Continuous mediator.}
#'   \item{Med2_1, Med2_2}{The second mediator has two continuous components.}
#'   \item{Med3}{Continuous mediator.}
#'   \item{Outcome1}{Continuous outcome.}
#'   \item{Outcome2}{Binary outcome, stored as integer 0 or 1.}
#' }
#' @details
#' The data use wide format and contain no missing values.
#' @source Simulated by Xiaxian Ou for flexPaths.
#' @seealso [pathsInfo()], [pathsEffect()], [flexMean()]
#' @examples
#' data("singTreat")
#' head(singTreat)
#' @keywords datasets
"singTreat"


#' Simulated Data with Multiple Binary Treatments
#'
#' Data simulated for the multiple-treatment examples in flexPaths.
#' Each row represents one simulated individual, with three binary treatments,
#' six ordered mediator stages, a baseline covariate, and a continuous outcome.
#'
#' @format A data frame with 1000 rows and 11 variables:
#' \describe{
#'   \item{X}{Continuous baseline covariate.}
#'   \item{T1, T2, T3}{Binary treatment variables, stored as integer 0 or 1.}
#'   \item{M1, M2, M3, M4, M5, M6}{Causally ordered continuous mediators.}
#'   \item{Y}{Continuous outcome.}
#' }
#' @details
#' The data use wide format and contain no missing values.
#' @source Simulated by Xiaxian Ou for flexPaths.
#' @seealso [pathsInfo()], [mflexMean()], [flexEffect()]
#' @examples
#' data("multiTreat")
#' head(multiTreat)
#' @keywords datasets
"multiTreat"


#' Medical Expenditure Panel Survey 2004 Sample
#'
#' A sample from the 2004 Medical Expenditure Panel Survey (MEPS),
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
#' @keywords datasets
"meps"
