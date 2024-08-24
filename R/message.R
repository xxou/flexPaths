.onLoad <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("purrr", quietly = TRUE)
  })
}
