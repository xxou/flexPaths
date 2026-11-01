## Submission

This is the first submission of flexPaths to CRAN.

## Test environments

* Local macOS aarch64, R 4.6.0

## R CMD check results

0 errors | 0 warnings | 4 notes

Checked with `R CMD check --as-cran` on 2026-11-01. All 497 test assertions,
examples and the PDF manual check passed.

After that check, `pathsInfo()` was updated to require at least one covariate.
All 533 local test assertions passed; a full check of this update is pending.

* New submission: this is the first submission of flexPaths to CRAN.
* Three existing notes remain pending: static code checks, Rd example line
  widths, and HTML validation skipped because the local HTML Tidy version is
  too old.
