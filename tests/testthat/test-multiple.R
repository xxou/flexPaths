test_that("multiple-treatment G respects temporal order and analytic means", {
  fit <- test_info(balanced_data(TRUE), A = c("A", "A2"), estimation = "G")
  active <- list(list(c(0, 0, 0), c(NA, 0, 0)),
                 list(c(1, 1, 1), c(NA, 1, 1)),
                 list(c(0, 1, 1), c(NA, 0, 1)))
  for (i in seq_along(active)) {
    out <- mflexPotential(fit, active[[i]])
    expect_s3_class(out, "mflexPotential")
    expect_true(is.na(out$active[2, 1]))
    expect_equal(out$results$value, c(3, 16, 12)[i], tolerance = 1e-8)
  }
})

test_that("fixed two-treatment regimes agree with independent sequential formulas", {
  data <- test_data()
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, A = c("A", "A2"), estimation = estimation)
    for (a1 in 0:1) for (a2 in 0:1) {
      out <- mflexPotential(fit, list(rep(a1, 3), c(NA, a2, a2)))
      ref <- reference_multi(data, a1, a2, estimation)
      expect_equal(as.numeric(out$potential_data), ref, tolerance = 1e-8)
      expect_equal(out$results$value, mean(ref), tolerance = 1e-8)
      if (estimation == "EIF") {
        expect_equal(out$results$SE, sd(ref) / sqrt(nrow(data)), tolerance = 1e-8)
      } else {
        expect_named(out$results, c("active", "value"))
      }
    }
  }
})

test_that("single and multiple implementations agree for one treatment", {
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(estimation = estimation)
    for (a in list(c(0, 0, 0), c(1, 1, 1), c(0, 1, 1), c(1, 0, 1))) {
      single <- flexPotential(fit, a)
      multiple <- mflexPotential(fit, list(a))
      expect_equal(multiple$potential_data, single$potential_data, tolerance = 1e-8)
      expect_equal(multiple$results, single$results, tolerance = 1e-8)
    }
  }
})
