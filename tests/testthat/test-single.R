test_that("G agrees with analytic means for all two-mediator paths", {
  fit <- test_info(balanced_data(), estimation = "G")
  active <- expand.grid(M1 = 0:1, M2 = 0:1, Y = 0:1)
  for (i in seq_len(nrow(active))) {
    a <- unlist(active[i, ], use.names = FALSE)
    out <- flexPotential(fit, a)
    expect_s3_class(out, "flexPotential")
    expect_equal(out$results$value, 3 + 2 * a[1] + 4 * a[2] + 2 * a[3], tolerance = 1e-8)
  }
})

test_that("single-treatment scores agree with independent regressions and weights", {
  data <- test_data()
  active <- list(c(0, 0, 0), c(1, 1, 1), c(1, 0, 1), c(0, 1, 0))
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, estimation = estimation)
    for (a in active) {
      out <- flexPotential(fit, a)
      ref <- reference_single(data, a, estimation)
      expect_true(all(is.finite(out$potential_data)))
      expect_equal(as.numeric(out$potential_data), ref, tolerance = 1e-8)
      expect_equal(out$results$value, mean(ref), tolerance = 1e-8)
      expect_equal(out$results$SE, sd(ref) / sqrt(nrow(data)), tolerance = 1e-8)
    }
  }
})
