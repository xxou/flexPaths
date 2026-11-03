test_that("G agrees with analytic means for all two-mediator paths", {
  fit <- test_info(balanced_data(), estimation = "G")
  assign <- expand.grid(M1 = 0:1, M2 = 0:1, Y = 0:1)
  for (i in seq_len(nrow(assign))) {
    a <- unlist(assign[i, ], use.names = FALSE)
    out <- flexMean(fit, a)
    expect_s3_class(out, "flexMean")
    expect_equal(out$results$value, 3 + 2 * a[1] + 4 * a[2] + 2 * a[3], tolerance = 1e-8)
  }
})

test_that("single-treatment scores agree with independent regressions and weights", {
  data <- test_data()
  assign <- list(c(0, 0, 0), c(1, 1, 1), c(1, 0, 1), c(0, 1, 0))
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, estimation = estimation)
    for (a in assign) {
      out <- flexMean(fit, a)
      ref <- reference_single(data, a, estimation)
      expect_true(all(is.finite(out$potential_data)))
      expect_equal(as.numeric(out$potential_data), ref, tolerance = 1e-8)
      expect_equal(out$results$value, mean(ref), tolerance = 1e-8)
      if (estimation == "EIF") {
        expect_equal(out$results$SE, sd(ref) / sqrt(nrow(data)), tolerance = 1e-8)
      } else {
        expect_named(out$results, c("assign", "value"))
      }
    }
  }
})
