test_that("G and IPW setup reminds users to bootstrap", {
  reminder <- paste("Please use bootstrapping to estimate standard errors",
                    "in the mean or effect estimation functions.")
  for (estimation in c("G", "IPW", "EIF")) {
    messages <- capture.output(info <- pathsInfo("A", "Y", "X", test_data(),
      list("M1", "M2"), estimation,
      model.treatment = if (estimation != "G") ~glm(family = binomial()),
      model.outcome = if (estimation != "IPW") ~lm()), type = "message")
    expect_identical(reminder %in% messages, estimation != "EIF")
  }
})

test_that("direct effect inference is limited to EIF differences", {
  data <- test_data()
  data$Y <- 0.3 + 0.01 * data$Y
  for (estimation in c("G", "IPW", "EIF")) {
    for (A in list("A", c("A", "A2"))) {
      fit <- test_info(data, A = A, estimation = estimation)
      if (length(A) == 1) {
        p0 <- flexPotential(fit, c(0, 0, 0))
        p1 <- flexPotential(fit, c(1, 1, 1))
      } else {
        p0 <- mflexPotential(fit, list(c(0, 0, 0), c(NA, 0, 0)))
        p1 <- mflexPotential(fit, list(c(1, 1, 1), c(NA, 1, 1)))
      }
      for (scale in c("diff", "risk", "oddsratio")) {
        inference <- if (estimation == "EIF" && scale == "diff")
          c("SE", "CI.lower", "CI.upper", "P.value")
        out <- flexEffect(p1, p0, scale = scale)
        expect_named(out$results, c("active", "Effect", inference))
        ref <- reference_effect(p1$potential_data, p0$potential_data, scale,
                                estimation = estimation)
        expect_equal(unlist(out$results[names(ref)], use.names = FALSE),
                     unname(ref), tolerance = 1e-8)
        if (length(A) == 1) for (decomposition in c("refer0", "sequential")) {
          out <- pathsEffect(fit, decomposition = decomposition, scale = scale)
          expect_named(out$results, c("Path", "Effect", inference))
        }
      }
    }
  }
})

test_that("G total effects omit the degenerate naive SE and retain bootstrap inference", {
  data <- test_data()
  fit <- test_info(data, estimation = "G")
  p0 <- flexPotential(fit, c(0, 0, 0))
  p1 <- flexPotential(fit, c(1, 1, 1))
  set.seed(704)
  out <- flexEffect(p1, p0, nboot = 10)
  set.seed(704)
  ref <- replicate(10, {
    boot <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
    unname(coef(lm(Y ~ A + X, data = boot))["A"])
  })
  expect_named(out$results, c("active", "Effect"))
  expect_equal(out$results$Effect, unname(coef(lm(Y ~ A + X, data))["A"]))
  expect_gt(out$boot_results$boot.SE, 0)
  expect_equal(out$boot_results$boot.SE, sd(ref), tolerance = 1e-8)
})
