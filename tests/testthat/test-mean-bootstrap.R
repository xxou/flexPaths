test_that("potential mean bootstrap agrees with independent refits", {
  withr::local_seed(31)
  data <- test_data()
  for (estimation in c("G", "IPW", "EIF")) {
    for (A in list("A", c("A", "A2"))) {
      fit <- test_info(data, A = A, estimation = estimation)
      set.seed(31)
      if (length(A) == 1) {
        out <- flexPotential(fit, c(1, 0, 1), nboot = 10)
        plain <- flexPotential(fit, c(1, 0, 1))
      } else {
        out <- mflexPotential(fit, list(c(1, 1, 1), c(NA, 0, 0)), nboot = 10)
        plain <- mflexPotential(fit, list(c(1, 1, 1), c(NA, 0, 0)))
      }
      set.seed(31)
      ref <- replicate(10, {
        boot <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
        mean(if (length(A) == 1) reference_single(boot, c(1, 0, 1), estimation)
        else reference_multi(boot, 1, 0, estimation))
      })
      expect_equal(out$results, plain$results)
      expect_equal(out$potential_data, plain$potential_data)
      expect_equal(out$pathsInfo, fit)
      expect_null(plain$boot_results)
      expect_named(out$boot_results, c("active", "value", "boot.SE", "nboot"))
      expect_identical(out$boot_results$active, out$results$active)
      expect_equal(out$boot_results$value, out$results$value)
      expect_equal(out$boot_results$boot.SE, sd(ref), tolerance = 1e-8)
      expect_equal(out$boot_results$nboot, 10)
      expect_output(print(out), "boot.SE")
      expect_output(print(out), paste0("estimation: ", estimation))
    }
  }
})

test_that("invalid mean bootstrap counts are rejected", {
  for (model_fun in list(flexPotential, mflexPotential)) {
    for (nboot in list(0, 1, -1, 2.5, NA_real_, NaN, Inf, "10", TRUE,
      c(2, 3), numeric(), 2 + 1i))
      expect_error(model_fun(test_info(), c(1, 0, 1), nboot = nboot), "'nboot'")
  }
})

test_that("mean bootstrap does not nest inside effect bootstrap", {
  withr::local_seed(32)
  for (A in list("A", c("A", "A2"))) {
    fit <- test_info(A = A, estimation = "G")
    if (length(A) == 1) {
      p0 <- flexPotential(fit, c(0, 0, 0), nboot = 3)
      p1 <- flexPotential(fit, c(1, 1, 1), nboot = 3)
      q0 <- flexPotential(fit, c(0, 0, 0))
      q1 <- flexPotential(fit, c(1, 1, 1))
    } else {
      p0 <- mflexPotential(fit, list(c(0, 0, 0), c(NA, 0, 0)), nboot = 3)
      p1 <- mflexPotential(fit, list(c(1, 1, 1), c(NA, 1, 1)), nboot = 3)
      q0 <- mflexPotential(fit, list(c(0, 0, 0), c(NA, 0, 0)))
      q1 <- mflexPotential(fit, list(c(1, 1, 1), c(NA, 1, 1)))
    }
    set.seed(32)
    out <- flexEffect(p1, p0, nboot = 10)
    set.seed(32)
    ref <- flexEffect(q1, q0, nboot = 10)
    expect_equal(out$results, ref$results)
    expect_equal(out$boot_results, ref$boot_results)
  }
})

test_that("a custom model works in potential mean bootstrap", {
  withr::local_seed(33)
  local_custom_model()
  custom <- test_info(model.outcome = ~ testlm())
  builtin <- test_info()
  set.seed(33)
  out <- flexPotential(custom, c(1, 0, 1), nboot = 10)
  set.seed(33)
  ref <- flexPotential(builtin, c(1, 0, 1), nboot = 10)
  expect_equal(out$results, ref$results, tolerance = 1e-8)
  expect_equal(out$boot_results, ref$boot_results, tolerance = 1e-8)
})
