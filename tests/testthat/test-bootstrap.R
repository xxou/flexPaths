test_that("flexEffect bootstrap agrees with an independent row-resampling loop", {
  withr::local_seed(21)
  data <- test_data()
  data$Y <- 0.3 + 0.01 * data$Y
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, estimation = estimation)
    p0 <- flexMean(fit, c(0, 0, 0))
    p1 <- flexMean(fit, c(1, 0, 1))
    p2 <- flexMean(fit, c(0, 1, 0))
    for (scale in c("diff", "risk", "oddsratio")) {
      set.seed(21)
      fit_boot <- flexEffect(list(p1, p2), p0, scale = scale, nboot = 10)
      if (estimation != "EIF" || scale != "diff")
        expect_named(fit_boot$results, c("assign", "Effect"))
      out <- fit_boot$boot_results
      expect_identical("boot.P.value" %in% names(out), scale == "diff")
      set.seed(21)
      ref <- replicate(10, {
        boot <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
        q0 <- reference_single(boot, c(0, 0, 0), estimation)
        vapply(list(c(1, 0, 1), c(0, 1, 0)), function(a) {
          q1 <- reference_single(boot, a, estimation)
          unname(reference_effect(q1, q0, scale, estimation = estimation)["Effect"])
        }, numeric(1))
      })
      labels <- c("101 vs 000", "010 vs 000")
      for (i in 1:2) {
        row <- out[out$assign == labels[i], ]
        expect_equal(row$boot.SE, sd(ref[i, ]), tolerance = 1e-8)
        expect_equal(c(row$boot.CI.lower, row$boot.CI.upper),
          unname(quantile(ref[i, ], c(0.025, 0.975))), tolerance = 1e-8)
        expect_equal(row$nboot, 10)
      }
    }
  }
})

test_that("pathsEffect bootstrap agrees with independent path resampling", {
  withr::local_seed(22)
  data <- test_data()
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, estimation = estimation)
    for (decomposition in c("refer0", "sequential")) {
      set.seed(22)
      out <- pathsEffect(fit, decomposition = decomposition, nboot = 10)
      if (estimation != "EIF") expect_named(out$results, c("Path", "Effect"))
      set.seed(22)
      ref <- replicate(10, {
        boot <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
        reference_paths(boot, estimation, decomposition)[, "Effect"]
      })
      for (i in 1:4) {
        row <- out$boot_results[out$boot_results$Path == out$results$Path[i], ]
        expect_equal(row$Effect, out$results$Effect[i], tolerance = 1e-8)
        expect_equal(row$boot.SE, sd(ref[i, ]), tolerance = 1e-8)
        expect_equal(c(row$boot.CI.lower, row$boot.CI.upper),
          unname(quantile(ref[i, ], c(0.025, 0.975))), tolerance = 1e-8)
      }
    }
  }
})

test_that("multiple-treatment bootstrap agrees with independent refits", {
  data <- test_data()
  data$Y <- 0.3 + 0.01 * data$Y
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, A = c("A", "A2"), estimation = estimation)
    p0 <- mflexMean(fit, list(c(0, 0, 0), c(NA, 0, 0)))
    p1 <- mflexMean(fit, list(c(1, 1, 1), c(NA, 1, 1)))
    for (scale in c("diff", "risk", "oddsratio")) {
      set.seed(26)
      out <- flexEffect(p1, p0, scale = scale, nboot = 10)
      if (estimation != "EIF" || scale != "diff")
        expect_named(out$results, c("assign", "Effect"))
      set.seed(26)
      ref <- replicate(10, {
        boot <- data[sample(seq_len(nrow(data)), replace = TRUE), ]
        q0 <- reference_multi(boot, 0, 0, estimation)
        q1 <- reference_multi(boot, 1, 1, estimation)
        unname(reference_effect(q1, q0, scale, estimation = estimation)["Effect"])
      })
      expect_equal(out$boot_results$boot.SE, sd(ref), tolerance = 1e-8)
      expect_equal(c(out$boot_results$boot.CI.lower, out$boot_results$boot.CI.upper),
        unname(quantile(ref, c(0.025, 0.975))), tolerance = 1e-8)
      expect_identical("boot.P.value" %in% names(out$boot_results), scale == "diff")
    }
  }
})

test_that("repeated comparisons do not change bootstrap uncertainty", {
  withr::local_seed(23)
  data <- test_data()
  data$Y <- 0.3 + 0.01 * data$Y
  fit <- test_info(data)
  p0 <- flexMean(fit, c(0, 0, 0))
  p1 <- flexMean(fit, c(1, 0, 1))
  p2 <- flexMean(fit, c(0, 1, 0))
  for (scale in c("diff", "risk", "oddsratio")) {
    for (pairs in list(list(p1, p1), list(p1, p2, p1))) {
      set.seed(23)
      ref <- flexEffect(unique(pairs), p0, scale = scale, nboot = 10)$boot_results
      set.seed(23)
      out <- flexEffect(pairs, p0, scale = scale, nboot = 10)$boot_results
      expected <- ref[match(out$assign, ref$assign), ]
      rownames(expected) <- NULL
      expect_equal(nrow(out), length(pairs))
      expect_equal(out, expected, tolerance = 1e-8)
      expect_false("comparison_id" %in% names(out))
    }
  }
})

test_that("a user model also works in the bootstrap refits", {
  withr::local_seed(24)
  local_custom_model()
  custom <- test_info(model.outcome = ~ testlm())
  builtin <- test_info()
  p0 <- flexMean(custom, c(0, 0, 0))
  p1 <- flexMean(custom, c(1, 0, 1))
  q0 <- flexMean(builtin, c(0, 0, 0))
  q1 <- flexMean(builtin, c(1, 0, 1))
  set.seed(24)
  custom_boot <- flexEffect(p1, p0, nboot = 10)$boot_results
  set.seed(24)
  builtin_boot <- flexEffect(q1, q0, nboot = 10)$boot_results
  expect_equal(custom_boot, builtin_boot, tolerance = 1e-8)
})

test_that("bootstrap point estimates stay paired with labels for ten mediators", {
  withr::local_seed(25)
  data <- test_data()
  data$Y <- 0.3 + 0.01 * data$Y
  M <- paste0("M", 1:10)
  for (i in 3:10) data[[M[i]]] <- sin(seq_len(nrow(data)) * i) + 0.1 * i * data$A
  fit <- test_info(data, M.list = as.list(M), estimation = "G")
  for (decomposition in c("refer0", "sequential")) {
    for (scale in c("diff", "risk", "oddsratio")) {
      out <- pathsEffect(fit, decomposition = decomposition, scale = scale, nboot = 10)
      index <- match(out$results$Path, out$boot_results$Path)
      expect_equal(out$boot_results$Effect[index], out$results$Effect, tolerance = 1e-8)
      if (scale == "diff") {
        expected <- round(2 * pnorm(-abs(out$results$Effect) / out$boot_results$boot.SE[index]), 4)
        expect_equal(out$boot_results$boot.P.value[index], expected, tolerance = 1e-8)
      }
    }
  }
})
