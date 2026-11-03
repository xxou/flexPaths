test_that("effect differences, SEs and intervals agree with individual scores", {
  fit <- test_info()
  p0 <- flexPotential(fit, c(0, 0, 0))
  p1 <- flexPotential(fit, c(1, 0, 1))
  out <- flexEffect(p1, p0, CI_level = 0.90)
  ref <- reference_effect(p1$potential_data, p0$potential_data, CI_level = 0.90)
  expect_s3_class(out, "flexEffect")
  expect_equal(unlist(out$results[names(ref)], use.names = FALSE), unname(ref), tolerance = 1e-8)
  expect_equal(flexEffect(p1, p1)$results$Effect, 0, tolerance = 1e-8)
})

test_that("ratios use potential means and odds on the probability scale", {
  fit <- test_info(Y = "B", estimation = "G", model.outcome = ~ glm(family = binomial()))
  p0 <- flexPotential(fit, c(0, 0, 0))
  p1 <- flexPotential(fit, c(1, 0, 1))
  expect_true(all(c(p0$results$value, p1$results$value) > 0.1 &
    c(p0$results$value, p1$results$value) < 0.9))
  for (scale in c("risk", "oddsratio")) {
    ref <- reference_effect(p1$potential_data, p0$potential_data, scale)
    expect_equal(flexEffect(p1, p0, scale = scale)$results$Effect, unname(ref), tolerance = 1e-8)
  }
})

test_that("list comparisons retain labels, values and intervals after reordering", {
  fit <- test_info()
  p0 <- flexPotential(fit, c(0, 0, 0))
  p1 <- flexPotential(fit, c(1, 0, 1))
  p2 <- flexPotential(fit, c(0, 1, 0))
  out <- flexEffect(list(p1, p2, p1), p0)$results
  expect_equal(nrow(out), 3)
  for (p in list(p1, p2)) {
    ref <- flexEffect(p, p0)$results
    rows <- out[out$active == ref$active, ]
    for (name in c("Effect", "SE", "CI.lower", "CI.upper"))
      expect_equal(rows[[name]], rep(ref[[name]], nrow(rows)), tolerance = 1e-8)
  }
  expect_equal(flexEffect(list(p1, p2), list(p0, p0))$results,
    flexEffect(list(p2, p1), p0)$results, tolerance = 1e-8)
  expect_equal(flexEffect(p0, list(p1, p2))$results$Effect,
    -flexEffect(list(p2, p1), p0)$results$Effect, tolerance = 1e-8)
})

test_that("both decompositions agree with independently calculated path contrasts", {
  data <- test_data()
  labels <- c("A->M1->...->Y", "A->M2->...->Y", "A->Y", "total effect: A->...->Y")
  for (estimation in c("G", "IPW", "EIF")) {
    fit <- test_info(data, estimation = estimation)
    for (decomposition in c("refer0", "sequential")) {
      out <- pathsEffect(fit, decomposition = decomposition)
      ref <- reference_paths(data, estimation, decomposition)
      expect_s3_class(out, "pathsEffect")
      expect_identical(out$results$Path, labels)
      expect_named(out$results, c("Path", colnames(ref),
        if (estimation == "EIF") "P.value"))
      expect_equal(unname(as.matrix(out$results[, colnames(ref)])), unname(ref), tolerance = 1e-8)
      if (decomposition == "sequential")
        expect_equal(sum(out$results$Effect[1:3]), out$results$Effect[4], tolerance = 1e-8)
    }
  }
})

test_that("decomposition ratios agree with independent means in the unit interval", {
  data <- balanced_data()
  data$Y <- 0.3 + 0.01 * data$Y
  fit <- test_info(data, estimation = "G")
  for (decomposition in c("refer0", "sequential")) for (scale in c("risk", "oddsratio")) {
    out <- pathsEffect(fit, decomposition = decomposition, scale = scale)
    expect_equal(out$results$Effect,
      as.numeric(reference_paths(data, "G", decomposition, scale)), tolerance = 1e-8)
  }
})
