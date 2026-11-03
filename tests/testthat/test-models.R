test_that("lm and glm wrappers match stats on new data", {
  data <- test_data()
  newdata <- data[1:7, ]; newdata$X <- newdata$X + 0.3
  for (X in list("X", c("X", "F"))) {
    fit <- fl_lm(data, X, "Y", quote(lm()))
    ref <- lm(reformulate(X, "Y"), data)
    expect_s3_class(fit, "fl_lm")
    expect_equal(predict(fit, newdata), predict(ref, newdata), tolerance = 1e-8)
    fit <- fl_glm(data, X, "B", quote(glm(family = binomial())))
    ref <- glm(reformulate(X, "B"), data, family = binomial())
    expect_equal(predict(fit, newdata), predict(ref, newdata, type = "response"), tolerance = 1e-8)
  }
  fit <- fl_lm(data, c("X", "F"), "Y", quote(lm(formula = unused ~ X * F)))
  expect_equal(predict(fit, newdata), predict(lm(Y ~ X * F, data), newdata), tolerance = 1e-8)
  fit <- fl_glm(data, "X", "B", quote(glm(formula = unused ~ X + I(X ^ 2), family = binomial())))
  expect_equal(predict(fit, newdata),
    predict(glm(B ~ X + I(X ^ 2), data, family = binomial()), newdata, type = "response"), tolerance = 1e-8)
  fit <- fl_glm(data, "X", "Y", quote(glm(family = gaussian())))
  expect_equal(predict(fit, newdata), predict(lm(Y ~ X, data), newdata), tolerance = 1e-8)
})

test_that("BART returns finite continuous and binary predictions", {
  withr::local_seed(11)
  data <- test_data()
  for (Y in c("Y", "B")) {
    fit <- fl_bart(data, c("X", "M1"), Y,
      quote(bart(ntree = 5, nskip = 10, ndpost = 20, nthread = 1, verbose = FALSE)))
    pred <- predict(fit, data[1:7, c("X", "M1")])
    expect_s3_class(fit, "fl_bart")
    expect_length(pred, 7)
    expect_true(all(is.finite(pred)))
    if (Y == "B") expect_true(all(pred >= 0 & pred <= 1))
  }
})

test_that("SuperLearner returns finite continuous and binary predictions", {
  withr::local_seed(12)
  data <- test_data()
  for (Y in c("Y", "B")) for (X in list("X", c("X", "F"))) {
    call <- if (Y == "Y")
      quote(SuperLearner(SL.library = c("SL.mean", "SL.glm"), family = gaussian(), cvControl = list(V = 2))) else
      quote(SuperLearner(SL.library = c("SL.mean", "SL.glm"), family = binomial(), cvControl = list(V = 2)))
    fit <- fl_SuperLearner(data, X, Y, call)
    pred <- predict(fit, data[1:7, X, drop = FALSE])
    expect_s3_class(fit, "fl_SuperLearner")
    expect_length(pred, 7)
    expect_true(all(is.finite(pred)))
    if (Y == "B") expect_true(all(pred >= 0 & pred <= 1))
  }
})

test_that("BART and SuperLearner run through potential outcomes and effects", {
  withr::local_seed(13)
  models <- list(
    ~ bart(ntree = 5, nskip = 10, ndpost = 20, nthread = 1, verbose = FALSE),
    ~ SuperLearner(SL.library = c("SL.mean", "SL.glm"), family = gaussian(), cvControl = list(V = 2)))
  for (model in models) {
    fit <- test_info(model.outcome = model)
    p0 <- flexPotential(fit, c(0, 0, 0))
    p1 <- flexPotential(fit, c(1, 0, 1))
    out <- flexEffect(p1, p0)
    expect_s3_class(out, "flexEffect")
    expect_true(all(is.finite(p1$potential_data)))
    expect_true(all(is.finite(out$results$Effect)))
    expect_equal(out$results$Effect, mean(p1$potential_data - p0$potential_data), tolerance = 1e-8)
  }
})

test_that("a sourced user model matches built-in lm", {
  local_custom_model()
  data <- test_data()
  newdata <- data[1:7, ]; newdata$X <- newdata$X + 0.3
  fit <- fl_testlm(data, c("X", "M1"), "Y", quote(testlm()))
  expect_equal(predict(fit, newdata), predict(lm(Y ~ X + M1, data), newdata), tolerance = 1e-8)
  custom <- test_info(data, model.outcome = ~ testlm())
  builtin <- test_info(data)
  p0 <- flexPotential(custom, c(0, 0, 0))
  p1 <- flexPotential(custom, c(1, 0, 1))
  q0 <- flexPotential(builtin, c(0, 0, 0))
  q1 <- flexPotential(builtin, c(1, 0, 1))
  expect_equal(p0$potential_data, q0$potential_data, tolerance = 1e-8)
  expect_equal(p1$potential_data, q1$potential_data, tolerance = 1e-8)
  expect_equal(flexEffect(p1, p0)$results, flexEffect(q1, q0)$results, tolerance = 1e-8)
})

test_that("a sourced user model works in pathsEffect and bootstrap", {
  withr::local_seed(14)
  local_custom_model()
  custom <- test_info(model.outcome = ~ testlm())
  builtin <- test_info()
  for (decomposition in c("refer0", "sequential")) {
    set.seed(14)
    out <- pathsEffect(custom, decomposition = decomposition, nboot = 10)
    set.seed(14)
    ref <- pathsEffect(builtin, decomposition = decomposition, nboot = 10)
    expect_equal(out$results, ref$results, tolerance = 1e-8)
    expect_equal(out$boot_results, ref$boot_results, tolerance = 1e-8)
  }
})
