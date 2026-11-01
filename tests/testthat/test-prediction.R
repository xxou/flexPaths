test_that("SuperLearner aligns new data with screened training predictors", {
  withr::local_seed(15)
  data <- test_data()
  data$id <- seq_len(nrow(data))
  newdata <- data[1:7, c("X", "M1", "M2")]
  for (X in list(c("X", "M1", "M2"), c("id", "X", "M1", "M2"))) {
    call <- quote(SuperLearner(SL.library = list(c("SL.glm", "screen.corRank")),
      family = gaussian(), cvControl = list(V = 2)))
    if ("id" %in% X) call$id <- "id"
    fit <- fl_SuperLearner(data, X, "Y", call)
    ref <- predict(fit$model, newdata, onlySL = TRUE)$pred
    expect_equal(predict(fit, newdata), ref, tolerance = 1e-8)
    expect_equal(predict(fit, newdata[, c("M2", "M1", "X")]), ref, tolerance = 1e-8)
    expect_equal(predict(fit, data[1:7, c("id", names(data)[names(data) != "id"])]),
      ref, tolerance = 1e-8)
    expect_equal(predict(fit, fit$dataX[1:7, , drop = FALSE]), ref, tolerance = 1e-8)
    expect_equal(predict(fit, newdata[1, , drop = FALSE]), ref[1, , drop = FALSE], tolerance = 1e-8)
    expect_equal(predict(fit), fit$model$SL.predict, tolerance = 1e-8)
    expect_error(predict(fit, newdata[, "X", drop = FALSE]), "Missing predictors")
  }
})

test_that("SuperLearner supports id with one training predictor", {
  withr::local_seed(16)
  data <- test_data()
  data$id <- seq_len(nrow(data))
  for (Y in c("Y", "B")) {
    call <- quote(SuperLearner(id = "id", SL.library = "SL.mean",
      family = gaussian(), cvControl = list(V = 2)))
    if (Y == "B") call$family <- quote(binomial())
    fit <- fl_SuperLearner(data, c("id", "X"), Y, call)
    expect_equal(as.numeric(predict(fit, data[1:7, ])), rep(mean(data[[Y]]), 7), tolerance = 1e-8)
  }
})

test_that("single-treatment G predictions match analytical values", {
  data <- balanced_data()
  fit <- test_info(data, estimation = "G")
  active <- expand.grid(a1 = 0:1, a2 = 0:1, ay = 0:1)
  for (i in seq_len(nrow(active))) {
    a <- as.numeric(active[i, ])
    ref <- 3 + 2 * a[1] + 4 * a[2] + 2 * a[3] + 5 * data$X
    expect_equal(flexPotential(fit, a)$potential_data, ref, tolerance = 1e-8)
    expect_equal(mflexPotential(fit, list(a))$potential_data, ref, tolerance = 1e-8)
  }
  for (decomposition in c("refer0", "sequential"))
    expect_equal(pathsEffect(fit, decomposition = decomposition)$results$Effect,
      c(2, 4, 2, 8), tolerance = 1e-8)
})

test_that("multiple-treatment G predictions match analytical values", {
  data <- balanced_data(multiple = TRUE)
  fit <- test_info(data, A = c("A", "A2"), estimation = "G")
  for (a1 in 0:1) for (a2 in 0:1) {
    out <- mflexPotential(fit, list(rep(a1, 3), c(NA, a2, a2)))
    expect_equal(out$potential_data, 3 + 8 * a1 + 5 * a2 + 5 * data$X, tolerance = 1e-8)
  }
})

test_that("the installed template can be sourced as a working user model", {
  local_custom_model(
    file = system.file("templates", "fl_model_template.R", package = "flexPaths", mustWork = TRUE),
    models = c("fl_template", "predict.fl_template"))
  data <- test_data()
  newdata <- data[1:7, ]; newdata$X <- newdata$X + 0.3
  model <- fl_template(data, c("X", "M1"), "Y", quote(template()))
  expect_equal(predict(model, newdata), predict(lm(Y ~ X + M1, data), newdata), tolerance = 1e-8)
  fit <- test_info(data, estimation = "G", model.outcome = ~template())
  expect_equal(flexPotential(fit, c(1, 0, 1))$potential_data,
    reference_single(data, c(1, 0, 1), "G"), tolerance = 1e-8)
  expect_equal(pathsEffect(fit)$results$Effect,
    reference_paths(data, "G", "refer0")[, "Effect"], tolerance = 1e-8)
})
