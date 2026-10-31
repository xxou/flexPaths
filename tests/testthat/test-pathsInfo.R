test_that("pathsInfo retains variables and accepts formulas or model lists", {
  data <- test_data()
  info <- test_info(data)
  listed <- suppressMessages(pathsInfo("A", "Y", "X", data,
    M.list = list(M2 = "M2", M1 = "M1"),
    model.treatment = list(~glm(family = binomial())), model.outcome = list(~lm())))
  expect_s3_class(info, "pathsInfo")
  expect_identical(info$data, transform(data, A = as.numeric(A)))
  expect_identical(info$A, "A")
  expect_identical(info$Y, "Y")
  expect_identical(info$cov_x, "X")
  expect_identical(info$M.list, listed$M.list)
  expect_equal(lapply(info$Omodel.lists, deparse), lapply(listed$Omodel.lists, deparse))
  expect_length(info$Pmodel.lists, 3)
  expect_length(info$Imodel.lists, 2)
  expect_identical(info$Imodel.source, "default")

  for (A in list("A", c("A", "A2"))) {
    for (M in list(list(M1 = "M1"), list(M1 = c("M1", "M2")))) {
      fit <- test_info(data, A = A, M.list = M)
      expect_identical(fit$A, A)
      expect_identical(fit$M.list, M)
      expect_true(all(is.finite(if (length(A) == 1)
        flexPotential(fit, c(0, 1))$potential_data else
        mflexPotential(fit, list(c(0, 1), c(NA, 1)))$potential_data)))
    }
  }
})

test_that("per-stage and explicit iterative models are retained", {
  data <- test_data()
  fit <- suppressMessages(pathsInfo("A", "Y", "X", data, list(M1 = "M1", M2 = "M2"),
    model.treatment = list(M2 = ~glm(family = binomial()), cov_x = ~glm(family = binomial()),
                           M1 = ~glm(family = binomial())),
    model.outcome = list(M2 = ~lm(), cov_x = ~lm(), M1 = ~lm()),
    model.iter = list(M2 = ~lm(), M1 = ~glm())))
  expect_named(fit$Omodel.lists, c("cov_x", "M1", "M2"))
  expect_named(fit$Imodel.lists, c("M1", "M2"))
  expect_identical(fit$Imodel.source, "user input")
  expect_equal(flexPotential(fit, c(1, 0, 1))$potential_data,
               flexPotential(test_info(data), c(1, 0, 1))$potential_data, tolerance = 1e-8)
  fit <- suppressMessages(pathsInfo(c("A", "A2"), "Y", "X", data,
    list("M1", "M2"), model.treatment = ~glm(), model.outcome = ~lm(), model.iter = ~lm()))
  expect_length(fit$Imodel.lists, 1)
  expect_identical(fit$type, "multiple treatments")
})

test_that("invalid variables, estimation and model counts are rejected", {
  data <- test_data()
  expect_error(test_info(data, A = "absent"), "must exist")
  expect_error(test_info(data, Y = "absent"), "must exist")
  expect_error(test_info(data, Y = "F"), "must be numeric")
  expect_error(test_info(data, M.list = list("absent")), "must exist")
  expect_error(test_info(data, estimation = "unknown"), "must be one of")
  expect_error(pathsInfo("A", "Y", "absent", data), "must exist")
  expect_error(test_info(data, model.outcome = list(~lm(), ~lm())), "model.outcome")
  expect_error(test_info(data, A = c("A", "A2"), model.outcome = list(~lm(), ~lm())), "model.outcome")
  expect_error(suppressMessages(pathsInfo("A", "Y", "X", data, list("M1", "M2"),
    estimation = "G", model.outcome = ~lm(), model.iter = list(~lm(), ~lm(), ~lm()))), "model.iter")
  expect_error(suppressMessages(pathsInfo("A", "Y", "X", data, list("M1", "M2"),
    estimation = "IPW", model.treatment = list(~glm(), ~glm()))), "model.treatment")
  fit <- test_info(data)
  expect_error(flexPotential(fit, c(0, 1)), "length of 'active'", fixed = TRUE)
})

test_that("mediators and per-stage models retain their numeric time order", {
  data <- test_data()
  M <- paste0("M", 1:12)
  for (i in 3:12) data[[M[i]]] <- data$M2 + i
  mediators <- setNames(as.list(M), M)
  for (k in c(1, 2, 9, 10, 12)) {
    expected <- mediators[seq_len(k)]
    fit <- test_info(data, M.list = rev(expected))
    expect_identical(fit$M.list, expected)
  }

  models <- setNames(lapply(0:12, function(i)
    as.formula(paste0("~glm(control = glm.control(maxit = ", i + 25, "))"))),
    c("cov_x", M))
  fit <- suppressMessages(pathsInfo("A", "Y", "X", data, rev(mediators),
    model.treatment = rev(models), model.outcome = rev(models),
    model.iter = rev(models[-1])))
  expect_identical(fit$M.list, mediators)
  expect_identical(fit$Pmodel.lists, models)
  expect_identical(fit$Omodel.lists, models)
  expect_identical(fit$Imodel.lists, models[-1])

  unnamed <- unname(rev(mediators))
  expect_identical(test_info(data, M.list = unnamed)$M.list, unnamed)
})

test_that("single-treatment active values must be numeric 0/1", {
  fit <- test_info()
  for (value in c(2, -1, 0.5, NA_real_, NaN, Inf))
    expect_error(flexPotential(fit, c(0, value, 1)), "only 0 and 1")
  for (a in list(c("0", "1", "1"), factor(c(0, 1, 1)), c(FALSE, TRUE, TRUE),
                 as.complex(c(0, 1, 1)), matrix(c(0, 1, 1), nrow = 1), list(0, 1, 1)))
    expect_error(flexPotential(fit, a), "numeric vector")
  for (a in list(numeric(), c(0, 1), c(0, 1, 0, 1)))
    expect_error(flexPotential(fit, a), "length of 'active'", fixed = TRUE)
  expect_equal(flexPotential(fit, c(0L, 1L, 1L))$potential_data,
               flexPotential(fit, c(0, 1, 1))$potential_data, tolerance = 1e-8)
})

test_that("multiple-treatment active inputs have valid values and dimensions", {
  fit <- test_info(A = c("A", "A2"))
  a1 <- c(0, 1, 1); a2 <- c(NA, 0, 1)
  for (a in list(a1, rbind(a1, a2), data.frame(a1, a2)))
    expect_error(mflexPotential(fit, a), "must be a list")
  for (a in list(list(), list(a1), list(a1, a2, a2)))
    expect_error(mflexPotential(fit, a), "number of treatments")
  for (a in list(list(c(0, 1), c(NA, 1)), list(a1, c(NA, 1)),
                 list(c(0, 1, 1, 1), a2), list(numeric(), a2)))
    expect_error(mflexPotential(fit, a), "length of 'M.list' + 1", fixed = TRUE)
  for (value in c(2, -1, 0.5, NaN, Inf))
    expect_error(mflexPotential(fit, list(a1, c(NA, value, 1))), "only 0, 1 or NA")
  for (a in list(c(NA, "0", "1"), factor(a2), c(NA, FALSE, TRUE),
                 as.complex(a2), matrix(a2, nrow = 1), list(NA, 0, 1)))
    expect_error(mflexPotential(fit, list(a1, a)), "numeric vector")
  out <- mflexPotential(fit, list(c(0L, 1L, 1L), c(NA_integer_, 0L, 1L)))
  expect_true(is.na(out$active[2, 1]))
  expect_equal(out$potential_data, mflexPotential(fit, list(a1, a2))$potential_data,
               tolerance = 1e-8)
})
