treatment_data <- function(multiple = FALSE) {
  data <- expand.grid(A = 0:1, X = c(-1, 1), U1 = c(-1, 1),
    U2 = c(-1, 1), E = c(-1, 1))
  if (multiple) {
    data <- rbind(transform(data, A2 = 0), transform(data, A2 = 1))
  }
  data$M1 <- 0.25 * data$A + data$X + data$U1
  data$M2 <- 0.5 * data$A + 0.25 * data$M1 + data$X + data$U2
  if (multiple) data$M2 <- data$M2 + 0.5 * data$A2
  data$Y <- 3 + 2 * data$A + data$M1 + 2 * data$M2 + data$X + data$E
  if (multiple) data$Y <- data$Y + 3 * data$A2
  data
}

treatment_info <- function(data, A = "A", estimation = "G") {
  info <- suppressMessages(pathsInfo(
    A = A, Y = "Y", cov_x = "X", data = data,
    M.list = list(M1 = "M1", M2 = "M2"), estimation = estimation,
    model.treatment = if (estimation != "G") ~ glm(family = binomial()),
    model.outcome = if (estimation != "IPW") ~ lm()
  ))
  info
}

test_that("numeric and integer 0/1 inputs retain their values", {
  for (values in list(rep(c(0, 1), 16), rep(0:1, 16))) {
    data <- treatment_data()
    data$A <- values
    original <- data
    info <- treatment_info(data)
    expect_identical(info$data$A, rep(c(0, 1), 16))
    expect_identical(info$data[names(data) != "A"], data[names(data) != "A"])
    expect_identical(data, original)
  }
})

test_that("binary factors use declared level order, including ordered factors", {
  cases <- list(
    factor(rep(c("control", "treated"), 16)),
    factor(rep(c(0, 2), 16)),
    factor(rep(c("low", "high"), 16), levels = c("low", "high")),
    ordered(rep(c("low", "high"), 16), levels = c("low", "high"))
  )
  for (values in cases) {
    data <- treatment_data()
    data$A <- values
    original <- data
    info <- treatment_info(data)
    expect_identical(info$data$A, rep(c(0, 1), 16))
    expect_identical(info$data[names(data) != "A"], data[names(data) != "A"])
    expect_identical(data, original)
  }
  data$A <- factor(rep(c("low", "high"), 16), levels = c("high", "low"))
  expect_identical(treatment_info(data)$data$A, rep(c(1, 0), 16))
})

test_that("all treatment columns are normalized independently", {
  data <- treatment_data(multiple = TRUE)
  data$A <- factor(data$A, levels = 0:1, labels = c("control", "treated"))
  data$A2 <- as.integer(data$A2)
  info <- treatment_info(data, A = c("A", "A2"))
  expect_identical(info$data$A, rep(c(0, 1), 32))
  expect_identical(info$data$A2, rep(c(0, 1), each = 32))

  data$A2 <- ordered(data$A2, levels = c(1, 0))
  info <- treatment_info(data, A = c("A", "A2"))
  expect_identical(info$data$A, rep(c(0, 1), 32))
  expect_identical(info$data$A2, rep(c(1, 0), each = 32))
})

test_that("invalid numeric treatment values are rejected", {
  for (value in c(2, -1, 0.5, Inf, -Inf)) {
    data <- treatment_data()
    data$A[1] <- value
    expect_error(treatment_info(data), "Numeric treatment 'A' must contain only 0 and 1", fixed = TRUE)
  }
  data <- treatment_data(multiple = TRUE)
  data$A2 <- 2 * data$A2
  expect_error(treatment_info(data, A = c("A", "A2")), "Numeric treatment 'A2'", fixed = TRUE)
})

test_that("other treatment types and nonbinary factors are rejected", {
  data <- treatment_data()
  for (values in list(as.character(data$A), as.logical(data$A),
    as.complex(data$A), as.list(data$A), matrix(data$A, ncol = 1))) {
    data$A <- values
    expect_error(treatment_info(data), "Treatment 'A' must be numeric 0/1 or a factor", fixed = TRUE)
  }
  for (values in list(factor(rep("only", 32)),
    factor(rep(c("a", "b"), 16), levels = c("a", "b", "c")),
    factor(rep(NA_character_, 32)),
    factor(rep(c("a", NA), 16), exclude = NULL))) {
    data$A <- values
    expect_error(treatment_info(data), "Treatment 'A' must be a factor with exactly two non-missing levels", fixed = TRUE)
  }
})

test_that("missing rows are removed without changing factor encoding", {
  for (use_factor in c(FALSE, TRUE)) {
    data <- treatment_data()
    if (use_factor) data$A <- factor(data$A, levels = 0:1, labels = c("low", "high"))
    data$A[1] <- NA
    data$X[2] <- NA
    data$M1[3] <- NA
    original <- data
    info <- treatment_info(data)
    expect_identical(info$data$A, rep(c(0, 1), 16)[-(1:3)])
    expect_identical(info$data[names(data) != "A"], data[-(1:3), names(data) != "A"])
    expect_identical(data, original)
  }
  data <- treatment_data()
  data$A <- factor(data$A, levels = 0:1, labels = c("low", "high"))
  data$X[data$A == "low"] <- NA
  # The remaining level must stay 1 even if no observations of level 0 remain.
  expect_identical(treatment_info(data)$data$A, rep(1, 16))
  data <- data[!is.na(data$X), ]
  expect_identical(treatment_info(data)$data$A, rep(1, 16))
})

test_that("invalid treatments cannot be hidden by missing values elsewhere", {
  data <- treatment_data()
  data$A[1] <- 2
  data$X[1] <- NA
  expect_error(treatment_info(data), "Numeric treatment 'A'", fixed = TRUE)
})

test_that("factor encoding preserves single-treatment estimates and effects", {
  numeric_data <- treatment_data()
  factor_data <- numeric_data
  factor_data$A <- factor(numeric_data$A, levels = 0:1, labels = c("low", "high"))
  for (estimation in c("G", "IPW", "EIF")) {
    numeric_info <- treatment_info(numeric_data, estimation = estimation)
    factor_info <- treatment_info(factor_data, estimation = estimation)
    numeric_p0 <- flexPotential(numeric_info, c(0, 0, 0))
    numeric_p1 <- flexPotential(numeric_info, c(1, 0, 1))
    factor_p0 <- flexPotential(factor_info, c(0, 0, 0))
    factor_p1 <- flexPotential(factor_info, c(1, 0, 1))
    expect_true(all(is.finite(factor_p1$potential_data)))
    expect_equal(factor_p0$potential_data, numeric_p0$potential_data, tolerance = 1e-12)
    expect_equal(factor_p1$potential_data, numeric_p1$potential_data, tolerance = 1e-12)
    expect_equal(flexEffect(factor_p1, factor_p0)$results,
      flexEffect(numeric_p1, numeric_p0)$results, tolerance = 1e-12)
    expect_equal(pathsEffect(factor_info)$results,
      pathsEffect(numeric_info)$results, tolerance = 1e-12)
    if (estimation == "G") {
      # From the additive data-generating equations: 3 + 0.375*a_M1 + a_M2 + 2*a_Y.
      expect_equal(factor_p0$results$value, 3, tolerance = 1e-12)
      expect_equal(factor_p1$results$value, 5.375, tolerance = 1e-12)
    }
  }
})

test_that("factor encoding preserves multiple-treatment estimates and effects", {
  numeric_data <- treatment_data(multiple = TRUE)
  factor_data <- numeric_data
  factor_data$A <- factor(numeric_data$A, levels = 0:1, labels = c("low", "high"))
  factor_data$A2 <- ordered(numeric_data$A2, levels = 0:1, labels = c("low", "high"))
  for (estimation in c("G", "IPW", "EIF")) {
    numeric_info <- treatment_info(numeric_data, A = c("A", "A2"), estimation = estimation)
    factor_info <- treatment_info(factor_data, A = c("A", "A2"), estimation = estimation)
    active0 <- list(c(0, 0, 0), c(NA, 0, 0))
    active1 <- list(c(1, 0, 1), c(NA, 1, 1))
    numeric_p0 <- mflexPotential(numeric_info, active0)
    numeric_p1 <- mflexPotential(numeric_info, active1)
    factor_p0 <- mflexPotential(factor_info, active0)
    factor_p1 <- mflexPotential(factor_info, active1)
    expect_true(all(is.finite(factor_p1$potential_data)))
    expect_equal(factor_p0$potential_data, numeric_p0$potential_data, tolerance = 1e-12)
    expect_equal(factor_p1$potential_data, numeric_p1$potential_data, tolerance = 1e-12)
    expect_equal(flexEffect(factor_p1, factor_p0)$results,
      flexEffect(numeric_p1, numeric_p0)$results, tolerance = 1e-12)
    if (estimation == "G") {
      # A2 contributes 1 through M2 and 3 directly to Y.
      expect_equal(factor_p0$results$value, 3, tolerance = 1e-12)
      expect_equal(factor_p1$results$value, 9.375, tolerance = 1e-12)
    }
  }
})
