# Small data and standard models shared by the tests.
test_data <- function(n = 120) {
  withr::local_seed(2026)
  data <- data.frame(X = rnorm(n))
  data$A <- rbinom(n, 1, plogis(0.2 * data$X))
  data$M1 <- 0.3 * data$A + 0.2 * data$X + rnorm(n)
  data$A2 <- rbinom(n, 1, plogis(0.2 * data$A + 0.2 * data$M1))
  data$M2 <- 0.3 * data$A + 0.3 * data$A2 + 0.2 * data$M1 + rnorm(n)
  data$Y <- 3 + 0.5 * data$A + 0.5 * data$A2 + 0.4 * data$M1 +
    0.3 * data$M2 + 0.2 * data$X + rnorm(n)
  data$B <- rbinom(n, 1, plogis(-0.3 + 0.3 * data$A + 0.2 * data$M1))
  data$F <- factor(rep(c("low", "high"), length.out = n))
  data
}

balanced_data <- function(multiple = FALSE) {
  data <- expand.grid(A = 0:1, X = c(-1, 1), U1 = c(-1, 1),
    U2 = c(-1, 1), E = c(-1, 1))
  if (multiple) data <- rbind(transform(data, A2 = 0), transform(data, A2 = 1))
  data$M1 <- data$A + data$X + data$U1
  data$M2 <- 2 * data$A + 0.5 * data$M1 + data$X + data$U2
  if (multiple) data$M2 <- data$M2 + data$A2
  data$Y <- 3 + 2 * data$A + data$M1 + 2 * data$M2 + data$X + data$E
  if (multiple) data$Y <- data$Y + 3 * data$A2
  data
}

test_info <- function(data = test_data(), A = "A", Y = "Y", estimation = "EIF",
                      M.list = list(M1 = "M1", M2 = "M2"), model.outcome = ~ lm()) {
  suppressMessages(pathsInfo(A, Y, "X", data, M.list, estimation,
    model.treatment = if (estimation != "G") ~ glm(family = binomial()),
    model.outcome = if (estimation != "IPW") model.outcome))
}

# Keep known failures visible; opt in only when reproducing or verifying a fix.
skip_known_issue <- function(reason) {
  if (Sys.getenv("FLEXPATHS_RUN_KNOWN_ISSUES") != "true") skip(reason)
}
