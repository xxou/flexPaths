# Independent calculations for two mediators, using only base R regressions.
# These deliberately do not use flexPaths model or estimation helpers.
reference_single <- function(data, assign, estimation) {
  a1 <- assign[1]; a2 <- assign[2]; ay <- assign[3]
  q2 <- predict(lm(Y ~ A + X + M1 + M2, data), transform(data, A = ay))
  q1 <- predict(lm(q2 ~ A + X + M1, data), transform(data, A = a2))
  q0 <- predict(lm(q1 ~ A + X, data), transform(data, A = a1))
  if (estimation == "G") return(as.numeric(q0))

  p0 <- predict(glm(A ~ X, data, family = binomial()), type = "response")
  p1 <- predict(glm(A ~ X + M1, data, family = binomial()), type = "response")
  p2 <- predict(glm(A ~ X + M1 + M2, data, family = binomial()), type = "response")
  prob <- function(p, a) if (a == 1) p else 1 - p
  w0 <- (data$A == a1) / prob(p0, a1)
  w1 <- (data$A == a2) / prob(p0, a2) *
    prob(p1, a1) / prob(p1, a2) * prob(p0, a2) / prob(p0, a1)
  w2 <- (data$A == ay) / prob(p0, ay) *
    prob(p1, a1) / prob(p1, ay) * prob(p0, ay) / prob(p0, a1) *
    prob(p2, a2) / prob(p2, ay) * prob(p1, ay) / prob(p1, a2)
  if (estimation == "IPW") return(as.numeric(w2 * data$Y))
  as.numeric(q0 + w0 * (q1 - q0) + w1 * (q2 - q1) + w2 * (data$Y - q2))
}

# A fixed two-treatment regime: A -> M1 -> A2 -> M2 -> Y.
reference_multi <- function(data, a1, a2, estimation) {
  q1 <- predict(lm(Y ~ X + A + M1 + A2, data), transform(data, A = a1, A2 = a2))
  q0 <- predict(lm(q1 ~ X + A, data), transform(data, A = a1))
  if (estimation == "G") return(as.numeric(q0))
  p1 <- predict(glm(A ~ X, data, family = binomial()), type = "response")
  p2 <- predict(glm(A2 ~ X + A + M1, data, family = binomial()), type = "response")
  if (a1 == 0) p1 <- 1 - p1
  if (a2 == 0) p2 <- 1 - p2
  w1 <- (data$A == a1) / p1
  w2 <- (data$A == a1 & data$A2 == a2) / (p1 * p2)
  if (estimation == "IPW") return(as.numeric(w2 * data$Y))
  as.numeric(q0 + w1 * (q1 - q0) + w2 * (data$Y - q1))
}

reference_effect <- function(p1, p0, scale = "diff", CI_level = 0.95, estimation = "EIF") {
  m1 <- mean(p1); m0 <- mean(p0)
  if (scale == "risk") return(c(Effect = m1 / m0))
  if (scale == "oddsratio") return(c(Effect = (m1 / (1 - m1)) / (m0 / (1 - m0))))
  if (estimation != "EIF") return(c(Effect = m1 - m0))
  SE <- sd(p1 - p0) / sqrt(length(p1))
  c(Effect = m1 - m0, SE = SE,
    CI.lower = m1 - m0 + qnorm((1 - CI_level) / 2) * SE,
    CI.upper = m1 - m0 + qnorm((1 + CI_level) / 2) * SE)
}

# The two documented path decompositions, ordered as M1, M2, direct, total.
reference_paths <- function(data, estimation, decomposition, scale = "diff") {
  assign <- if (decomposition == "refer0")
    list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1), c(0, 0, 0)) else
    list(c(1, 1, 1), c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
  p <- lapply(assign, function(a) reference_single(data, a, estimation))
  if (decomposition == "refer0") {
    p1 <- p[1:4]; p0 <- rep(p[5], 4)
  } else {
    p1 <- p[c(1, 2, 3, 1)]; p0 <- p[c(2, 3, 4, 4)]
  }
  do.call(rbind, lapply(seq_len(4), function(i)
    reference_effect(p1[[i]], p0[[i]], scale, estimation = estimation)))
}
