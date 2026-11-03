# ******************************************************************************
#    Potential outcome for flexible pathway in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(data = singTreat, A = "Treat", Y = "Outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "Med1", M2 = c('Med2_1', 'Med2_2'), M3 = 'Med3'),
                     estimation = "EIF",
                     model.outcome = ~ glm(family = gaussian()),
                     model.treatment = ~ glm(family = binomial())
)


# EIF potential means include SE.
potential_outcome0 <- flexPotential(pathsInfo = EIF_fit, active = c(0, 0, 0, 0))
potential_outcome1 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 1, 1))
potential_outcome2 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 1))
potential_outcome3 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 0))
potential_outcome1$results

# G and IPW potential means omit SE.
G_fit <- pathsInfo(
  data = singTreat[1:200, ], A = "Treat", Y = "Outcome1",
  cov_x = c("X1", "X2"), M.list = list(M1 = "Med1"),
  estimation = "G", model.outcome = ~lm()
)
flexPotential(G_fit, c(1, 0))$results
# Use nboot for bootstrap SE; 10 is for illustration only.
set.seed(123)
flexPotential(G_fit, c(1, 0), nboot = 10)
