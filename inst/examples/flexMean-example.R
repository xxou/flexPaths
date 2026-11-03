# ******************************************************************************
#    Potential outcome for flexible pathway in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(data = singTreat, A = "Treat", Y = "Outcome1", cov_x = c("X1", "X2"),
                    M.list = list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"),
                                  M3 = "Med3"),
                    estimation = "EIF",
                    model.outcome = ~glm(family = gaussian()),
                    model.treatment = ~glm(family = binomial())
)


# EIF potential means include SE.
# assign gives treatment values for (M1, M2, M3, Y), not mediator values.
# M2 is one stage containing both Med2_1 and Med2_2.
# c(0, 0, 0, 0) uses treatment 0 throughout.
# c(1, 0, 1, 1) uses treatment 1 for M1, M3 and Y, and 0 for M2.
potential_outcome0 <- flexMean(pathsInfo = EIF_fit, assign = c(0, 0, 0, 0))
potential_outcome1 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 1, 1))
potential_outcome2 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 0, 1))
potential_outcome3 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 0, 0))

# G and IPW potential means omit SE.
G_fit <- pathsInfo(
  data = singTreat[1:200, ], A = "Treat", Y = "Outcome1",
  cov_x = c("X1", "X2"), M.list = list(M1 = "Med1"),
  estimation = "G", model.outcome = ~lm()
)
# With one mediator stage, c(1, 0) uses treatment 1 for M1 and 0 for Y.
flexMean(G_fit, c(1, 0))
# Use nboot for bootstrap SE; 10 is for illustration only.
set.seed(123)
flexMean(G_fit, c(1, 0), nboot = 10)
