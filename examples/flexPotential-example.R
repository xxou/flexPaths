# ******************************************************************************
#    Potential outcome for flexible pathway in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                     estimation = "EIF",
                     model.outcome = ~ glm(family = gaussian()),
                     model.treatment = ~ glm(family = binomial())
)


potential_outcome0 <- flexPotential(pathsInfo = EIF_fit, active = c(0, 0, 0, 0))
potential_outcome1 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 1, 1))
potential_outcome2 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 1))
potential_outcome3 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 0))
