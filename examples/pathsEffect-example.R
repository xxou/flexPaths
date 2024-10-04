# ******************************************************************************
#    Path specific effects via each mediator in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                     estimation = "EIF",
                     model.outcome = ~ glm(family = gaussian()),
                     model.propensity = ~ glm(family = binomial())
)

# 0 reference decomposition
results_refer0 <- pathsEffect(pathsInfo = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95)
results_refer0

# sequential sequential
results_seq <- pathsEffect(pathsInfo = EIF_fit, decomposition = "sequential", scale = "diff", CI_level = 0.95)
results_seq

# bootstrap and parallel computation
results_boot <- pathsEffect(pathsInfo = EIF_fit, decomposition = "sequential",
                            scale = "diff", CI_level = 0.95, nboot = 10, m.cores = NULL)

results_parallel <- pathsEffect(pathsInfo = EIF_fit, decomposition = "sequential",
                                scale = "diff", CI_level = 0.95, nboot = 10, m.cores = 3)


