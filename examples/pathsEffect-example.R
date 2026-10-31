# ******************************************************************************
#    Path specific effects via each mediator in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(data = singTreat, A = "Treat", Y = "Outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "Med1", M2 = c('Med2_1', 'Med2_2'), M3 = 'Med3'),
                     estimation = "EIF",
                     model.outcome = ~ glm(family = gaussian()),
                     model.treatment = ~ glm(family = binomial())
)

# 0 reference decomposition
results_refer0 <- pathsEffect(pathsInfo = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95)
results_refer0

# sequential sequential
results_seq <- pathsEffect(pathsInfo = EIF_fit, decomposition = "sequential", scale = "diff", CI_level = 0.95)
results_seq

# For a small bootstrap example, see help("flexEffect").


