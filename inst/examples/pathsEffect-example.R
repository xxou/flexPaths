# ******************************************************************************
#    Path specific effects via each mediator in one single treatment
# ******************************************************************************

data("singTreat")

# input information for PSE
EIF_fit <- pathsInfo(
  data = singTreat,
  A = "Treat",
  Y = "Outcome1",
  cov_x = c("X1", "X2"),
  M.list = list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"), M3 = "Med3"),
  estimation = "EIF",
  model.outcome = ~glm(family = gaussian()),
  model.treatment = ~glm(family = binomial())
)

# EIF differences include SE, confidence intervals and p-values.
# 0 reference decomposition
# Compare each single-stage assignment of 1 with assign = c(0, 0, 0, 0).
results_refer0 <- pathsEffect(
  pathsInfo = EIF_fit,
  decomposition = "refer0",
  scale = "diff", CI_level = 0.95
)

# sequential decomposition
# Compare adjacent assignments: 1111, 0111, 0011, 0001, then 0000,
# where the four positions represent (M1, M2, M3, Y).
results_seq <- pathsEffect(
  pathsInfo = EIF_fit,
  decomposition = "sequential",
  scale = "diff", CI_level = 0.95
)

# G/IPW and risk/oddsratio results contain only point estimates without nboot.
# Use nboot for bootstrap SE and intervals with any estimation method.
# Bootstrap p-values are returned only on the difference scale.
# For a small bootstrap example, see help("flexEffect").
