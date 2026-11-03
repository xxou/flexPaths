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
  M.list = list(M1 = "Med1",
                M2 = c("Med2_1", "Med2_2"),
                M3 = "Med3"),
  estimation = "EIF",
  model.outcome = ~glm(family = gaussian()),
  model.treatment = ~glm(family = binomial())
)


# 0 reference decomposition: natural path-specific effects
# In flexMean(), assign follows (M1, M2, M3, Y); 1000 means c(1, 0, 0, 0).
# On the difference scale, the mean contrasts are:
#   M1: 1000 - 0000; M2: 0100 - 0000; M3: 0010 - 0000; direct: 0001 - 0000.
# For example, the M1 target contrast can be expressed as:
# flexEffect(p1 = flexMean(EIF_fit, assign = c(1, 0, 0, 0)),
#            p0 = flexMean(EIF_fit, assign = c(0, 0, 0, 0)), scale = "diff")
results_refer0 <- pathsEffect(
  pathsInfo = EIF_fit,
  decomposition = "refer0",
  scale = "diff", CI_level = 0.95
)

# sequential decomposition: cumulative path-specific effects
#   M1: 1111 - 0111; M2: 0111 - 0011; M3: 0011 - 0001; direct: 0001 - 0000.
# These contrasts sum to the total effect, 1111 - 0000.
# Here the M1 comparison keeps the remaining entries at 1, rather than 0:
# flexEffect(p1 = flexMean(EIF_fit, assign = c(1, 1, 1, 1)),
#            p0 = flexMean(EIF_fit, assign = c(0, 1, 1, 1)), scale = "diff")
results_seq <- pathsEffect(
  pathsInfo = EIF_fit,
  decomposition = "sequential",
  scale = "diff", CI_level = 0.95
)

