# ******************************************************************************
#    Potential outcome for flexible pathway in multiple treatments
# ******************************************************************************
data("multiTreat")

mfit <- pathsInfo(data = multiTreat,
                  Y = "Y", A = c("T1", "T2", "T3"), cov_x = "X",
                  M.list = list(M1 = "M1", M2 = "M2", M3 = "M3",
                                M4 = "M4", M5 = "M5", M6 = "M6"),
                  estimation = "EIF",
                  model.treatment = list(~glm(family = binomial())),
                  model.outcome = list(
                    ~SuperLearner(SL.library = "SL.mean", family = gaussian())
                  ),
                  model.iter = list(~glm(family = gaussian()))
)
# EIF potential means include SE.
# assign contains one vector for each treatment: T1, T2 and T3, in that order.
# Each vector gives treatment values for (M1, M2, M3, M4, M5, M6, Y).
# NA marks stages before that treatment; here T3 starts after M3.
mp1 <- mflexMean(assign = list(a1 = c(0, 1, 0, 0, 1, 0, 0),
                                   a2 = c(1, 0, 1, 1, 1, 0, 0),
                                   a3 = c(NA, NA, NA, 1, 0, 0, 0)),
                      mfit)
# mp2 changes T1's assignment for M2 and T2's assignments for M1 and M3.
mp2 <- mflexMean(assign = list(a1 = c(0, 0, 0, 0, 1, 0, 0),
                                   a2 = c(0, 0, 0, 1, 1, 0, 0),
                                   a3 = c(NA, NA, NA, 1, 0, 0, 0)),
                      mfit)
# G/IPW omit direct SE; use nboot for bootstrap SE with any method.
# A small IPW example; 10 bootstrap samples are for illustration only.
IPW_fit <- pathsInfo(
  data = multiTreat[1:200, ], A = c("T1", "T2"), Y = "Y", cov_x = "X",
  M.list = list(M1 = "M1", M2 = "M2"), estimation = "IPW",
  model.treatment = ~glm(family = binomial()))
set.seed(123)
# Each vector is ordered as (M1, M2, Y); both treatments are set to 1 throughout.
mflexMean(IPW_fit, list(c(1, 1, 1), c(1, 1, 1)), nboot = 10)
