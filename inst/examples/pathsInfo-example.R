# ***************************************
#    single treatment example
# ***************************************
data("singTreat")

# EIF estimation
# A single model specification is reused across the relevant stages.
EIF_fit <- pathsInfo(
  data = singTreat,
  A = "Treat",
  Y = "Outcome1",
  cov_x = c("X1", "X2"),
  M.list = list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"), M3 = "Med3"),
  estimation = "EIF",
  model.outcome = ~SuperLearner(
    SL.library = c("randomForest", "xgboost"),
    family = gaussian()
  ),
  model.treatment = ~bart(verbose = FALSE, ndpost = 200)
)

# Named entries select different treatment models for each conditioning set.
# cov_x uses baseline covariates; M2 also includes stages M1 and M2.
EIF_fit <- pathsInfo(
  data = singTreat,
  A = "Treat",
  Y = "Outcome1",
  cov_x = c("X1", "X2"),
  M.list = list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"), M3 = "Med3"),
  estimation = "EIF",
  model.outcome = ~SuperLearner(
    SL.library = c("randomForest", "xgboost"),
    family = gaussian()
  ),
  model.treatment = list(
    cov_x = ~bart(verbose = FALSE, ndpost = 200),
    M1 = ~glm(family = binomial()),
    M2 = ~SuperLearner(
      SL.library = c("randomForest", "xgboost"),
      family = binomial()
    ),
    M3 = ~glm(family = binomial())
  )
)

# Named formulas allow different predictor terms at each stage.
EIF_fit <- pathsInfo(
  data = singTreat,
  A = "Treat",
  Y = "Outcome1",
  cov_x = c("X1", "X2"),
  M.list = list(M1 = "Med1", M2 = c("Med2_1", "Med2_2"), M3 = "Med3"),
  estimation = "EIF",
  model.outcome = list(
    cov_x = ~glm(
      formula = Outcome1 ~ X1 * X2,
      family = gaussian()
    ),
    M1 = ~glm(
      formula = Outcome1 ~ Med1 + X1 * X2,
      family = gaussian()
    ),
    M2 = ~glm(
      formula = Outcome1 ~ Med1 + Med2_1 * Med2_1 + X1 * X2,
      family = gaussian()
    ),
    M3 = ~glm(
      formula = Outcome1 ~ Med1 + Med2_1 * Med2_1 + Med3 + X1 * X2,
      family = gaussian()
    )
  ),
  model.treatment = list(~bart(verbose = FALSE, ndpost = 200))
)

# ***************************************
#   multiple treatments fit example
# ***************************************
data("multiTreat")

mfit <- pathsInfo(
  data = multiTreat,
  Y = "Y",
  A = c("T1", "T2", "T3"),
  cov_x = "X",
  M.list = list(
    M1 = "M1",
    M2 = "M2",
    M3 = "M3",
    M4 = "M4",
    M5 = "M5",
    M6 = "M6"
  ),
  estimation = "EIF",
  model.treatment = list(~glm(family = binomial())),
  model.outcome = list(
    ~SuperLearner(SL.library = "SL.mean", family = gaussian())
  ),
  model.iter = list(~glm(family = gaussian()))
)
