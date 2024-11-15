# ***************************************
#    single treatment  example
# ***************************************
data("singTreat")
# EIF estimation
EIF_fit <- pathsInfo(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                     estimation = "EIF",
                     model.outcome = ~SuperLearner(SL.library = c('randomForest','xgboost'),family = gaussian()),
                     model.treatment = ~ bart(verbose = FALSE, ndpost = 200)
)



EIF_fit <- pathsInfo(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                     estimation = "EIF",
                     model.outcome = ~ SuperLearner(SL.library = c('randomForest','xgboost'),family = gaussian()),
                     model.treatment = list(cov_x = ~ bart(verbose = FALSE, ndpost = 200),
                                             M1 = ~ glm(family = binomial()),
                                             M2 = ~ SuperLearner(SL.library = c('randomForest','xgboost'),family = binomial()),
                                             M3 = ~ glm(family = binomial()))
)


EIF_fit <- pathsInfo(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                     M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                     estimation = "EIF",
                     model.outcome = list(cov_x = ~ glm(formula = outcome1 ~X1*X2 ,family = gaussian()),
                                          M1 = ~ glm(formula = outcome1 ~med1+X1*X2, family = gaussian()),
                                          M2 = ~ glm(formula = outcome1 ~med1+med2_1*med2_1+X1*X2, family = gaussian()),
                                          M3 = ~ glm(formula = outcome1 ~med1+med2_1*med2_1+med3+X1*X2, family = gaussian())),
                     model.treatment = list( ~ bart(verbose = FALSE, ndpost = 200))
)

# ***************************************
#   multiple treatments fit example
# ***************************************
data("multiTreat")
mfit<- pathsInfo(data=multiTreat,
                Y = "Y",
                A = c("t1","t2","t3"),
                cov_x = "X",
                M.list = list(
                  M1 = 'm1',
                  M2 = 'm2',
                  M3 = 'm3',
                  M4 = 'm4',
                  M5 = 'm5',
                  M6 = 'm6'
                ),
                estimation = "EIF",
                model.treatment =list( ~  glm(family = binomial())),
                model.outcome = list( ~SuperLearner(SL.library = "SL.mean",family = gaussian())),
                model.iter  = list(~glm(family = gaussian()))
)
