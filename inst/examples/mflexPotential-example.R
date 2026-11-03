# ******************************************************************************
#    Potential outcome for flexible pathway in multiple treatments
# ******************************************************************************

data("multiTreat")

mfit<- pathsInfo(data = multiTreat,
                 Y = "Y",
                 A = c("T1","T2","T3"),
                 cov_x = "X",
                 M.list = list(
                   M1 = 'M1',
                   M2 = 'M2',
                   M3 = 'M3',
                   M4 = 'M4',
                   M5 = 'M5',
                   M6 = 'M6'
                 ),
                 estimation = "EIF",
                 model.treatment =list( ~  glm(family = binomial())),
                 model.outcome = list( ~SuperLearner(SL.library = "SL.mean",family = gaussian())),
                 model.iter  = list(~glm(family = gaussian()))
)

# EIF potential means include SE.
mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp1$results
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
# With estimation = "G" or "IPW", the results omit SE.
# Use nboot for bootstrap SE with any estimation method.
# A small IPW example; 10 bootstrap samples are for illustration only.
IPW_fit <- pathsInfo(
  data = multiTreat[1:200, ], A = c("T1", "T2"), Y = "Y", cov_x = "X",
  M.list = list(M1 = "M1", M2 = "M2"), estimation = "IPW",
  model.treatment = ~glm(family = binomial())
)
set.seed(123)
mflexPotential(IPW_fit, list(c(1, 1, 1), c(1, 1, 1)), nboot = 10)
