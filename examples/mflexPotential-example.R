# ******************************************************************************
#    Potential outcome for flexible pathway in multiple treatments
# ******************************************************************************

data("multiTreat")

mfit<- pathsInfo(data = multiTreat,
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
                 model.propensity =list( ~  glm(family = binomial())),
                 model.outcome = list( ~SuperLearner(SL.library = "SL.mean",family = gaussian())),
                 model.iter  = list(~glm(family = gaussian()))
)

mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
