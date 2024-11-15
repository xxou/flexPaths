# ****************************************************************************************
#    Path-specific effect for flexible pathway in one single/multiple treatment(s)
# ****************************************************************************************

###########################
# single treatment
###########################
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

# one comparison
flex_results1 <- flexEffect(p1 = potential_outcome1,
                            p0 = potential_outcome0, scale = "diff", CI_level = 0.95,nboot = 100, m.cores = 3)
flex_results1

# multiple comparisons
# PSE: potential_outcome1 - potential_outcome0
# PSE: potential_outcome2 - potential_outcome3
flex_results <- flexEffect(p1 = list(potential_outcome1, potential_outcome2),
                           p0 = list(potential_outcome0, potential_outcome3), scale = "diff",
                           CI_level = 0.95,nboot = 100, m.cores = 4)
flex_results

###########################
# multiple treatments
###########################
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
                  model.treatment =list( ~  glm(family = binomial())),
                  model.outcome = list( ~SuperLearner(SL.library = "SL.mean",family = gaussian())),
                  model.iter  = list(~glm(family = gaussian()))
)

mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)

flexEffect(p1 = mp1, p0 = mp2, scale = "diff", CI_level = 0.95, nboot =100 , m.cores = 6)
