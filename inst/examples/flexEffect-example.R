# ****************************************************************************************
#    Path-specific effect for flexible pathway in one single/multiple treatment(s)
# ****************************************************************************************

###########################
# single treatment
###########################
data("singTreat")
# input information for PSE
EIF_fit <- pathsInfo(data = singTreat[1:200, ], A = "Treat", Y = "Outcome1",
                     cov_x = c("X1", "X2"),
                     M.list = list(M1 = "Med1", M2 = c('Med2_1', 'Med2_2'), M3 = 'Med3'),
                     estimation = "EIF",
                     model.outcome = ~ glm(family = gaussian()),
                     model.treatment = ~ glm(family = binomial())
)


potential_outcome0 <- flexPotential(pathsInfo = EIF_fit, active = c(0, 0, 0, 0))
potential_outcome1 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 1, 1))
potential_outcome2 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 1))
potential_outcome3 <- flexPotential(pathsInfo = EIF_fit, active = c(1, 0, 0, 0))

# one comparison; 10 bootstrap samples are used only to illustrate the syntax
set.seed(123)
flex_results1 <- flexEffect(p1 = potential_outcome1,
                            p0 = potential_outcome0, scale = "diff",
                            CI_level = 0.95, nboot = 10, m.cores = NULL)
flex_results1

# EIF differences include SE, confidence intervals and p-values.
# multiple comparisons
# PSE: potential_outcome1 - potential_outcome0
# PSE: potential_outcome2 - potential_outcome3
flex_results <- flexEffect(p1 = list(potential_outcome1, potential_outcome2),
                           p0 = list(potential_outcome0, potential_outcome3), scale = "diff",
                           CI_level = 0.95)
flex_results

# G and IPW effects return point estimates without direct SE or intervals.
G_fit <- pathsInfo(
  data = singTreat[1:200, ], A = "Treat", Y = "Outcome1",
  cov_x = c("X1", "X2"), M.list = list(M1 = "Med1"),
  estimation = "G", model.outcome = ~lm()
)
g0 <- flexPotential(G_fit, c(0, 0))
g1 <- flexPotential(G_fit, c(1, 1))
flexEffect(g1, g0)
# Use nboot for bootstrap SE and intervals; 10 is for illustration only.
set.seed(123)
flexEffect(g1, g0, nboot = 10)
# On risk and oddsratio scales, use nboot even with EIF estimation.
# Bootstrap p-values are returned only on the difference scale.

###########################
# multiple treatments
###########################
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

mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)

flexEffect(p1 = mp1, p0 = mp2, scale = "diff", CI_level = 0.95)
