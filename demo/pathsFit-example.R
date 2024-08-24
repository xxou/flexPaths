# demo/pathsFit-demo.R

#######################################
#    single treatment fit example
#######################################
data("singTreat")
# EIF estimation
EIF_fit<-pathsFit(
  data = singTreat, A  = "treat", Y="outcome1", cov_x = c("X1","X2"),
  M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'),M3 = 'med3'),
  estimation = "EIF",
  model.outcome = list( ~ glm(family = gaussian()),
                        ~ SuperLearner(SL.library = c("SL.mean"),family = gaussian()),
                        ~ glm(family = gaussian()),
                        ~ bart(family = gaussian())),
  model.propensity = ~ bart(verbose = F,ndpost = 200),
  model.iter = list( ~ glm(family = gaussian()),
                     ~ SuperLearner(SL.library = c("SL.mean"),family = gaussian()),
                     ~ glm(family = gaussian()) )
)

# IPW estimation
IPW_fit<-pathsFit(
  data = singTreat, A  = "treat", Y="outcome1", cov_x = c("X1","X2"),
  M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'),M3 = 'med3'),
  estimation = "IPW",
  model.propensity =list(cov_x = ~ bart(verbose = F,ndpost = 200), # P(A|X)
                         M1 = ~ glm(family = binomial()),          # P(A|M1, X)
                         M2 = ~ glm(family = binomial()),          # P(A|M1,M2, X)
                         M3 = ~ glm(family = binomial()) )         # P(A|M1,M2,M3, X)
)

# G estimation
G_fit<-pathsFit(
  data = singTreat, A  = "treat", Y="outcome1", cov_x = c("X1","X2"),
  M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'),M3 = 'med3'),
  estimation = "G",
  model.outcome = list(~ SuperLearner(SL.library = c("SL.mean"),family = gaussian()))
)


#######################################
#   multiple treatments fit example
#######################################
data("multiTreat")
mfit<- pathsFit(data=multiTreat,
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
