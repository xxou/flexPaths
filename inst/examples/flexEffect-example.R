# ******************************************************************************
#    Path-specific effect(s) for flexible pathway(s)
#    in single or multiple treatment settings
# ******************************************************************************

# ................
# single treatment
# ................
data("singTreat")
EIF_fit <- pathsInfo(data = singTreat[1:200, ], A = "Treat", Y = "Outcome1",
                    cov_x = c("X1", "X2"),
                    M.list = list(M1 = "Med1",
                                  M2 = c("Med2_1", "Med2_2"),
                                  M3 = "Med3"),
                    estimation = "EIF",
                    model.outcome = ~glm(family = gaussian()),
                    model.treatment = ~glm(family = binomial())
)


# assign specifies the treatment assignment for (M1, M2, M3, Y).
# For example, c(0, 0, 0, 0) assigns treatment 0 to all components,
# while c(1, 0, 1, 1) assigns treatment 1 to M1, M3, and Y,
# and treatment 0 to M2.

potential_outcome0 <- flexMean(pathsInfo = EIF_fit, assign = c(0, 0, 0, 0))
potential_outcome1 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 1, 1))
potential_outcome2 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 0, 1))
potential_outcome3 <- flexMean(pathsInfo = EIF_fit, assign = c(1, 0, 0, 0))

# one comparison; 10 bootstrap samples are used only to illustrate the syntax
set.seed(123)
flex_results1 <- flexEffect(p1 = potential_outcome1,
                            p0 = potential_outcome0,
                            scale = "diff", CI_level = 0.95,
                            nboot = 10, m.cores = NULL)


# multiple comparisons
# PSE: potential_outcome1 - potential_outcome0
# PSE: potential_outcome2 - potential_outcome3
flex_results <- flexEffect(p1 = list(potential_outcome1, potential_outcome2),
                           p0 = list(potential_outcome0, potential_outcome3),
                           scale = "diff", CI_level = 0.95)


# ....................
# multiple treatments
# ....................
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


# Data structure: T1, T2, M1, M2, M3, T3, M4, M5, M6, Y.
# assign is a list of treatment assignment vectors, one for each
# treatment (T1, T2, T3), in that order.
# Each vector specifies the treatment values for (M1, M2, M3, M4, M5, M6, Y).
# NA indicates mediator stages that occur before the corresponding treatment.
mp1 <- mflexMean(assign = list(a1 = c(0,  1,  0,  0, 1, 0, 0),
                                    a2 = c(1,  0,  1,  1, 1, 0, 0),
                                    a3 = c(NA, NA, NA, 1, 0, 0, 0)),
                      mfit)
mp2 <- mflexMean(assign = list(a1 = c(0,  0,  0,  0, 1, 0, 0),
                                    a2 = c(0,  0,  0,  1, 1, 0, 0),
                                    a3 = c(NA, NA, NA, 1, 0, 0, 0)),
                      mfit)

flexEffect(p1 = mp1, p0 = mp2,
           scale = "diff", CI_level = 0.95)
