README
================
Xiaxian Ou
2024-08-24

- [1. single treatment](#1-single-treatment)
  - [path fit](#path-fit)
  - [mediation effects through each
    mediator](#mediation-effects-through-each-mediator)
  - [flexible paths for single
    treatment](#flexible-paths-for-single-treatment)
- [multiple treatments](#multiple-treatments)
  - [flexible paths for multiple
    treatments](#flexible-paths-for-multiple-treatments)

**flexPaths**: flexible causal path-specific effect

- two decomposition strategies

- flexible paths

- flexible models

- flexible estimators

Install this package:

``` r
install.packages("devtools")

library(devtools)

install_github("xxou/flexPaths")
library(flexPaths)
```

# 1. single treatment

## path fit

``` r
data("singTreat")
EIF_fit <- pathsFit(data = singTreat, A = "treat", Y = "outcome1", cov_x = c("X1", "X2"),
                    M.list = list(M1 = "med1", M2 = c('med2_1', 'med2_2'), M3 = 'med3'),
                    estimation = "EIF",
                    model.outcome = list(~ glm(family = gaussian())),
                    model.propensity = ~ bart(verbose = FALSE, ndpost = 200)
)
```

## mediation effects through each mediator

``` r
effect_results1 <- pathsEffect(pathsFit = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95)
effect_results2 <- pathsEffect(pathsFit = EIF_fit, decomposition = "refer0", scale = "diff", CI_level = 0.95, nboot = 100, m.cores = 6)

effect_results3 <- pathsEffect( pathsFit = EIF_fit, decomposition = "sequential", scale = "diff", CI_level = 0.95)
```

## flexible paths for single treatment

potential outcome:

``` r
potential_outcome0 <- flexPotential(pathsFit = EIF_fit, active = c(0, 0, 0, 0))
potential_outcome1 <- flexPotential(pathsFit = EIF_fit, active = c(1, 0, 1, 1))
potential_outcome2<- flexPotential(pathsFit = EIF_fit, active = c(1, 0, 0, 1))
```

use potential outcome to calculate PSE

``` r
flex_results1 <- flexEffect(p1 = list(potential_outcome1, potential_outcome2),
                            p0 = potential_outcome0, scale = "diff", CI_level = 0.95,nboot = 5)
```

# multiple treatments

## flexible paths for multiple treatments

pathfit

``` r
data("multiTreat")
mfit<- pathsFit(data = multiTreat,
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
```

potential outcome for active setting

``` r
mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
```

PSE

``` r
flexEffect(p1 = mp1, p0 = mp2, scale = "diff", CI_level = 0.95, nboot =2 , m.cores = 8)
```
