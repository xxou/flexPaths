README
================
Xiaxian Ou
2024-08-26

- [1. single treatment](#1-single-treatment)
  - [path fit](#path-fit)
  - [mediation effects through each
    mediator](#mediation-effects-through-each-mediator)
  - [flexible paths for single
    treatment](#flexible-paths-for-single-treatment)
- [2. multiple treatments](#2-multiple-treatments)
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

# 2. multiple treatments

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

The length of the active list corresponds to the number of treatments.
In active1 example, each vector in $a_1, a_2, a_3$ represents values
assigned to $M_1, M_2, ..., M_6$ and $Y$. The first three values are
NAs, indicating that the $a_3$ is between $M_3$ and $M_4$. The longest
path in this scenario is
$A_1 \rightarrow A_2 \rightarrow M_1 \rightarrow M_2 \rightarrow M_3 \rightarrow A_3 \rightarrow M_4 \rightarrow M_5 \rightarrow M_6 \rightarrow Y$.

For active2, the longest path is
$A_1  \rightarrow M_1 \rightarrow M_2 \rightarrow A_2 \rightarrow M_3  \rightarrow M_4 \rightarrow A_3 \rightarrow M_5 \rightarrow M_6 \rightarrow Y$.

``` r
active1 = list(a1=c(0,1,0,0,1,0,0),
               a2=c(1,0,1,1,1,0,0),
               a3=c(NA,NA,NA,1,0,0,0))



active2 = list(a1=c(0,1,0,0,1,0,0),
               a2=c(NA,NA,1,1,1,0,0),
               a3=c(NA,NA,NA,NA,0,0,0))
```

Potential outcome for active1 setting.

``` r
mp1<-mflexPotential(active = list(a1=c(0,1,0,0,1,0,0),
                                  a2=c(1,0,1,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
mp2<-mflexPotential(active = list(a1=c(0,0,0,0,1,0,0),
                                  a2=c(0,0,0,1,1,0,0),
                                  a3=c(NA,NA,NA,1,0,0,0)),mfit)
```

PSE : mp1 - mp2

``` r
flexEffect(p1 = mp1, p0 = mp2, scale = "diff", CI_level = 0.95, nboot =2 , m.cores = 8)
```
