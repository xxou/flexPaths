% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mflexPotential.R
\name{mflexPotential}
\alias{mflexPotential}
\title{Causal Path-Specific Potential Outcome Calculation for Multiple Treatments}
\usage{
mflexPotential(pathsInfo, active)
}
\arguments{
\item{pathsInfo}{An object of class \code{pathsInfo}, which is the output from the \code{pathsInfo} function, containing the necessary fitted models.}

\item{active}{A list of vectors specifying the active values for each mediator (M1 to Mk) and the direct value setting for the outcome (Y).
Each vector in the list corresponds to the setting of one treatment for the mediators. The order of treatments in the list should match the order of mediators.}

\item{data}{A data frame containing all the variables required for the analysis.}
}
\value{
An object of class \code{mflexPotential}, which includes the following components:
\describe{
\item{\code{call}}{The matched call to the \code{mflexPotential} function.}
\item{\code{pathsInfo}}{The original \code{pathsInfo} object used for the analysis.}
\item{\code{active}}{A matrix representing the active values used for each mediator and the direct outcome.}
\item{\code{results}}{A data frame containing the average potential outcome (value) and its standard error (SE).}
\item{\code{potential_data}}{A vector containing the calculated potential outcomes for each observation in the dataset.}
}
}
\description{
This function calculates causal path-specific potential outcomes for multiple treatments.
Users can specify different levels of each treatment for each mediator (M1 to Mk) and directly define the outcome (Y) to represent the potential outcome for a single path or several combined paths.
The function supports various estimation methods, including EIF (Efficient Influence Function), IPW (Inverse Probability Weighting), and G-computation.
}
\examples{
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
}
