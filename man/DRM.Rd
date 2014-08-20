\name{DRM}
\alias{DRM}

\title{
Estimation of dichotomous logistic Rasch model (Rasch, 1960)
}
\description{
Estimation of dichotomous logistic Rasch model (Rasch, 1960).

Parameters are estimated by CML.
}

\usage{
DRM(data, desmat, start)
}

\arguments{
  \item{data}{
Data matrix or data frame; rows represent observations (persons), columns represent the items. 
}
\item {desmat}{
Design matrix. If missing, a design matrix for the Rasch model is automatically computed.
}
  \item{start}{
Starting values for parameter estimation. If missing, a vector of 0 is used as starting values.
}
}


\details{
.
}
\value{
\item{data}{data matrix according to the input}
\item{data_p}{data matrix with data transformed to a response interval between 0 and 1}
\item{itempar}{estimated item parameters}
\item{itempar_se_low}{estimated lower boundary for standard errors of estimated item parameters}
\item{itempar_se_up}{estimated upper boundary for standard errors of estimated item parameters}
\item{itempar_se}{estimated mean standard errors of estimated item parameters}
\item{distrpar}{estimated distribution parameter}
\item{distrpar_se_low}{estimated lower boundary for standard errors of estimated distribution parameter}
\item{distrpar_se_up}{estimated upper boundary for standard errors of estimated distribution parameter}
\item{itempar_se}{estimated mean standard errors of estimated distribution parameter}
\item{iterations}{Number of Newton-Raphson iterations for each item pair}
\item{call}{call of the CRSM function}

}
\references{
Mueller, H. (1987). A Rasch model for continuous ratings. Psychometrika, 52, 165-181.


Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und kontinuierliche Ratingskalen. [Probabilistic models for discrete and continuous rating scales]. Bern: Huber.

}
\author{
Christine Hohensinn
}

\examples{
#estimate CRSM item parameters
#data(example1)
#res_crsm <- CRSM(example1, min=0, max=1)

#summary(res_crsm)

}
\keyword{continuous rating scale model}

