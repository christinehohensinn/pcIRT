\name{CRSM}
\alias{CRSM}
\alias{summary.CRSM}
\alias{print.CRSM}

\title{
Estimation of Continuous Rating Scale Model (Mueller, 1987)
}
\description{
Estimation of the Rating Scale Model for continuous data by Mueller (1987).


\deqn{P_{vi}(a \leq X \leq b) = \frac{\int_a^b exp[x \mu + x(2c-x) \theta] dx}{\int_{c-\frac{d}{2}}^{c+\frac{d}{2}} exp[t \mu + t(2c-t) \theta] dt}}


Parameters are estimated by a pairwise algorithm.
}

\usage{
CRSM(data, start, min, max)

\method{summary}{CRSM}(object, \dots)
\method{print}{CRSM}(x, \dots)
}

\arguments{
  \item{data}{
Data matrix or data frame; rows represent observations (persons), columns represent the items. 
}
  \item{start}{
Starting values for parameter estimation. If missing, a vector of 0 is used as starting values.
}
  \item{min}{
The minimum value of the response scale (on which the data are based).
}
  \item{max}{
The maximum value of the response scale (on which the data are based).
}
  \item{object}{object of class \code{CRSM}}
  \item{x}{{object of class \code{CRSM}}}  
  \item{\dots}{\dots}
}


\details{
The parameters of the Continuous Rating Scale Model are estimated by a pairwise algorithm (the algorithm is described in detail in Mueller, 1999) using Newton-Raphson iterations for optimizing. For the pairwise algorithm no assumption on the person parameter distribution is necessary. 
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
data(example1)
res_crsm <- CRSM(example1, min=0, max=1)

summary(res_crsm)

}
\keyword{continuous rating scale model}

