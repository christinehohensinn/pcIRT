% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRSM.R, R/print.CRSM.R, R/summary.CRSM.R
\name{CRSM}
\alias{CRSM}
\alias{summary.CRSM}
\alias{print.CRSM}
\alias{print.CRSM}
\alias{summary.CRSM}
\title{Estimation of Continuous Rating Scale Model (Mueller, 1987)}
\usage{
CRSM(data, low, high, start, conv = 1e-04)

\method{print}{CRSM}(x, ...)

\method{summary}{CRSM}(object, ...)
}
\arguments{
\item{data}{Data matrix or data frame; rows represent observations
(persons), columns represent the items.}

\item{low}{The minimum value of the response scale (on which the data are
based).}

\item{high}{The maximum value of the response scale (on which the data are
based).}

\item{start}{Starting values for parameter estimation. If missing, a vector
of 0 is used as starting values.}

\item{conv}{Convergence criterium for parameter estimation.}

\item{x}{object of class \code{CRSM}}

\item{\dots}{\dots{}}

\item{object}{object of class \code{CRSM}}
}
\value{
\item{data}{data matrix according to the input} \item{data_p}{data
matrix with data transformed to a response interval between 0 and 1}
\item{itempar}{estimated item parameters} \item{itempar_se_low}{estimated
lower boundary for standard errors of estimated item parameters}
\item{itempar_se_up}{estimated upper boundary for standard errors of
estimated item parameters} \item{itempar_se}{estimated mean standard errors
of estimated item parameters} \item{disppar}{estimated dispersion
parameter} \item{disppar_se_low}{estimated lower boundary for standard
errors of estimated dispersion parameter} \item{disppar_se_up}{estimated
upper boundary for standard errors of estimated dispersion parameter}
\item{itempar_se}{estimated mean standard errors of estimated item
parameter} \item{disp_est}{estimated dispersion parameters for all item pairs}\item{iterations}{Number of Newton-Raphson iterations for each
item pair} \item{low}{minimal data value entered in call}\item{high}{maximal data value entered in call}\item{call}{call of the CRSM function}
}
\description{
Estimation of the Rating Scale Model for continuous data by Mueller (1987).
}
\details{
\deqn{P_{vi}(a \leq X \leq b) = \frac{\int_a^b exp[x \mu + x(2c-x) \theta]
dx}{\int_{c-\frac{d}{2}}^{c+\frac{d}{2}} exp[t \mu + t(2c-t) \theta] dt}}

Parameters are estimated by a pairwise conditional likelihood estimation (a pseudo-likelihood approach, described in Mueller, 1999).

The parameters of the Continuous Rating Scale Model are estimated by a pairwise cml approach using Newton-Raphson iterations for optimizing.
}
\examples{

#estimate CRSM item parameters
data(analog)
res_crsm <- CRSM(extraversion, low=-10, high=10)

summary(res_crsm)

}
\references{
Mueller, H. (1987). A Rasch model for continuous ratings.
Psychometrika, 52, 165-181.

Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und
kontinuierliche Ratingskalen. [Probabilistic models for discrete and
continuous rating scales]. Bern: Huber.
}
\author{
Christine Hohensinn
}
\keyword{continuous}
\keyword{model}
\keyword{rating}
\keyword{scale}
