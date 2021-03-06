% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/person_par.CRSM.R, R/person_par.MPRM.R,
%   R/person_par.R
\name{person_par.CRSM}
\alias{person_par.CRSM}
\alias{person_par.MPRM}
\alias{person_par}
\alias{person_par.MPRM}
\alias{person_par.CRSM}
\title{Estimation of person parameters}
\usage{
\method{person_par}{CRSM}(object, ...)

\method{person_par}{MPRM}(object, ..., set0 = FALSE)

person_par(object, ...)
}
\arguments{
\item{object}{Object of class \code{MPRM} or \code{CRSM}}

\item{\dots}{\dots{}}

\item{set0}{if \code{set0=TRUE} for those raw scores patterns with 0
observations (except in the reference category) the person parameter value
is set minimal. With this procedure it is possible to estimate at least the
remaining person parameters of these raw score pattern.  Note: only relevant
for person parameter estimation of MPRM. The person parameters for each raw score vector are constrained to sum zero}
}
\value{
\item{ptable}{table showing for each (observed) raw score the
corresponding estimated person parameter and standard error}
\item{pparList}{for each person raw score, estimated person parameter and
the standard error is displayed} \item{fun_calls}{number of function calls}
\item{call}{function call}
}
\description{
This function performs the estimation of person parameters for the
multidimensional polytomous Rasch model or the continuous Rating Scale
model.
}
\details{
The estimation is performed by Maximum Likelihood Estimation. Thus,
parameters for extreme scores are not calculated!
}
\examples{

#estimate CRSM for the first four items
data(analog)
res_cr <- CRSM(extraversion, low=-10, high=10)

#estimate person parameters for CRSM
pp <- person_par(res_cr)


}
\references{
Fischer, G. H. (1974). Einfuehrung in die Theorie
psychologischer Tests [Introduction to test theory]. Bern: Huber.

Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und
kontinuierliche Ratingskalen. [Probabilistic models for discrete and
continuous rating scales]. Bern: Huber.
}
\seealso{
\code{\link{CRSM}}
}
\author{
Christine Hohensinn
}
\keyword{parameter}
\keyword{person}
