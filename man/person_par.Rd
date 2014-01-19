\name{person_par}
\alias{person_par}
\alias{person_par.MPRM}
\alias{person_par.CRSM}

\title{
Estimation of person parameters
}
\description{
This function performs the estimation of person parameters for the multidimensional polytomous Rasch model or the continuous Rating Scale model.
}
\usage{
person_par(object, \dots)

\method{person_par}{MPRM}(object, \dots, set0=FALSE)
\method{person_par}{CRSM}(object, \dots)

}

\arguments{
  \item{object}{
Object of class \code{MPRM} or \code{CRSM}
}
  \item{set0}{
if \code{set0=TRUE} for those raw scores patterns with 0 observations (except in the reference category) the person parameter value is set minimal. With this procedure it is possible to estimate at least the remaining person parameters of these raw score pattern.  Note: only relevant for person parameter estimation of MPRM.
}
   \item{\dots}{\dots}
}

\details{
The estimation is performed by Maximum Likelihood Estimation. Thus, parameters for extreme scores are not calculated!
}
\value{
\item{ptable}{table showing for each (observed) raw score the corresponding estimated person parameter and standard error}
\item{pparList}{for each person raw score, estimated person parameter and the standard error is displayed}
\item{fun_calls}{number of function calls}
\item{call}{function call}
}


\references{
Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.


Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und kontinuierliche Ratingskalen. [Probabilistic models for discrete and continuous rating scales]. Bern: Huber.

}
\author{
Christine Hohensinn
}

\seealso{
\code{\link{CRSM}}
}
\examples{
#estimate CRSM for the first four items
data(example1)
ex_cr <- CRSM(example1[,1:4], min=0, max=1)

#estimate person parameters for CRSM
pp <- person_par(ex_cr)

}
\keyword{ person parameter}

