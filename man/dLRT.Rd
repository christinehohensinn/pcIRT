\name{dLRT}
\alias{dLRT}
\alias{summary.dLR}
\alias{print.dLR}

\title{
Dimensionality test for the multidimensional polytomous Rasch model
}
\description{
This function tests whether the multidimensional polytomous Rasch model can be reduced to a unidimensional polytomous model. 
}
\usage{
dLRT(MPRMobj)

\method{summary}{dLR}(object, \dots)
\method{print}{dLR}(x, \dots)

}
\arguments{
  \item{MPRMobj}{
Object of class \code{MPRM}
}
  \item{object}{object of class \code{dLR}}
  \item{x}{{object of class \code{dLR}}}  
  \item{\dots}{\dots}
}

\details{
For this test, a unidimensional model assuming the categories as linearly dependent is computed. Subsequently a Likelihood Ratio test is conducted.
}
\value{
\item{emp_Chi2}{\eqn{\chi^2} distributed value of the Likelihood Ratio test}
\item{df}{degrees of freedom of the test statistic}
\item{pval}{p value of the test statistic}
}

\references{
Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.
}
\author{
Christine Hohensinn
}

\seealso{
\code{\link{MPRM}}
\code{\link{LRT}}
}
\examples{
#simulate data set
simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)

#estimate MPRM item parameters
res_mprm <- MPRM(simdat$datmat)

res_dlrt <- dLRT(res_mprm)
summary(res_dlrt)

}
\keyword{ dimensionality }
\keyword{ model test}
