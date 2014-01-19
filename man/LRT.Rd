\name{LRT}
\alias{LRT}
\alias{summary.aLR}
\alias{print.aLR}

\title{
Computes Andersen's Likelihood Ratio Test for the multidimensional polytomous Rasch model
}
\description{
Andersen's Likelihood Ratio Test is a model test for Rasch models (based on CML estimation) and splits the data set into subsamples to test the person homogeneity

}
\usage{
LRT(MPRMobj, splitcrit = "score")

\method{summary}{aLR}(object, \dots)
\method{print}{aLR}(x, \dots)

}

\arguments{
  \item{MPRMobj}{
Object of class \code{MPRM}
}
  \item{splitcrit}{
Vector or the character vector \code{"score"} to define the split criterion. The default split criterion \code{"score"} splits the sample according to the median of the raw score. Vector can be numeric, factor or character. (see details)
}
  \item{object}{object of class \code{aLR}}
  \item{x}{{object of class \code{aLR}}} 
  \item{\dots}{\dots}
}

\details{
The default split criterion \code{"score"} computes the raw score of every person according to the category values in the data set. The sample is split by the median of this raw score. 
}
\value{
\item{emp_Chi2}{\eqn{\chi^2} distributed value of the Likelihood Ratio test}
\item{df}{degrees of freedom of the test statistic}
\item{pval}{p value of the test statistic}
\item{itempar}{estimated item parameters for each subsample}
\item{item_se}{estimated standard errors for the item parameters for each subsample}
}
\references{
Andersen, E. B. (1973). A goodness of fit test for the Rasch model. Psychometrika, 38, 123-
140.

Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.


}
\author{
Christine Hohensinn
}

\seealso{
\code{\link{MPRM}}
\code{\link{dLRT}}
}
\examples{
#simulate data set
simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)

#estimate MPRM item parameters
res_mprm <- MPRM(simdat$datmat)

#compute Andersen's Likelihood Ratio test
res_lrt <- LRT(res_mprm)
summary(res_lrt)
  }

\keyword{Likelihood Ratio test}
\keyword{model test}
