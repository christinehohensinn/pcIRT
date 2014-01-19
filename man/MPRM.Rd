\name{MPRM}
\alias{MPRM}
\alias{summary.MPRM}
\alias{print.MPRM}

\title{
Estimation of Multidimensional Polytomous Rasch model (Rasch, 1961)
}
\description{
This function estimates the multidimensional polytomous Rasch model by Rasch (1961). 
The model estimates item category parameters \eqn{\beta} for each item and each category and takes each category of data as another dimension. 


\deqn{P_{vih} = \frac{exp(\theta_{vh} + \beta_{ih})}{\sum_{j=0}^{m} exp(\theta_{vj} + \beta_{ij})}}


Parameters are estimated by CML method.
}

\usage{
MPRM(data, desmat, start)

\method{summary}{MPRM}(object, \dots)
\method{print}{MPRM}(x, \dots)

}

\arguments{
  \item{data}{
Data matrix or data frame; rows represent observations (persons), columns represent the items 
}
  \item{desmat}{
Design matrix; if missing, the design matrix will be created automatically.
}
  \item{start}{
Starting values for parameter estimation. If missing, a vector of 0 is used as starting values.
}
  \item{object}{object of class \code{MPRM}}
  \item{x}{{object of class \code{MPRM}}}  
  \item{\dots}{\dots}

}
\details{
The parameters of the multidimensional polytomous Rasch model (Rasch, 1961) are estimated by CML estimation. For the CML estimation no assumption on the person parameter distribution is necessary. 
}
\value{
\item{data}{data matrix according to the input}
\item{design}{design matrix either according to the input or according to the automatically generated matrix}
\item{logLikelihood}{conditional log-likelihood}
\item{estpar}{estimated basic item category parameters}
\item{estpar_se}{estimated standard errors for basic item category parameters}
\item{itempar}{estimated item category parameters}
\item{itempar_se}{estimated standard errors for item category parameters}
\item{hessian}{Hessian matrix}
\item{convergence}{convergence of solution (see help files in \code{\link{optim}})}
\item{fun_calls}{number of function calls (see help files in \code{\link{optim}})}
}

\references{
Andersen, E. B. (1995). Polytomous Rasch models and their estimation. In G. H. Fischer and I. Molenaar (Eds.). Rasch Models - Foundations, Recent Developements, and Applications. Springer.


Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.


Rasch, G. (1961). On general laws and the meaning of measurement in psychology, Proceedings Fourth Berekely Symposium on Mathematical Statistiscs and Probability 5, 321-333.
}

\author{
Christine Hohensinn
}

\examples{
#simulate data set
simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)

#estimate MPRM item parameters
res_mprm <- MPRM(simdat$datmat)

summary(res_mprm)

}
\keyword{multidimensional polytomous Rasch model}
