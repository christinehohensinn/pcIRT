\name{MPRMl}
\alias{MPRMl}

\title{
Estimation of Multidimensional Polytomous Rasch model with linear restrictions
}

\description{
Estimation of the multidimensional polytomous Rasch model with setting linear restrictions on the item category parameters \eqn{\beta}. 

Parameter estimations is done by CML method.

}

\usage{
MPRMl(data, desmat, ldes, lp, start)
}


\arguments{
  \item{data}{
Data matrix or data frame; rows represent observations (persons), columns represent the items 
}
  \item{desmat}{
Design matrix
}
  \item{ldes}{
a numeric vector of the same length as the number of item category parameters indicating which parameters are set linear dependent of which other parameters (see details)
  }
  \item{lp}{
a numeric vector with length equal to the number of item parameters set linear dependent. The vector indicates the number of scoring parameters (see details)  
  }
  \item{start}{
Starting values for parameter estimation. If missing, a vector of 0 is used as starting values.
}

}


\details{
With this function linear restrictions can be set on the general multidimensional polytomous Rasch model. Item category parameters can be set as being linear dependent to other item category parameters and the scoring parameter (as the multiple of the linear dependen parameters) is estimated. The restrictions are set by defining the arguments \code{ldes} and \code{lp}. \code{ldes} is a numerical vector of the same length as item category parameters in the general MPRM. A 0 in this vector indicates that no restriction is set. Putting in another number sets the item category parameter according to the vector position as linear dependent to that item category parameter with the position of the number included. For example, if item category parameter of item 1 and category 2 (that is position 2 in the vector \code{ldes}) should be linear dependent to the item category parameter of item 1 and category 1 (that is position 1 in the vector \code{ldes}), than the number 1 has to be on the second element of vector \code{ldes}. With the vector \code{lp} it is set, how many different scoring parameters have to be estimated and (if there are more than two) which of them should be equal. For example if 5 item category parameters are set linear dependent (by \code{ldes}) and according to the \code{ldes} vector the first, third and fourth have the same scoring parameters and the second and fifth have another scoring parameter, than \code{lp} must be a vector \code{lp = c(1,2,1,1,2)}. 

It is necessary that the design matrix is specified in accordance with the restrictions in \code{ldes} and \code{lp}.



}

\value{
\item{data}{data matrix according to the input}
\item{design}{design matrix according to the input}
\item{logLikelihood}{conditional log-likelihood}
\item{estpar}{estimated basic item category parameters}
\item{estpar_se}{estimated standard errors for basic item category parameters}
\item{itempar}{estimated item category parameters}
\item{itempar_se}{estimated standard errors for item category parameters}
\item{linpar}{estimated scoring parameters}
\item{linpar_se}{estimated standard errors for scoring parameters}
\item{hessian}{Hessian matrix}
\item{convergence}{convergence of solution (see help files in \code{\link{optim}})}
\item{fun_calls}{number of function calls (see help files in \code{\link{optim}})}


}
\references{
Andersen, E. B. (1974). Das mehrkategorielle logistische Testmodell [The polytomous logistic test model] In. W. F. Kempf (Ed.), Probabilistische Modelle in der Sozialpsychologie [Probabilistic model in social psychology]. Bern: Huber.


Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.


Rasch, G. (1961). On general laws and the meaning of measurement in psychology, Proceedings Fourth Berekely Symposium on Mathematical Statistiscs and Probability 5, 321-333.
}
\author{
Christine Hohensinn
}


\seealso{
\code{\link{MPRM}}}

\examples{
#simulate data set according to the general MPRM
simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)

#estimate the general MPRM
res_mprm <- MPRM(simdat$datmat)

#estimate a MPRM with linear restrictions; 
#for item 1 and 2 the second category is set linear dependent to the first category
ldes1 <- rep(0,length(res_mprm$itempar))
ldes1[c(2,5)] <- c(1,4)
lp1 <- rep(1,2)
#take the design matrix from the general MPRM and modify it according to the linear restriction
design1 <- res_mprm$design
design1[2,1] <- 1
design1[5,3] <- 1
design1[11,c(1,3)] <- -1
design1 <- design1[,-c(2,4)]

res_mprml <- MPRMl(simdat$datmat, desmat=design1, ldes=ldes1, lp=lp1)

summary(res_mprml)
}
\keyword{multidimensional polytomous Rasch model}
\keyword{linear restriction}
