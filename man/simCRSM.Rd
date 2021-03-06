% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simCRSM.R
\name{simCRSM}
\alias{simCRSM}
\title{simulate data according to CRSM}
\usage{
simCRSM(itempar, disp, perspar, mid = 0.5, len = 1, seed = NULL)
}
\arguments{
\item{itempar}{a numerical vector with item parameters}

\item{disp}{a number setting the dispersion parameter for the item set}

\item{perspar}{a numerical vector with the person parameters}

\item{mid}{the midpoint of the response scale (on which the data set is
generated)}

\item{len}{the length of the response scale (on which the data set is
generated)}

\item{seed}{a seed for the random number generated can optionally be set}
}
\value{
\item{datmat}{simulated data set} \item{true_itempar}{the fixed item
parameters according to the input} \item{true_disppar}{the fixed
dispersion parameter according to the input} \item{true_perspar}{the fixed
person parameters according to the input}
}
\description{
With this function data sets according to the Continous Rating Scale Model
are simulated
}
\details{
The midpoint and the length of the response scale define the interval of the
data set generated. The default of the function generates data according to
a response scale between 0 and 1 - that is midpoint 0.5 and length 1.
}
\examples{

#set item parameters
item_p <- c(-1.5,-0.5,0.5,1)

#set dispersion parameter for items
dis_p <- 5

#generate person parameters by a standard normal dispersion
pp <- rnorm(50, 0,1)

#simulate data set
#this is only an illustrating example for simulating data!
#In practice, a sample size of n=50 will be too small for most application
#demands
simdatC <- simCRSM(item_p, dis_p, pp)


}
\references{
Mueller, H. (1987). A Rasch model for continuous ratings.
Psychometrika, 52, 165-181.
}
\seealso{
\code{\link{simMPRM}}
}
\author{
Christine Hohensinn
}
\keyword{continuous}
\keyword{model}
\keyword{rating}
\keyword{scale}
\keyword{simulation}
