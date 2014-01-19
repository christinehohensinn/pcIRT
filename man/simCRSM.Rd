\name{simCRSM}
\alias{simCRSM}

\title{
simulate data according to CRSM
}

\description{
With this function data sets according to the Continous Rating Scale Model are simulated
}
\usage{
simCRSM(itempar, distr, perspar, mid = 0.5, len = 1, seed=NULL)
}
\arguments{
  \item{itempar}{
a numerical vector with item parameters 
}
  \item{distr}{
a number setting the distribution parameter for the item set
}
  \item{perspar}{
a numerical vector with the person parameters
}
  \item{mid}{
the midpoint of the response scale (on which the data set is generated)
}
  \item{len}{
the length of the response scale (on which the data set is generated)
}
  \item{seed}{
a seed for the random number generated can optionally be set
}

}
\details{
The midpoint and the length of the response scale define the interval of the data set generated. The default of the function generates data according to a response scale between 0 and 1 - that is midpoint 0.5 and length 1.

}
\value{
\item{datmat}{simulated data set}
\item{true_itempar}{the fixed item parameters according to the input}
\item{true_distrpar}{the fixed distribution parameter according to the input}
\item{true_perspar}{the fixed person parameters according to the input}
}

\references{
Mueller, H. (1987). A Rasch model for continuous ratings. Psychometrika, 52, 165-181.
}

\author{
Christine Hohensinn
}

\seealso{
\code{\link{simMPRM}}
}

\examples{
#set item parameters
item_p <- c(-1.5,-0.5,0.5,1)

#set distribution parameter for items
dis_p <- 5

#generate person parameters by a standard normal distribution
pp <- rnorm(50, 0,1)

#simulate data set
#this is only an illustrating example for simulating data!
#In practice, a sample size of n=50 will be too small for most application demands
simdatC <- simCRSM(item_p, dis_p, pp)

}
\keyword{continuous rating scale model}
\keyword{simulation}
