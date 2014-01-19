\name{pcIRT-package}
\alias{pcIRT-package}
\alias{pcIRT}
\docType{package}
\title{
IRT models for polytomous and continuous item responses
}
\description{
This package estimates the multidimensional polytomous Rasch model (Rasch, 1961) and provides functions to set linear restrictions on the item category parameters of this models. With this functions it is possible to test whether item categories can be collapsed or set as linear dependent. Thus it is also possible to test whether the multidimensional model can be reduced to a unidimensional model that is whether item categories represent a unidimensional continuum. For this case the scoring parameter of the categories is estimated.


This package estimates the Continuous Rating Scale model by Mueller (1987). It is an extension of the Rating Scale Model by Andrich (1978) on continuous responses (e.g. taken by a visual analog scale).
}
\details{
\tabular{ll}{
Package: \tab pcIRT\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-11-13\cr
License: \tab GPL-3\cr
}

}
\author{
Christine Hohensinn
Maintainer: Christine Hohensinn <christine.hohensinn@univie.ac.at>
}
\references{
Andersen, E. B. (1995). Polytomous Rasch models and their estimation. In G. H. Fischer and I. Molenaar (Eds.). Rasch Models - Foundations, Recent Developements, and Applications. Springer.


Fischer, G. H. (1974). Einfuehrung in die Theorie psychologischer Tests [Introduction to test
theory]. Bern: Huber.

Mueller, H. (1987). A Rasch model for continuous ratings. Psychometrika, 52, 165-181.


Rasch, G. (1961). On general laws and the meaning of measurement in psychology, Proceedings Fourth Berekely Symposium on Mathematical Statistiscs and Probability 5, 321-333.

}

\keyword{ package }
\keyword{ IRT }
\keyword{ Item Response Theory }
\keyword{ psychometrics }
\keyword{ multidimensional polytomous Rasch model }
\keyword{ Continuous Rating Scale Model }

\seealso{
\code{\link{MPRM}}
\code{\link{CRSM}}
}
\examples{
#simulate data set according to the multidimensional polytomous Rasch model (MPRM)
simdat <- simMPRM(rbind(matrix(c(-1.5,0.5,0.5,1,0.8,-0.3, 0.2,-1.2), ncol=4),0), 500)

#estimate MPRM item parameters
res_mprm <- MPRM(simdat$datmat)

summary(res_mprm)

}
