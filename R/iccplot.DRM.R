#'@rdname iccplot
#'@method iccplot DRM
#'@export

iccplot.DRM <- function(object, ...){

pp <- seq(-5,5, by=0.5)

y <- sapply(1:length(object$itempar), function(l){

        resp <- exp(pp-object$itempar[l]) / (1+exp(pp-object$itempar[l]))
        par(ask=TRUE)
        plot(pp, resp, type="l", xlab="theta", ylab="response probability", main=paste0(names(object$itempar)[l]))
        par(ask=FALSE)
  })
}

