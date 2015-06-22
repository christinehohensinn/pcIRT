#'@rdname iccplot
#'@method iccplot CRSM
#'@export

iccplot.CRSM <- function(object, ...){

  pp <- seq(-30,30, by=1)
  rs <- seq(0,1, by=0.05)
  pp_rs <- expand.grid(theta = pp, resp = rs)
  func <- function(x, per, itpar){exp(x*(per-itpar) + x*(1-x)*object$disppar)}

  par(ask=TRUE)
  z2 <- sapply(1:length(object$itempar), function(l){
      resp <- sapply(1:nrow(pp_rs), function(z1){
        if(pp_rs[z1,2] == 0){
          z <- integrate(func, per=pp_rs[z1,1], itpar=object$itempar[l], lower=0.0000001, upper=0.01)$value
        }  else {
          z <- integrate(func, per=pp_rs[z1,1], itpar=object$itempar[l], lower=(pp_rs[z1,2]-0.01), upper=pp_rs[z1,2])$value
        }
          n <- integrate(func, per=pp_rs[z1,1], itpar=object$itempar[l],lower=0, upper=1)$value
      z/n
      })
      mat1 <- matrix(resp, nrow=length(pp))
      oldscale <- rs*(object$high-object$low)+object$low
      persp(pp, oldscale, mat1, xlab="theta", ylab="response scale", zlab="response probability", theta=35, phi=25, ticktype="detailed", main=names(object$itempar)[l])
  })
  par(ask=FALSE)

}
