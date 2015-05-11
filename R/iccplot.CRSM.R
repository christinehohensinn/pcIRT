#'@rdname iccplot
#'@method iccplot CRSM
#'@export

iccplot.CRSM <- function(object, ...){

  pp <- seq(-30,30, by=1)
  rs <- seq(0.1,1, by=0.1)
  pp_rs <- expand.grid(theta = pp, resp = rs)
  func <- function(x, per, itpar){exp(x*(per-itpar) + x*(1-x)*object$distrpar)}

  par(ask=TRUE)
  z <- sapply(1:length(object$itempar), function(l){
      resp <- sapply(1:nrow(pp_rs), function(z1){
          z <- integrate(func, per=pp_rs[z1,1], itpar=object$itempar[l], lower=(pp_rs[z1,2]-0.1), upper=pp_rs[z1,2])$value
          n <- integrate(func, per=pp_rs[z1,1], itpar=object$itempar[l],lower=0, upper=1)$value
      z/n
      })
      mat1 <- matrix(resp, nrow=length(pp))
      persp(pp, rs, mat1, xlab="theta", ylab="response scale", zlab="response probability", theta=35, phi=25, ticktype="detailed", main=names(object$itempar)[l])
  })
  par(ask=FALSE)

}
