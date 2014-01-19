print.CRSM <-
function(x, ...){
  
  cat("\n Call: ", deparse(x$call), "\n\n")
  
    parall <- rbind(cbind("item estimates"=x$itempar, "SE"=x$se.item.mean), "lambda"=c(x$distrpar, x$se.distr.mean))
      parall <- rbind(cbind("item estimates"=x$itempar, "SE"=x$itempar_se),"lambda"=c(x$distrpar, x$distrpar_se) )
  
  cat("Parameter estimates: \n")
  print(parall)  
  
}
