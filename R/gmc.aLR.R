gmc.aLR <-
function(object, ...){
  if(length(object$itempar) > 2) {stop("there are more than two subsamples! Only two subsamples can be plotted")}
  
  labs <- as.vector(t(outer(1:ncol(object$itempar[[1]]), paste0("-", 1:(nrow(object$itempar[[1]])-1)), FUN="paste0")))
  
  plot(object$itempar[[1]][-nrow(object$itempar[[1]]),], object$itempar[[2]][-nrow(object$itempar[[2]]),], type="p", pch=20, xlab="item par. group 1", ylab="item par. group 2", xlim=c(-ceiling(max(abs(unlist(object$itempar)))), ceiling(max(abs(unlist(object$itempar))))), ylim=c(-ceiling(max(abs(unlist(object$itempar)))),ceiling(max(abs(unlist(object$itempar))))))
  abline(0,1)
  text(object$itempar[[1]], object$itempar[[2]], labels=labs, pos=1)
  
}
