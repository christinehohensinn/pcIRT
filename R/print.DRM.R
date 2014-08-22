#' @rdname drm
#'
#' @export
#'
#' @method print DRM
#'
#' @param x object of class \code{DRM}

print.DRM <-
  function(object, ...){
    
    parall <- rbind(cbind("item estimates"=object$itempar, "SE"=object$itempar_se, "SE low"= object$itempar_se_low, "SE up"=object$itempar_se_up),"distrpar"=c(object$distrpar, object$distrpar_se, object$distrpar_se_low, object$distrpar_se_up) )
    
    cat("Parameter estimates: \n")
    print(parall)  
    
  }
