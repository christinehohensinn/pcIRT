#' Estimation of dichotomous logistic Rasch model (Rasch, 1960)
#' 
#' Estimation of dichotomous logistic Rasch model (Rasch, 1960).
#' 
#' Parameters are estimated by CML.
#' 
#' .
#' 
#' @param data Data matrix or data frame; rows represent observations
#' (persons), columns represent the items.
#' @return \item{data}{data matrix according to the input} \item{data_p}{data
#' matrix with data transformed to a response interval between 0 and 1}
#' \item{itempar}{estimated item parameters} \item{itempar_se_low}{estimated
#' lower boundary for standard errors of estimated item parameters}
#' \item{itempar_se_up}{estimated upper boundary for standard errors of
#' estimated item parameters} \item{itempar_se}{estimated mean standard errors
#' of estimated item parameters} \item{distrpar}{estimated distribution
#' parameter} \item{distrpar_se_low}{estimated lower boundary for standard
#' errors of estimated distribution parameter} \item{distrpar_se_up}{estimated
#' upper boundary for standard errors of estimated distribution parameter}
#' \item{itempar_se}{estimated mean standard errors of estimated distribution
#' parameter} \item{iterations}{Number of Newton-Raphson iterations for each
#' item pair} \item{call}{call of the CRSM function}
#' @author Christine Hohensinn
#' @references Mueller, H. (1987). A Rasch model for continuous ratings.
#' Psychometrika, 52, 165-181.
#' 
#' Mueller, H. (1999). Probabilistische Testmodelle fuer diskrete und
#' kontinuierliche Ratingskalen. [Probabilistic models for discrete and
#' continuous rating scales]. Bern: Huber.
#' @keywords continuous rating scale model
#' @examples
#' 
#' #estimate CRSM item parameters
#' #data(example1)
#' #res_crsm <- CRSM(example1, min=0, max=1)
#' 
#' #summary(res_crsm)
#' 
#' 
#' @export DRM
DRM <- 
  function(data, desmat, start){

data <- as.matrix(data)
rv <- rowSums(data)
si <- colSums(data)
dat.r <- data[!(rv %in% c(0,ncol(data))),]

rv.r <- rowSums(dat.r)
si.r <- colSums(dat.r) 
fr.r <- tabulate(rv.r, nbins=(ncol(data)-1))

if(missing(start)){
  para <- rep(0,(ncol(data)-1))  
} else if(start == "improved"){
  para <-(nrow(dat.r)-si.r[1:(ncol(dat.r)-1)])/si.r[1:(ncol(dat.r)-1)]-(sum(log((nrow(dat.r)-si.r)/si.r))/length(si.r))
} else if(start != "improved" & is.numeric(start)==FALSE){
  stop("start values have to be numeric!")
} else {
  para <- start
}

if(missing(desmat)){
  desmat <- rbind(-1,diag(1,nrow=length(para)))
}

cmlLoglik <- function(sigs, si, roh, desmat){
  sigs2 <- desmat %*%sigs
  term1 <- as.vector(si%*%sigs2)
  epsi  <- exp(sigs2)
  fgmat <- matrix(0,ncol=length(epsi), nrow=length(epsi))
  fgmat[1,]<- cumsum(epsi)
  gam <- gamfunk(epsi, fgmat)[1:(ncol(data)-1)]
  term2 <- as.vector(roh%*%log(gam))
  term1-term2
}

res <- optim(para, fn=cmlLoglik, si=si.r, roh=fr.r,desmat=desmat, control=list(fnscale=-1, maxit=1000), hessian=TRUE, method="L-BFGS-B")

# normalize item parameters
itempar <- as.vector(desmat %*% res$par)

#SE
estpar_se <- sqrt(diag(solve(res$hessian*(-1))))

itempar_se <- sqrt(diag(desmat %*% solve(res$hessian*(-1)) %*% t(desmat)))
                   

res_all <- list(data=data, logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itempar*(-1), itempar_se=itempar_se, hessian=res$hessian*(-1), convergence=res$convergence, fun_calls=res$counts, call=call)
class(res_all) <- "DRM"
res_all
}


