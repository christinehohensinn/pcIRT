MPRMl <-
  function(data, desmat, ldes,lp, start){
    
    call <- match.call()
    
    if(any(diff(as.numeric(names(table(data))),lag=1) != 1)){stop("categories level must be consecutive numbers")}
    
    if(is.data.frame(data)) {data <- as.matrix(data)}
    
    if(min(data) != 0){
      data <- data - min(data)
    }
    
    kateg.zahl <- length(table(data))
    item.zahl <- ncol(data)
    
    # margin vector groups of persons
    
    row.table <- apply(data+1,1,function(x) sprintf("%04d",(tabulate(x,nbins=kateg.zahl))))
    
    pat   <- apply(row.table,2, function(n) paste0(n, collapse=""))
    patt  <- table(pat)
    
    #first term (last category left out because these item parameters are 0)
    
    col.table <- apply(data+1, 2, function(s) tabulate(s,nbins=kateg.zahl))
    
    #designmatrix
    if(missing(desmat)){
      desmat <- designMPRM(data)
    } else {desmat <- as.matrix(desmat)} 
    
    #starting values
    #improved starting values for MPRM
    if(missing(start)){
      startval <- rep(0, (ncol(desmat)+max(lp)))
   } else {
      if (!is.numeric(start)) stop("Error: starting values are not numeric!")
      if (length(start) != ncol(desmat)) stop("Error: incorrect number of starting values!")}
    
    #pattern
    
    patmat <- t(xsimplex(kateg.zahl,item.zahl))
    patmat.o <- patmat[order(patmat[,kateg.zahl], decreasing=T),]
    
    patt.c <- apply(patmat.o,1,function(p) paste0(sprintf("%04d",p),collapse=""))
    
    cL <- function(para=startval,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o, patt.c=patt.c, patt=patt, desmat=desmat, ldes=ldes, lp=lp){
      
      fmat <- desmat %*% para[-c((length(para)-max(lp)+1):length(para))]
      for(k in seq_along(lp)){
        fmat[ldes!=0][k] <- fmat[ldes[ldes!=0][k]]*(para[c((length(para)-max(lp)+1):length(para))][lp[k]])
      }
      
      eps <- matrix(exp(fmat), nrow=kateg.zahl)     
      fir <- sum(col.table* log(eps))
      
      
      cf.g <- combfunc(kateg.zahl, item.zahl, eps.mat=t(eps), patmat.o)
      
      #search for gamma for each patt
      
      ord <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
      cfg.o <- log(cf.g$gammat[ord])
      sec <- sum(patt*cfg.o,na.rm=T)
      
      
      fir-sec
    }

res <- optim(startval, cL, kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o,patt.c=patt.c, patt=patt, desmat=desmat, ldes=ldes, lp=lp,method="BFGS", control=list(maxit=500, fnscale=-1), hessian=TRUE)
    
    estpar_se <- sqrt(diag(solve(res$hessian*(-1))))
    
    fmat <- desmat %*% res$par[-c((length(res$par)-max(lp)+1):length(res$par))]
    
    for(k in seq_along(lp)){
      fmat[ldes!=0][k] <- fmat[ldes[ldes!=0][k]]*(res$par[c((length(res$par)-max(lp)+1):length(res$par))][lp[k]])
    }
    
    itmat_se <- matrix(sqrt(diag(desmat %*% solve(res$hessian[-c((length(res$par)-max(lp)+1):length(res$par)),-c((length(res$par)-max(lp)+1):length(res$par))]*(-1)) %*% t(desmat))), nrow=kateg.zahl)
    
    itmat <- matrix(as.vector(fmat), nrow=kateg.zahl)
    
    if(!is.null(colnames(data))){
      colnames(itmat) <- paste("beta", colnames(data))
      colnames(itmat_se) <- paste("SE", colnames(data))
    } else {
      colnames(itmat) <- paste("beta item", 1:ncol(itmat))
      colnames(itmat_se) <- paste("SE item", 1:ncol(itmat))
    }
    
    rownames(itmat) <- paste("cat", 1:nrow(itmat))
    rownames(itmat_se) <- paste("cat", 1:nrow(itmat))

    linpar <- res$par[c((length(res$par)-max(lp)+1):length(res$par))]
    linpar_se <- estpar_se[c((length(res$par)-max(lp)+1):length(res$par))]
    
    res_all <- list(data=data, design=desmat, logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itmat*(-1), itempar_se=itmat_se, linpar=linpar, 
linpar_se=linpar_se, hessian=res$hessian, convergence=res$convergence, fun_calls=res$counts, call=call)
    class(res_all) <- "MPRMl"
    res_all
}

