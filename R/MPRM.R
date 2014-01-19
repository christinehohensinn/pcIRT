MPRM <-
function(data, desmat, start){

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

#Startwerte

#improved starting values for MPRM
if(missing(start)){
  startval <- rep(0, ncol(desmat))
  } else if(is.numeric(start)){
  startval <- start
  } else if(start=="improved"){
  diff <- sweep(col.table,2,(col.table[kateg.zahl,]), FUN="/")
  quotcs <- sweep((col.table)^(-1),2, col.table[kateg.zahl,], FUN="*")
  zweit <- (apply(quotcs, 2, function(ss){prod(ss)}))^(1/item.zahl)
  startval1 <- log(as.vector(t(t(diff)*zweit)[-kateg.zahl,-item.zahl]))
  startval <- as.vector(t(apply(matrix(as.vector(desmat %*% startval1), nrow=kateg.zahl)[-kateg.zahl, -item.zahl],1, function(st) scale(st, scale=F))))
  } else {
  #if (!is.numeric(start)) stop("Error: starting values are not numeric!")
  if (length(start) != ncol(desmat)) stop("Error: incorrect number of starting values!")
}
#erstellen der Matrix mit allen Randsummenpattern for combinatoric function

patmat <- t(xsimplex(kateg.zahl,item.zahl))
patmat.o <- patmat[order(patmat[,kateg.zahl], decreasing=T),]

patt.c <- apply(patmat.o,1,function(p) paste0(sprintf("%04d",p),collapse=""))

cL <- function(para=startval,kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o, patt.c=patt.c, patt=patt, desmat=desmat){
  eps <- matrix(exp(desmat %*% para), nrow=kateg.zahl)
  
  fir <- sum(col.table* log(eps))
  
  if(kateg.zahl > 2){
    cf.g <- combfunc(kateg.zahl, item.zahl, eps.mat=t(eps), patmat.o)
    
    ord <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
    cfg.o <- log(cf.g$gammat[ord])
    sec <- sum(patt*cfg.o,na.rm=T)
    
  } else {
    cf.g <- gamfunk(eps.mat=eps[1,])
    sec <- sum(patt[-1] * log(cf.g$gammat), na.rm=T)
  }
  #search for gamma for each patt

  fir-sec
}

der1 <- function(para=startval, col.table=col.table, kateg.zahl=kateg.zahl, item.zahl=item.zahl, patmat.o=patmat.o, patt=patt, patt.c=patt.c, desmat=desmat){
  eps <- matrix(exp(as.vector(desmat %*% para)), nrow=kateg.zahl)  

    eps.f <- list(eps)
  
    for (e in seq_len(item.zahl-1)){
      eps.f[[e+1]] <- cbind(eps.f[[e]][,2:item.zahl],eps.f[[e]][,1])
    }
  if(kateg.zahl > 2){
    cf.all <- lapply(eps.f, function(gg) {combfunc(kateg.zahl, item.zahl, eps.mat=t(gg), patmat.o)})
  
  #Vektor raussuchen wo die auftretenden pattern in patt.c vorkommen
  
    ord2 <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
    cf.o <- lapply(cf.all, function(ln) ln$gam.quot[ord2,])
  
    cf.oNR <- sapply(1:item.zahl, function(l2) {colSums(t(t(cf.o[[l2]]))*as.vector(patt), na.rm=TRUE)})
  
    cf.oNR2 <- cf.oNR*eps[-kateg.zahl,]
  } else {
    
    for (e in seq_len(item.zahl)){
      eps.f[[e+1]] <- cbind(eps.f[[e]][,2:item.zahl],eps.f[[e]][,1])
    }    
    eps.f[[1]] <- NULL

    cf.all <- lapply(eps.f, function(gg) {gamfunk(eps.mat=gg[1,])})
    
    #Vektor raussuchen wo die auftretenden pattern in patt.c vorkommen
    
    #ord2 <- sapply(names(patt), function(hpat) which(patt.c %in% hpat))
    cf.o <- lapply(cf.all, function(ln) ln$gam.quot)
    
    cf.oNR <- sapply(1:item.zahl, function(l2) {colSums(t(t(cf.o[[l2]]))*as.vector(patt[-c(1)]), na.rm=TRUE)})
    
    cf.oNR2 <- cf.oNR*eps[1,]
    
  }
  #von früher ohne Designmatrix
  #as.vector((col.table[-kateg.zahl,-item.zahl] - col.table[-kateg.zahl,item.zahl]) - (cf.oNR2[,-item.zahl] - cf.oNR2[,item.zahl]))
 
  t(desmat[-seq(kateg.zahl,item.zahl*kateg.zahl, by=kateg.zahl),]) %*% as.vector(col.table[-kateg.zahl,] - cf.oNR2)
}

res <- optim(startval, cL,gr=der1, kateg.zahl=kateg.zahl, item.zahl=item.zahl,col.table=col.table, patmat.o=patmat.o,patt.c=patt.c, patt=patt, desmat=desmat, method="BFGS", control=list(maxit=500, fnscale=-1), hessian=TRUE)

estpar_se <- sqrt(diag(solve(res$hessian*(-1))))

itmat <- matrix(as.vector(desmat %*% res$par), nrow=kateg.zahl)
itmat_se <- matrix(sqrt(diag(desmat %*% solve(res$hessian*(-1)) %*% t(desmat))), nrow=kateg.zahl)

if(!is.null(colnames(data))){
  colnames(itmat) <- paste("beta", colnames(data))
  colnames(itmat_se) <- paste("SE", colnames(data))
  } else {
  colnames(itmat) <- paste("beta item", 1:ncol(itmat))
  colnames(itmat_se) <- paste("SE item", 1:ncol(itmat))
}
   
rownames(itmat) <- paste("cat", 1:nrow(itmat))
rownames(itmat_se) <- paste("cat", 1:nrow(itmat))


res_all <- list(data=data, design=desmat, logLikelihood=res$value, estpar=res$par, estpar_se=estpar_se, itempar=itmat*(-1), itempar_se=itmat_se, hessian=res$hessian, convergence=res$convergence, fun_calls=res$counts, call=call)
class(res_all) <- "MPRM"
res_all
}
