simMPRM <-
function(itempar, persons=500, seed=NULL){

#normalization of person parameters
kateg <- nrow(itempar)

ppar <- rnorm(persons*kateg, 0,1)

lauf <- seq(1,length(ppar), by=kateg)
for (l in lauf) {
  ppar[l:(l+kateg-1)] <- ppar[l:(l+kateg-1)] - ppar[l+kateg-1]    
  }
items <- ncol(itempar)
itmat <- matrix(rep(itempar*(-1), persons), ncol=items*kateg, byrow = TRUE)

persmat <- matrix(ppar, nrow=persons, byrow=T)

persmat.ext <- matrix(rep(persmat,items), ncol=items*kateg, byrow=FALSE)

zahler <- exp(persmat.ext+itmat)

nenn <- sapply(seq(1,kateg*items, by=kateg), function(y) rowSums(zahler[,y:(y+kateg-1)], na.rm=TRUE))
nenn.est <- matrix(rep(nenn, each=kateg),ncol=items*kateg,byrow=TRUE)


prob.mat <- mapply(function(z,n) zahler[,z]/nenn[,n], z=1:(kateg*items), n=rep(1:items, each=kateg))

if (!is.null(seed)) {set.seed(seed)}
   
datmat <- t(apply(prob.mat, 1, function(ma){
         sapply(seq(1,kateg*items, by=kateg), function(ma2) {sample(0:(kateg-1),size=1,prob=ma[ma2:(ma2+kateg-1)])})
      }))

return(list(datmat=datmat, true_itempar=itempar, true_perspar=ppar))
}
