#gamma function for RM

gamfunk <- function(eps.mat){
  gam.mat <- matrix(NA,ncol=length(eps.mat), nrow=length(eps.mat))
  gam.mat[1,]<- cumsum(eps.mat)
  for (i in 2:nrow(gam.mat)){
    for (j in 2:ncol(gam.mat)){
      if(j>=i){
        if(i==j){
          gam.mat[i,j] <- gam.mat[i-1,j-1]*eps.mat[j]
        } else {
          gam.mat[i,j] <- gam.mat[i, j-1] + gam.mat[i-1, j-1]*eps.mat[j]
        }
        
      }
    }
  }
  quot <- matrix(c(c(1,gam.mat[1:(nrow(gam.mat)-1),(length(eps.mat)-1)]) / gam.mat[,length(eps.mat)]), ncol=1)
  return(list(gammat=gam.mat[,length(eps.mat)], gam.quot=quot))
}  
