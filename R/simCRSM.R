simCRSM <-
function(itempar, distr, perspar, mid=0.5, len=1, seed=NULL){
  
mat.intercept <- outer(perspar, itempar, "-")  
  
forint <- as.vector(mat.intercept)
  

funk.n <- function(st, int, mid, distr){exp(st*int+ st*(2*mid-st)*distr)} 

nenn  <- sapply(forint, function(i) {
  integrate(funk.n, int=i, mid=mid, distr=distr, lower=(mid-0.5*len), upper=(mid+0.5*len))$value
})

interv <- seq(mid-0.5*len,mid+0.5*len,length.out=500)  
zahl.int <- lapply(forint, function(o) {
    sapply(1:(length(interv)-1), function(j) {
      integrate(funk.n, int=o, mid=mid, distr=distr, lower=interv[j], upper=interv[j+1])$value
  })
})

pval <- lapply(1:length(forint), function(y) zahl.int[[y]]/nenn[y])

if (!is.null(seed)) {set.seed(seed)}

p.response.sampled <- sapply(1:length(forint), function(s) sample(1:(length(interv)-1), size=1, prob=pval[[s]]))
p.response <- interv[p.response.sampled]

datmat <- matrix(p.response, ncol=length(itempar))
list(datmat=datmat, true_itempar=itempar, true_distrpar=distr, true_perspar=perspar)
}
