#' Estimation of Continuous Rating Scale Model (Mueller, 1987)
#'
#' Estimation of the Rating Scale Model for continuous data by Mueller (1987).
#'
#' \deqn{P_{vi}(a \leq X \leq b) = \frac{\int_a^b exp[x \mu + x(2c-x) \theta]
#' dx}{\int_{c-\frac{d}{2}}^{c+\frac{d}{2}} exp[t \mu + t(2c-t) \theta] dt}}
#'
#' Parameters are estimated either by a pairwise conditional likelihood estimation (as pseudo-likelihood approach) or  a pairwise algorithm.
#'
#' The parameters of the Continuous Rating Scale Model are estimated by a cml approach. a
#' pairwise algorithm (the algorithm is described in detail in Mueller, 1999)
#' using Newton-Raphson iterations for optimizing.  Alternatively item parameters can be estimated by conditional maximum likelihood estimation. Because of the time-consuming estimation of high-dimensional integrals, a pairwise conditional likelihood estimation is used in this case. For both estimation methods no assumption on the person parameter distribution are necessary. Note, that "pcml" as estimation method is considerably slower than the "pair" method.
#'
#' @aliases CRSM summary.CRSM print.CRSM
#' @param data Data matrix or data frame; rows represent observations
#' (persons), columns represent the items.
#' @param min The minimum value of the response scale (on which the data are
#' based).
#' @param max The maximum value of the response scale (on which the data are
#' based).
#' @param method Estimation method is either maximization of the pseudo-conditional-maximum-likelihood ("pcml")
#' or a pairwise algorithm("pair"), default is "pair"
#' @param start Starting values for parameter estimation. If missing, a vector
#' of 0 is used as starting values.
#' @param conv Convergence criterium for the pairwise algorithm ("pair")
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
#' data(example1)
#' res_crsm <- CRSM(example1, min=0, max=1)
#'
#' summary(res_crsm)
#'#'
#' @export CRSM
#' @rdname crsm
CRSM <-
function(data, min, max, method="pair", start, conv=0.0001){

call <- match.call()

if(is.data.frame(data)) {data <- as.matrix(data)}
if(missing(min)){stop("Error: enter lowest possible value of items")}
if(missing(min)){stop("Error: enter highest possible value of items")}

# bring data in interval (0,1)
data_p <- (data - min)/(max-min)

#check if few oberservations for item / person
if (any(apply(data_p, 2, function(dc) sum(dc != 0 & dc != 1)) %in% c(0,1)))
{cat("warning: there is at least 1 item with only 1 observation containing no extreme rawscore.
     This will probably cause estimation problems")}
if (any(apply(data_p, 1, function(dc) sum(dc != 0 & dc != 1)) %in% c(0,1)))
{cat("warning: there is at least 1 person with only 1 observation containing no extreme rawscore.
     This will probably cause estimation problems")}

if(missing(start)){
  para1     <- rep(0, 2)
} else {
  para1     <- start
}

combis <- combn(ncol(data_p),2)

if(method=="pcml"){
  #pseudo CML

  if(nrow(data) > 1000 | ncol(data) > 15){cat("warning: pcml estimation method can take some time for a data matrix of this size. Note that the pair method is faster.")}

  cloglik_crsm <- function(para,dataset){
    si <- colSums(dataset)[-1]   #Itemscores
    omega <- sum(dataset^2)

    zahler <- as.vector(-si%*%para[1:(ncol(dataset)-1)]-omega*para[ncol(dataset)])

    hilfsvec <- ncol(dataset):2

    nenn <- apply(dataset,1,function(roh){
      rvec <- cumsum(rev(roh))
      upperb <- rvec[length(roh)] - c(0,rvec[-c(length(roh),(length(roh)-1))])-(hilfsvec-1)*0
      lowerb <- rvec[length(roh)] - c(0,rvec[-c(length(roh),(length(roh)-1))])-(hilfsvec-1)*1
      lowerb2 <- rev(ifelse(lowerb<0, 0,lowerb))
      upperb2 <- rev(ifelse(upperb>1, 1,upperb))

      integrand <- function(arg){
        tval <- arg[1:(ncol(dataset)-1)]
        nennteil <- exp(-tval%*%para[1:(ncol(dataset)-1)]- (rvec[ncol(dataset)]-sum(tval))^2*para[ncol(dataset)] - para[ncol(dataset)]*sum(tval^2))
        return(nennteil)
      }

      stest <- function(t, paraI){exp(-t*paraI[1]- (rvec[2]- t)^2*paraI[2] - paraI[2]*(t^2))}

      integrate(stest, paraI=para, lower=lowerb2, upper=upperb2, stop.on.error=F)$value
    }
    )

    #nenner <- sapply(nenn, function(v) v$integral)

    zahler - sum(log(nenn))

  }
  parlist <- lapply(1:ncol(combis), function(m) {
      zwi <- data_p[,c(combis[1,m],combis[2,m])]
      zwi2 <- zwi[rowSums(zwi)!=0 & rowSums(zwi)!=2,, drop=FALSE] #delete extreme scores
      optim(para1, fn=cloglik_crsm, dataset=zwi2, control=list(fnscale=-1, maxit=200), method="BFGS", hessian=TRUE)
  })

  beta       <- sapply(parlist, function(l) l$par[1])*(-1)
  lambda.all <- sapply(parlist, function(la) la$par[2])
  iterations <- sapply(parlist, function(it) it$counts)

} else if(method=="pair"){ #pairwise

iterations <- vector(mode="numeric",length=ncol(combis))

S0n <- function(t, paraI)      {exp(-(t/2)*paraI[1] - (t^2/2)*paraI[2])}
S1n <- function(t, paraI) {t*   exp(-(t/2)*paraI[1] - (t^2/2)*paraI[2])}
S2n <- function(t, paraI) {t^2* exp(-(t/2)*paraI[1] - (t^2/2)*paraI[2])}
S3n <- function(t, paraI) {t^3* exp(-(t/2)*paraI[1] - (t^2/2)*paraI[2])}
S4n <- function(t, paraI) {t^4* exp(-(t/2)*paraI[1] - (t^2/2)*paraI[2])}



parlist <- lapply(1:ncol(combis), function(m) {
  zwi <- data_p[,c(combis[1,m],combis[2,m])]
  zwi2 <- zwi[rowSums(zwi)!=0 & rowSums(zwi)!=2,, drop=FALSE] #delete extreme scores

  uij <- zwi2[,1]-zwi2[,2]
  sij <- sum(uij)

  bound       <- round(1-abs(rowSums(zwi2)-1),6) #interval width
  bound.v     <- sort(unique(bound))
  bound.n     <- unname(table(bound))

  iter        <- 0
  para        <- NA

  while( is.na(para) || max(abs(para1-para)) > conv){
    para <- para1

    S0n.i <- sapply(bound.v, function(b) integrate(S0n,paraI=para,lower=-b, upper=b, stop.on.error=F)$value)
    S1n.i <- sapply(bound.v, function(b) integrate(S1n,paraI=para,lower=-b, upper=b, stop.on.error=F)$value)
    S2n.i <- sapply(bound.v, function(b) integrate(S2n,paraI=para,lower=-b, upper=b, stop.on.error=F)$value)
    S3n.i <- sapply(bound.v, function(b) integrate(S3n,paraI=para,lower=-b, upper=b, stop.on.error=F)$value)
    S4n.i <- sapply(bound.v, function(b) integrate(S4n,paraI=para,lower=-b, upper=b, stop.on.error=F)$value)


    abl1    <- vector(length=2)
    abl1[1] <- -sij/2 + 0.5*sum((S1n.i/S0n.i)*bound.n)
    abl1[2] <- -sum(uij^2)/2 + 0.5*sum((S2n.i/S0n.i)*bound.n)

    abl2        <- matrix(NA,ncol=2,nrow=2)
    abl2[1,1]   <- -0.25 * sum(((S2n.i/S0n.i)-(S1n.i/S0n.i)^2)*bound.n)
    abl2[1,2]   <- -0.25 * sum(((S3n.i/S0n.i)-((S1n.i*S2n.i)/S0n.i^2))*bound.n)
    abl2[2,1]   <- abl2[1,2]
    abl2[2,2]   <- -0.25 * sum(((S4n.i/S0n.i)-(S2n.i/S0n.i)^2)*bound.n)

    para1    <- as.vector(para - abl1%*%solve(abl2))
    iter <- iter+1
  }
  return(list(para1=para1, iterations=iter, hessian=abl2))
}
)


beta       <- sapply(parlist, function(l) l$para1[1])
lambda.all <- sapply(parlist, function(la) la$para1[2])
iterations <- sapply(parlist, function(it) it$iterations)
} else {
  stop("error: enter correct estimation method or leave it empty for using the default")
}

#calculate item parameters

lambda     <- mean(lambda.all)

betas <- sapply(1:ncol(data_p), function(be) {
  sum(beta[combis[1,]==be],beta[combis[2,]==be]*(-1))/ncol(data_p)
}
)
betas <- betas - mean(betas)
names(betas) <- paste(rep("beta ", length(betas)), colnames(data))

#standard errors for items
varDiff.item <- sapply(parlist, function(s) diag(solve(s$hessian)*(-1))[1])
sdDiff.item <- sapply(parlist, function(s) sqrt(diag(solve(s$hessian))*(-1))[1])
varDiff.distr <- sapply(parlist, function(s) diag(solve(s$hessian)*(-1))[2])
sdDiff.distr <- sapply(parlist, function(s) sqrt(diag(solve(s$hessian)*(-1))[2]))

se.item.up <- sapply(1:ncol(data_p), function(di) {
  sqrt((sum(sdDiff.item[combis[1,]==di],sdDiff.item[combis[2,]==di]))^2/(ncol(data_p))^2)
}
)
names(se.item.up) <- paste(rep("SE(beta) up ", length(betas)), colnames(data))

se.item.low <- sapply(1:ncol(data_p), function(di) {
  sqrt(sum(varDiff.item[combis[1,]==di],varDiff.item[combis[2,]==di])/(ncol(data_p))^2)
}
)
names(se.item.low) <- paste(rep("SE(beta) low ", length(betas)), colnames(data))

se.item.mean <- (se.item.up+se.item.low)/2
names(se.item.mean) <- paste(rep("SE(beta) ", length(betas)), colnames(data))

se.distr.low <- sqrt(sum(varDiff.distr)/(ncol(combis))^2)
names(se.distr.low) <- "SE(lambda) low"

indmat <- apply(combis, 2, function(co){
  colSums(matrix(combis %in% co,nrow=2))
})
indmat[indmat==2] <- 1
indmat.n <- t(indmat*sdDiff.distr)

se.distr.up <- sqrt(sum(indmat.n*sdDiff.distr)/ncol(combis)^2)
names(se.distr.up) <- "SE(lambda) up"

se.distr.mean <- mean(c(se.distr.low, se.distr.up))
names(se.distr.mean) <- "SE(lambda)"




res_all   <- list(data=data, data_p=data_p, itempar=betas,itempar_se_low=se.item.low, itempar_se_up=se.item.up, itempar_se=se.item.mean, distrpar=lambda, distrpar_se_low=se.distr.low, distrpar_se_up=se.distr.up,  distrpar_se=se.distr.mean, iterations=iterations, call=call)

class(res_all) <- "CRSM"

res_all

}
