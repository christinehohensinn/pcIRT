weight_test <-
function(MPRMobj, score_param){
  
  call <- match.call()
  
  if(length(score_param) != (length(table(MPRMobj$data))-2)){stop("Error: wrong number of score parameters!")}
    
  #LR-Test for unidimensionality of item parameters
  
  ep_res  <- EPRM_red(MPRMobj$data)
  ep_resS <- EPRM_red(MPRMobj$data, score_par=score_param)
  
    chi2 <- -2*(ep_resS$logLikelihood-ep_res$logLikelihood)
    df <- length(ep_res$estpar)-length(ep_resS$estpar)
    pvalue <- 1-pchisq(chi2,df)
    
  res <- list(emp_Chi2=chi2, df=df, pval=pvalue, unconstrLoglikelihood=ep_res$logLikelihood, constrLogLikelihood=ep_resS$logLikelihood, unconstrNrPar=length(ep_res$estpar),constrNrPar=length(ep_resS$estpar), unconstrItempar=ep_res$itempar*(-1), constrItempar=ep_resS$itempar*(-1), unconstrScoreParameter=ep_res$score_par)
  class(res) <- "wt"
  res  
}
