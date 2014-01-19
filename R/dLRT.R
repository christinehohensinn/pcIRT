dLRT <-
function(MPRMobj){
  
  eprm <- EPRM_red(MPRMobj$data)
  
  emp_Chi2 <- -2*(eprm$logLikelihood-MPRMobj$logLikelihood)
  df <- length(MPRMobj$estpar)-length(eprm$estpar)
  pvalue <- 1-pchisq(emp_Chi2,df)
  
  res <- list(emp_Chi2=emp_Chi2, df=df, pvalue=pvalue)
  class(res) <- "dLR"
  res
}
