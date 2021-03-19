############################
# Written by: Stephanie Chen (stchen3@ncsu.edu)
# Purpose: Calculates the error variance based on a covariance estimate (input: C)
# Updated: Aug 4, 2018

calc.sigsq <- function(data,C,times){
  ymat <- irreg2mat.mod(data,times) # dense mat
  nsubj <- nrow(ymat)
  
  Rii <- t(sapply(1:nsubj,function(x) ymat[x,]^2)) # Yij * Yij
  diag.cov <- colMeans(Rii, na.rm=TRUE) # diag of sample cov
  diag.cov[is.nan(diag.cov)]<-0
  ind.low<-which(abs(times-quantile(times,0.25))==min(abs(times-quantile(times,0.25))))
  ind.high<-which(abs(times-quantile(times,0.75))==min(abs(times-quantile(times,0.75))))
  #ind.low=1
  #ind.high=41
  return(2*max(0,sintegral(times[ind.low:ind.high],(diag.cov-diag(C))[ind.low:ind.high])$value))
}
