#########################fit null hypothesis covariance matrix coefficience##########################
fitNull.iid<-function(data){
  
  sigsq0 <- var(data$.value)
  #sigsq1 <- as.numeric(VarCorr(temp)[2,1])
  return(sigsq0)
}
####################################################################################################


##############################resample function for bootstrap#######################################
resample.iid<-function(data,mu,sigsq){
  subj <- data$.id
  usubj <- unique(subj)
  nsubj<-length(unique(data$.id))
  mi <- unlist(lapply(1:nsubj,function(i) length(which(subj==usubj[i]))))
  
  xi_ij <- rnorm(sum(mi), mean=0, sd=sqrt(sigsq))
  
  #mean.null + null.variance
  data.bs<-data.frame(.value=mu+xi_ij,.index=data$.index,.id=data$.id)
  return(data.bs)
}