#########################fit null hypothesis covariance matrix coefficience##########################
fitNull.exch<-function(data){
  try(lme(.value ~ .id, random = ~ 1|.id, method="ML", data=data))
  #try(rlmer(.value ~ 1 + (1|.id), method="DASvar", data=data))
  #VarCorr(temp)$.id[1]
  #try(lmer(.value~ (1|.id),data=data, REML = FALSE)
  #sigsq0 <- as.data.frame(VarCorr(temp))[1,4]
}


####################################################################################################


###########################generate null covariance matrix##########################################
calc.R0.exch<-function(fit.null,data){
  subj <- data$.id
  usubj <- unique(subj)
  nsubj <- length(usubj)
  mi <- unlist(lapply(1:nsubj,function(i) length(which(subj==usubj[i]))))
  m <- sum(mi*(mi-1))

  sigsq0 <- as.numeric(VarCorr(fit.null)[1,1])
  Rbar0 <- rep(sigsq0,m)
  return(Rbar0)
}
#####################################################################################################


##############################resample function for bootstrap#######################################
resample.exch<-function(data,mu,fit.null,sigsq){
  subj <- data$.id
  usubj <- unique(subj)
  nsubj<-length(unique(data$.id))
  mi <- unlist(lapply(1:nsubj,function(i) length(which(subj==usubj[i]))))
  
  sigsq0 <- as.numeric(VarCorr(fit.null)[1,1])
  xi_coef <- rnorm(n = nsubj, mean = 0, sd = sqrt(sigsq0))
  xi_i <- rep(xi_coef, mi)
  xi_ij <- rnorm(sum(mi), mean=0, sd=sqrt(sigsq))
  
  #mean.null + null.variance
  data.bs<-data.frame(.value=mu+xi_i+xi_ij,.index=data$.index,.id=data$.id)
  return(data.bs)
}
########################################################################################################
