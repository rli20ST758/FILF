################################build function to test exchangeable covariance matrix###########################

test.cov.exch <- function(data, numLongiPoints, nbs=1000, nb=10){
  
  source("calc.mean.R")
  source("select.knots.R")
  source("trunc.mat.R")
  source("ts.fit.R")
  source("calc.RA.R")
  source("calc.sigsq.R")
  source("irreg2mat.mod.R")
  source("p.bs.R")
  source("exch.R")
  
  times<-seq(0,1,length.out=numLongiPoints)
  mu <- calc.mean(data,nb)
  
  data.demean <- data.frame(.value=data$.value-mu,.index=data$.index,.id=data$.id)
  fit.null <- fitNull.exch(data.demean) # null fit
  if('try-error' %in% class(fit.null)){
    #issue with null fit
    next
  }
  
  sigsq0 <- as.numeric(VarCorr(fit.null)[1,1])
  
  b.fit<-ts.fit(data.demean,times=times,H=nb)
  C.alt<-trunc.mat(b.fit,b.fit$R,times) # Smooth alt cov
  data.sigsq <- data.frame(y=data$.value-mu,argvals=data$.index,subj=data$.id)
  sigsq<-(face.sparse(data.sigsq))$sigma2
  
  Rbar0.fit<-calc.R0.exch(fit.null,data)
  C.null<-trunc.mat(b.fit,Rbar0.fit,times) # Smooth null cov
  Tn <- norm(C.alt-C.null,type='F')
  
  #bootstrap
  bs.stats<-c()
  bs.success=0
  while(bs.success<nbs){
    this.bs<-resample.exch(data, mu, fit.null, sigsq)
    mu.bs <- calc.mean(this.bs,nb)
    data.demean.bs <- data.frame(.value=this.bs$.value-mu.bs,.index=this.bs$.index,.id=this.bs$.id)
    fit.null.bs<-fitNull.exch(data.demean.bs) # null fit
    if('try-error' %in% class(fit.null)){
      #issue with null fit
      next
    }
    RbarA<-calc.RA(data.demean.bs) # R for alt
    Rbar0.fit.bs<-calc.R0.exch(fit.null.bs,data.demean.bs)
    C.alt.bs<-trunc.mat(b.fit,RbarA,times) #slow here
    C.null.bs<-trunc.mat(b.fit,Rbar0.fit.bs,times) #slow here
    Tn.bs<-norm(C.alt.bs-C.null.bs,type='F')
    bs.stats<-c(bs.stats,Tn.bs) #save bs stats
    bs.success=bs.success+1
  }
  Tn.stats<-p.bs(Tn,unlist(bs.stats))
  
  return(list(C.alt=C.alt,C.null=C.null,Tn=Tn,p=Tn.stats$p,bs.approx=bs.stats,sigsq=sigsq,sigsq0=sigsq0))
}

