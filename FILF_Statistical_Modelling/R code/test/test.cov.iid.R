################################build function to test independent covariance matrix###########################

test.cov.iid <- function(data, numLongiPoints, nbs=1000, nb=10){
  
  source("calc.mean.R")
  source("select.knots.R")
  source("trunc.mat.R")
  source("ts.fit.R")
  source("calc.RA.R")
  source("calc.sigsq.R")
  source("irreg2mat.mod.R")
  source("p.bs.R")
  source("resample.iid.R")
  
  times<-seq(0,1,length.out=numLongiPoints)
  mu<-calc.mean(data,nb)

  data.demean <- data.frame(.value=data$.value-mu,.index=data$.index,.id=data$.id)
  b.fit<-ts.fit(data.demean,times=times,H=nb)
  C.alt<-trunc.mat(b.fit,b.fit$R,times) # Smooth alt cov
  #sigsq<- calc.sigsq(data.demean,C.alt,times) # error var
  data.sigsq <- data.frame(y=data$.value-mu,argvals=data$.index,subj=data$.id)
  sigsq<-(face.sparse(data.sigsq))$sigma2
  
  C.null<-matrix(0,nc=numLongiPoints,nr=numLongiPoints) # null cov
  Tn <- norm(C.alt-C.null,type='F')
  
  
  #bootstrap
  bs.stats<-c()
  bs.success=0
  while(bs.success<nbs){
    this.bs<-resample.iid(data,mu,sigsq)
    mu.bs<-calc.mean(this.bs,nb)
    data.demean.bs <- data.frame(.value=this.bs$.value-mu.bs,.index=this.bs$.index,.id=this.bs$.id)
    #fit.null.bs<-fitNull.iid(data.demean.bs) # null fit
    RbarA<-calc.RA(data.demean.bs) # R for alt
    #Rbar0.fit.bs<-calc.R0.exch(fit.null.bs,data.demean.bs)
    C.alt.bs<-trunc.mat(b.fit,RbarA,times) #slow here
    C.null.bs<-matrix(0,nc=numLongiPoints,nr=numLongiPoints)
    Tn.bs<-norm(C.alt.bs-C.null.bs,type='F')
    bs.stats<-c(bs.stats,Tn.bs) #save bs stats
    bs.success=bs.success+1
  }
  Tn.stats<-p.bs(Tn,unlist(bs.stats))
  
  return(list(C.alt=C.alt,C.null=C.null,Tn=Tn,p=Tn.stats$p,bs.approx=bs.stats,sigsq=sigsq))
}

