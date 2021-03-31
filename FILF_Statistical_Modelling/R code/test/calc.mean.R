#######################
# Written by: Stephanie Chen (stchen3@ncsu.edu)
# Purpose: Calculates smooth mean using 10 thin plate regression splines
# Updated: Aug 4, 2018

calc.mean<-function(data,nb=10){ #smooth mean
  gam0<-gam(as.vector(data$.value)~s(data$.index,k=nb))
  #gam0<-gam(as.vector(data$.value)~s(data$.index,k=10))
  return(gam0$fitted.values)
}
