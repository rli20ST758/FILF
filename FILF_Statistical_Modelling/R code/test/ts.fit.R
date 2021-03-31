############################
# Written by: Stephanie Chen (stchen3@ncsu.edu)
# Purpose: Tensor-spline fit (using cubic B-splines) to data for smooth null and alternative covariance
# Updated: Aug 4, 2018

ts.fit <- function(data,times,H=10){
  p <- 3 ##cubic splines: degrees
  #knots <- select.knots(seq(-1,1,by=0.01),H-p)
  knots <- select.knots(seq(0,1,by=0.01),H-p)
  
  y <- data$.value
  t <- data$.index
  subj <- data$.id
  usubj <- unique(subj)
  n <- length(usubj)
  
  B <- matrix(,nrow=0,ncol=H^2)
  Time <- matrix(,nrow=0,ncol=2)
  Y <- c()
  R <- c()
  
  for(i in 1:n){
    index <- which(subj==usubj[i])
    mi <- length(index)
    ti <- t[index]
    yi <- y[index]
    
    Bi <- spline.des(knots=knots, x= ti, ord = p+1,outer.ok =TRUE)$design
    Btildei <- Bi%x%Bi
    Ji <- matrix(1,mi,mi)
    diag(Ji) <- 0
    di <- c(Ji)
    Bbari <- Btildei[which(di==1),]
    T1i <- (ti%x%rep(1,mi))[which(di==1)]
    T2i <- (rep(1,mi)%x%ti)[which(di==1)]
    Ri <- (yi%x%yi)[which(di==1)]
    
    B <- rbind(B,Bbari)
    Time <- rbind(Time,cbind(T1i,T2i))
    R <- c(R,Ri)
  }
  
  t0 <- times ##fine grid
  Bstar <- spline.des(knots=knots, x= t0, ord = 4,outer.ok =TRUE)$design
  temp1 <- t(Bstar)%*%Bstar
  Eigen1 <- eigen(temp1)
  A1 <- Eigen1$vectors%*%sqrt(diag(Eigen1$values))%*%t(Eigen1$vectors)
  temp2 <- t(B)%*%B
  Eigen2 <- eigen(temp2)
  BtB.inv <- Eigen2$vectors %*% tcrossprod(diag(1/Eigen2$values),Eigen2$vectors)
  
  return(list(B=B,Bstar=Bstar,Time=Time,R=R,BtB.inv=BtB.inv,Bstar.tensor=Bstar %x% Bstar))
}
