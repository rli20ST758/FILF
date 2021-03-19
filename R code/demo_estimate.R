
library(mgcv)
library(refund)
library(face)
library(MASS)
library(ggplot2)

#setwd("~/Dowdonw/code")
source("Functions_100217.R")
source("fgee.R")
source("fgee_init.R")
source("fgee_initBp.R")


set.seed(1)
    
numFunctPoints = 101
Nsubj = 100
numLongiPoints = 41
    
data.unsp <- GenerateData(Nsubj=Nsubj, numFunctPoints = numFunctPoints, 
                          min_visit=3, max_visit=6, numLongiPoints = numLongiPoints, 
                          sigma_sq = 1.5, sigma_z11 = 3, sigma_z12 = 1.5, 
                          sigma_z21 = 2, sigma_z22 = 1, corstr = "unspecified")
    

########################################
#Method 1: initial estimate
########################################
    
start_time <- Sys.time()
init.estim <- fgee_init(formula = Y ~ X1 + X2, 
                        Y=data.unsp$data$Y, Cov=data.unsp$data$Cov, 
                        s=data.unsp$data$funcArg)
end_time <- Sys.time()
time.init <- difftime(end_time, start_time, units = "secs")  
    

Y.mean.init.ise <- mean.ise(init.estim$Y.mean.init, data.unsp$mean, numFunctPoints)
beta.init.acc <- beta.acc(beta.hat=init.estim$beta.init,beta.true=data.unsp$beta,numFunctPoints=numFunctPoints)
confidence.init <- eval.beta.ci(beta.hat=init.estim$beta.init,beta.se=init.estim$beta.init.se,beta.true=data.unsp$beta,numFunctPoints=numFunctPoints)
    
init <- c(time.init, Y.mean.init.ise, beta.init.acc, confidence.init$length.ave, confidence.init$coverage.rate)
    

###########################################################################
#Method 2: initial bootstrap estimate to get bootstrap confidence interval
###########################################################################

start_time <- Sys.time()
initBp.estim <- fgee_initBp(formula = Y ~ X1 + X2, 
                            Y=data.unsp$data$Y, Cov=data.unsp$data$Cov,
                            s=data.unsp$data$funcArg, subjID=data.unsp$data$subjID,
                            Tij=data.unsp$data$Tij, numLongiPoints = numLongiPoints, n.Bp=300)
end_time <- Sys.time()
time.initBp <- difftime(end_time, start_time, units = "secs")
    
Y.mean.initBp.ise <- mean.ise(initBp.estim$Y.mean.init,data.unsp$mean, numFunctPoints)
beta.initBp.acc <- beta.acc(beta.hat=initBp.estim$beta.init,beta.true=data.unsp$beta, numFunctPoints=numFunctPoints)
confidence.initBp <- eval.betaBp.ci(beta.lower=initBp.estim$beta.band.lower, beta.upper=initBp.estim$beta.band.upper,
                                        beta.true=data.unsp$beta, numFunctPoints=numFunctPoints)
initBp <- c(time.initBp, Y.mean.initBp.ise, beta.initBp.acc, confidence.initBp$length.ave, confidence.initBp$coverage.rate)
    
    

########################################
#iid estimate
########################################
start_time <- Sys.time()
iid.estim <- fgee(formula = Y ~ X1 + X2, Y=data.unsp$data$Y, Cov=data.unsp$data$Cov,
                  s=data.unsp$data$funcArg, subjID=data.unsp$data$subjID, 
                  Tij=data.unsp$data$Tij, numLongiPoints=numLongiPoints,
                  FPCA1.pve = 0.95, FPCA2.pve = 0.95, FPCA1.knots = 20, FPCA2.knots = 15,
                  corstr = "independent")
end_time <- Sys.time()
time.iid <- difftime(end_time, start_time, units = "secs")
    
result.iid <- eval.accuracy(phi.hat=list(phi_1=iid.estim$phi[,1],phi_2=iid.estim$phi[,2]), phi=data.unsp$phi,
                            numLongiPoints=41, numFunctPoints=101, 
                            Y.hat=iid.estim$Y.hat, Y.star=data.unsp$Y.star,
                            Y.mean=data.unsp$mean, Y.mean.hat=iid.estim$Y.mean,Y.mean.init=iid.estim$Y.mean.init,
                            beta.hat=iid.estim$beta, beta.init=iid.estim$beta.init, beta=data.unsp$beta)
    
confidence.iid <- eval.beta.ci(beta.hat=iid.estim$beta,beta.se=iid.estim$beta.se,beta.true=data.unsp$beta,numFunctPoints=numFunctPoints)
iid <- c(time.iid, result.iid$mean.insmpl.acc, result.iid$beta.acc, confidence.iid$length.ave, confidence.iid$coverage.rate)
    
    
########################################
#exchangable estimate
########################################

start_time <- Sys.time()
exch.estim <- fgee(formula = Y ~ X1 + X2, Y=data.unsp$data$Y, Cov=data.unsp$data$Cov,
                   s=data.unsp$data$funcArg, subjID=data.unsp$data$subjID, 
                   Tij=data.unsp$data$Tij, numLongiPoints=numLongiPoints,
                   FPCA1.pve = 0.95, FPCA2.pve = 0.95, FPCA1.knots = 20, FPCA2.knots = 15,
                   corstr = "exchangeable")
end_time <- Sys.time()
time.exch <- difftime(end_time, start_time, units = "secs")

result.exch <- eval.accuracy(phi.hat=list(phi_1=exch.estim$phi[,1],phi_2=exch.estim$phi[,2]), phi=data.unsp$phi,
                             numLongiPoints=41, numFunctPoints=101, 
                             Y.hat=exch.estim$Y.hat, Y.star=data.unsp$Y.star,
                             Y.mean=data.unsp$mean, Y.mean.hat=exch.estim$Y.mean,Y.mean.init=exch.estim$Y.mean.init,
                             beta.hat=exch.estim$beta, beta.init=exch.estim$beta.init, beta=data.unsp$beta)

confidence.exch <- eval.beta.ci(beta.hat=exch.estim$beta,beta.se=exch.estim$beta.se,beta.true=data.unsp$beta,numFunctPoints=numFunctPoints)
exch <- c(time.exch, result.exch$mean.insmpl.acc, result.exch$beta.acc, confidence.exch$length.ave, confidence.exch$coverage.rate)
    

########################################
#unspecified estimate
########################################
#with discrete
start_time <- Sys.time()
unsp.estim <- fgee(formula = Y ~ X1 + X2, Y=data.unsp$data$Y, Cov=data.unsp$data$Cov,
                   s=data.unsp$data$funcArg, subjID=data.unsp$data$subjID, 
                   Tij=data.unsp$data$Tij, numLongiPoints=numLongiPoints,
                   FPCA1.pve = 0.95, FPCA2.pve = 0.95, FPCA1.knots = 20, FPCA2.knots = 15,
                   corstr = "unspecified")
end_time <- Sys.time()
time.unsp <- difftime(end_time, start_time, units = "secs")

result.unsp <- eval.accuracy(phi.hat=list(phi_1=unsp.estim$phi[,1],phi_2=unsp.estim$phi[,2]), phi=data.unsp$phi,
                             numLongiPoints=41, numFunctPoints=101, 
                             Y.hat=unsp.estim$Y.hat, Y.star=data.unsp$Y.star,
                             Y.mean=data.unsp$mean, Y.mean.hat=unsp.estim$Y.mean,Y.mean.init=unsp.estim$Y.mean.init,
                             beta.hat=unsp.estim$beta, beta.init=unsp.estim$beta.init, beta=data.unsp$beta)

confidence.unsp <- eval.beta.ci(beta.hat=unsp.estim$beta,beta.se=unsp.estim$beta.se,beta.true=data.unsp$beta,numFunctPoints=numFunctPoints)

unsp <- c(time.unsp, result.unsp$mean.insmpl.acc, result.unsp$beta.acc, confidence.unsp$length.ave, confidence.unsp$coverage.rate)


