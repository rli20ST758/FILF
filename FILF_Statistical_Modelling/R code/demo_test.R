
library(refund)
library(face)
library(fields)
library(mgcv)
library(nlme) 
library(MASS)
library(matrixcalc)
library(Matrix)
library(Bolstad)
library(splines)



#setwd("~/Dowdonw/R code/test")
source("GenData.R")
source("calc.mean.R")
source("select.knots.R")
source("trunc.mat.R")
source("ts.fit.R")
source("calc.RA.R")
source("calc.sigsq.R")
source("irreg2mat.mod.R")
source("p.bs.R")
source("test.cov.exch.R")
source("exch.R")


set.seed(1)
    

numFunctPoints = 101
Nsubj = 100
numLongiPoints = 41
FPCA1.pve = 0.95 
FPCA1.knots = 20
## additional setups
k0 <- 10 ##number of basis functions for fixed effects functions
bs <- "ps"

data.exch <- GenData(Nsubj=Nsubj, numFunctPoints = numFunctPoints, min_visit=3, max_visit=6, 
                     numLongiPoints = 41, sigma_sq = 1.5, delta = 0,
                     sigma_z11 = 3, sigma_z12 = 1.5, 
                     sigma_z21 = 2, sigma_z22 = 1, model = "exchangable")


data <- data.exch$data
Y <- data$Y
s <- data$funcArg
Cov <- data$Cov
Tij <- data$Tij
X1 <- Cov[,1]
X2 <- Cov[,2]
formula <- Y ~ X1 + X2
model.formula <- as.character(formula)

#initial mean estimates
fit_init <- pffr(formula, yind=s, 
                 bs.yindex = list(bs = bs, k = k0, m = c(2, 2))) ## initial fit
Y.mean.init <- fitted(fit_init) ## fitted means
Y.resid <- as.matrix(Y - Y.mean.init)

#score
fpca_margin <- fpca.face(Y.resid,center = FALSE, argvals=s, knots=FPCA1.knots, pve=FPCA1.pve)

score <- fpca_margin$scores/sqrt(length(s))
score1 <- data.frame(.value=score[,1],.index=Tij,.id=data$subjID)
score2 <- data.frame(.value=score[,2],.index=Tij,.id=data$subjID)

    

###############test structure of G_k(T,T')########################
test1 <- test.cov.exch(score1,numLongiPoints=numLongiPoints,nbs=1000, nb=10)
test2 <- test.cov.exch(score2,numLongiPoints=numLongiPoints,nbs=1000, nb=10)
    

###############test structure of G_k(T,T')########################
test3 <- test.cov.iid(score1,numLongiPoints=numLongiPoints,nbs=1000, nb=10)
test4 <- test.cov.iid(score2,numLongiPoints=numLongiPoints,nbs=1000, nb=10)






