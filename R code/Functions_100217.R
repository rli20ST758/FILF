
###########################################################################
#function to generate data according to the approach described in 
#Notes 3 summary
###########################################################################
library(MASS)
############################################################################
GenerateData <- function(Nsubj=100, numFunctPoints = 101, min_visit=8, max_visit=12,
                         numLongiPoints = 41, sigma_sq = 1.5, 
                         sigma_z11 = 3, sigma_z12 = 1.5, 
                         sigma_z21 = 2, sigma_z22 = 1,
                         corstr = c("unspecified", "exchangeable", "independent", "laplace")){
  #Arguments:
  #Nsubj -- number of subjects, default to 100
  #numFunctPoints - number of points at which functions are observed
  #min_visit = 8 #(min # visits is 8) - minimum number of visits for each subject
  #max_visit = 12 #(max # visits is 12) - maximum number of visits for each subject
  #numLongiPoints -- number of points to evaluate inner eigenfunctions (for visits) on 
  #inner (visits) eigenfunction scores
  #sigma_sq -- white noise variance defaults to 1.5 for SNR = 5 (see Ana-Maria paper)
  #sigma_z11 <- 3
  #sigma_z12 <- 1.5
  #sigma_z21 <- 2
  #sigma_z22 <- 1
  #corstr = c("unspecified", "exchangeable", "independent", "laplace") -- type of covariance structure to generate data from
  ############################################################
  psi = NULL 
  xi = NULL
  zeta = NULL
  sigma = NULL
  #time points at which functional responses are collected
  s <- seq(0, 1, length = numFunctPoints)
  
  ########################################
  #Select time points for visits j
  ########################################
  #vector of visit argvalues T of length numLongiPoints
  T <- seq(0, 1, len = numLongiPoints) 
    
  #select number of visits per subject from uniform distribution
  runifdisc <- function(n, min=0, max=1) {sample(min:max, n, replace=TRUE)}
  m_i <- runifdisc(n=Nsubj, min = min_visit, max = max_visit)
  
  #vector of subject visits
  subjID <- rep(1:Nsubj, m_i)
  n <- length(subjID)
  
  #select for each subject i from all possible argvalues T, subject visits of length m_i
  Tij <- lapply(m_i, function(x){sort(sample(T, size = x, replace = FALSE))})
  T0ij <- Tij #for laplace corvariance matrix
  #Keep indicators for the visit number selected
  visitID <- lapply(Tij, function(x) which(T %in% x))
  
  #convert into a vector
  Tij <- unlist(Tij)
  visitID <- unlist(visitID)
  
  #create column id names 
  times <- paste("t.", 1:length(s), sep = "")
  
  #create row id names
  visits <- paste("S.", subjID, "v.", visitID , sep="")
  
  ########################################
  #Define k=2 outer basis functions 
  ########################################
  b_per <- 2
  phi_1 = function(s) {rep(1, length(s))}
  phi_2 = function(s, b_per=2) {sqrt(2) * sin (b_per*pi*s)}
  #phi values
  phi <- list()
  phi$phi_1 <- phi_1(s)
  phi$phi_2 <- phi_2(s)
  ########################################
  #Define covariates
  ########################################
  Cov1.pred <- rnorm(Nsubj)
  Cov1 <- rep(Cov1.pred, m_i)
  a <- 1
  pho <- 0.7
  ar.sim <- unlist(lapply(m_i,function(b){arima.sim(model=list(ar=pho),n=b)}))
  Cov2 <- a*Tij + ar.sim
  Cov <- cbind(Cov1,Cov2)
  rownames(Cov) <- visits
  colnames(Cov) <- c("X1", "X2")
  
  ########################################
  #Define mean response functions
  ########################################
  
  beta0 <- sqrt(2)*cos(3*pi*s) + 1
  beta1 <- 2 + cos(2*pi*s)
  beta2 <- 2 + sin(pi*s)
  
  ################################################
  #Define mean function at each time (s) and visit (T) point 
  ################################################
  
  mean <- Cov%*%t(cbind(beta1,beta2)) + matrix(rep(beta0,n),nr=n,byrow = TRUE) 
  rownames(mean) <- visits
  colnames(mean) <- times
  
  if(corstr == "unspecified"){
    ########################################
    #Define time varying components
    #xi_ik(T) = zeta_ik1 psi_k1(T) + zeta_ik2 psi_k2(T)
    ########################################
  
    #define psi_k1 and psi_k2
    psi_11 <- function(T) {sqrt(2) * cos (2*pi*T)}
    psi_12 <- function(T) {sqrt(2) * sin (2*pi*T)}
    psi_21 <- function(T) {sqrt(2) * cos (4*pi*T)}
    psi_22 <- function(T) {sqrt(2) * sin (4*pi*T)}
  
    #calculate eigenfunctions values at each point
    psi <- list()
    psi$psi_11 <- psi_11(T)
    psi$psi_12 <- psi_12(T) 
    psi$psi_21 <- psi_21(T)
    psi$psi_22 <- psi_22(T) 
  
    #generate zetas with n = #subjects from normal distrubution with mean 0 and variances defined at the beginning of the script
    #define zeta_k1 and zeta_k2
    zeta_11 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11-1))
    zeta_12 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z12))
    zeta_21 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21-0.5))
    zeta_22 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z22))
    #zeta values that will be used for evaluation, inner eigenfunctions scores
    zeta <- list()
    zeta$zeta_11 <- zeta_11
    zeta$zeta_12 <- zeta_12
    zeta$zeta_21 <- zeta_21
    zeta$zeta_22 <- zeta_22
    
    #create xi's
    xi_1 <-  rep(zeta_11, m_i)*psi$psi_11[visitID] + rep(zeta_12, m_i)*psi$psi_12[visitID]
    names(xi_1) <- visits
    xi_2 <-  rep(zeta_21, m_i)*psi$psi_21[visitID] + rep(zeta_22, m_i)*psi$psi_22[visitID]
    names(xi_2) <- visits
    #xi values that will be used for evaluation
    xi <-list()
    xi$xi_1 <- xi_1
    xi$xi_2 <- xi_2

    ########################################
    #define X_i(s, T_ij)
    #X_i(s, T_ij) = xi_1*phi_1 + xi_2*phi_2
    ########################################
    X <- xi_1%*%matrix(phi_1(s), nrow=1) + xi_2%*%matrix(phi_2(s), nrow=1)
    rownames(X) <- visits
    colnames(X) <- times
    
    ########################################
    #Generate random error terms
    #there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    #noise/signal = 1:5
    #inner.noise_1 <- rnorm(n = n, mean = 0, sd = sqrt((sigma_z11+sigma_z12)/6))
    #inner.noise_2 <- rnorm(n = n, mean = 0, sd = sqrt((sigma_z21+sigma_z22)/6))
    inner.noise_1 <- rnorm(n = n, mean = 0, sd = sqrt(1))
    inner.noise_2 <- rnorm(n = n, mean = 0, sd = sqrt(0.5))
    epsilon.inner <- (inner.noise_1)%*%matrix(phi_1(s), nrow=1) + (inner.noise_2)%*%matrix(phi_2(s), nrow=1)
    rownames(epsilon.inner) <- visits
    colnames(epsilon.inner) <- times
    
    #outer noise
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    #assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    #combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon.inner_ij(s) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon.inner + epsilon #response function used for evaluation (matrix form)
    Y.star = mean + X #used for accuracy evaluation
    }#end if corstr == unspecified
  
  else if (corstr == "exchangeable"){
    #k=1
    #subject specific random effects
    xi_i1_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11))
    names(xi_i1_coef) <- paste("S_", 1:Nsubj, sep="")
    xi_i1 <- rep(xi_i1_coef, m_i)
    names(xi_i1) <- visits
    #subject-visit specific effects
    xi_ij1 <- rnorm(sum(m_i), mean = 0, sd = sqrt(sigma_z12))
    names(xi_ij1) <- visits
    #k=2
    #define subject-visit specific effects
    xi_i2_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21))
    names(xi_i2_coef) <- paste("S_", 1:Nsubj, sep="")
    xi_i2 <- rep(xi_i2_coef, m_i)
    names(xi_i2) <- visits
    #subject-visit specific effects
    xi_ij2 <- rnorm(sum(m_i), mean = 0, sd = sqrt(sigma_z22))  
    names(xi_ij2) <- visits
    
    #xi values that will be used for evaluation
    #subject scores
    xi <-list()
    xi$xi_i1 <- xi_i1_coef
    xi$xi_i2 <- xi_i2_coef
    #subject-visit specific
    xi$xi_ij1 <- xi_ij1
    xi$xi_ij2 <- xi_ij2
    
    #define X_i(s, T_ij)
    #subject specific random process U_i
    U_i = t(matrix(phi_1(s))%*%xi_i1) + t(matrix(phi_2(s))%*%xi_i2)
    rownames(U_i) <- visits
    colnames(U_i) <- times
    
    #subject-visit specific random process V_ij
    V_ij <-  t(matrix(phi_1(s))%*%xi_ij1) + t(matrix(phi_2(s))%*%xi_ij2)
    rownames(V_ij) <- visits
    colnames(V_ij) <- times
    #response X = U_i + V_ij 
    X = U_i + V_ij 
    ########################################
    #Generate random error terms
    #there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    #assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    #combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon #response function used for evaluation (matrix form)
    Y.star = mean + X #used for accuracy evaluation
    }#end if corstr == exchangeable
  
  else if (corstr == "independent"){
    #k=1
    #subject-visit specific effects
    xi_ij1 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z11+sigma_z12))
    names(xi_ij1) <- visits
    #k=2
    #subject-visit specific effects
    xi_ij2 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z21+sigma_z22))  
    names(xi_ij2) <- visits
    
    #xi values that will be used for evaluation
    #subject-visit specific scores
    xi <-list()
    xi$xi_ij1 <- xi_ij1
    xi$xi_ij2 <- xi_ij2
    
    #define X_i(s, T_ij)
    #subject-visit specific random process V_ij
    V_ij <-  t(matrix(phi_1(s))%*%xi_ij1) + t(matrix(phi_2(s))%*%xi_ij2)
    rownames(V_ij) <- visits
    colnames(V_ij) <- times
    #response X = V_ij 
    X = V_ij 
    ########################################
    #Generate random error terms
    #there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    #assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    #combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon #response function used for evaluation (matrix form)
    Y.star = mean + X #used for accuracy evaluation
    }#end if corstr == independent
  
  else if (corstr == "laplace"){
    ########################################
    #Define time varying components
    #cov(xi_ik(T),xi_ik(T')) = sigma_zk1*exp(-|T-T'|)
    ########################################
    #create xi's
    xi_1 <- unlist(lapply(1:Nsubj,function(x){
      cov <- exp(-abs(outer(T0ij[[x]],T0ij[[x]],"-")))
      mvrnorm(1, rep(0,length(T0ij[[x]])), Sigma=(sigma_z11+sigma_z12)*cov)}))
    names(xi_1) <- visits
    xi_2 <- unlist(lapply(1:Nsubj,function(x){
      cov <- exp(-abs(outer(T0ij[[x]],T0ij[[x]],"-")))
      mvrnorm(1, rep(0,length(T0ij[[x]])), Sigma=(sigma_z21+sigma_z22)*cov)}))
    names(xi_2) <- visits
    
    #xi values that will be used for evaluation
    xi <-list()
    xi$xi_1 <- xi_1
    xi$xi_2 <- xi_2
    
    ########################################
    #define X_i(s, T_ij)
    #X_i(s, T_ij) = xi_1*phi_1 + xi_2*phi_2
    ########################################
    X <- xi_1%*%matrix(phi_1(s), nrow=1) + xi_2%*%matrix(phi_2(s), nrow=1)
    rownames(X) <- visits
    colnames(X) <- times
    
    ########################################
    #Generate random error terms
    #there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
    ########################################
    epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                      nrow = n, ncol = numFunctPoints)
    #assign row and column names to error terms matrix epsilon_ij
    rownames(epsilon) <- visits
    colnames(epsilon) <- times
    
    ########################################
    #combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
    ########################################
    Y = mean + X + epsilon #response function used for evaluation (matrix form)
    Y.star = mean + X #used for accuracy evaluation
    }#end if corstr == laplace
  

  return(list(data=list(subjID = subjID, Tij = Tij, visitID= visitID, funcArg = s,
               Y = Y, Cov = Cov), X = X, Y.star = Y.star, mean = mean, beta = cbind(beta0,beta1,beta2), 
              phi = phi,  xi = xi, psi = psi, zeta = zeta))
  
}#end of function GenerateData

#################
#function that calculates ISE for eigenfuctions and corrects for the sign to compare scores
#################

IMSE_fn <- function(fn_hat, fn, len){
  flip=FALSE
  dif1 <- sum((fn - fn_hat)^2)/len#if signs are the same
  dif2 <- sum((fn + fn_hat)^2)/len#if signs are opposite
  Ind <- which(c(dif1,dif2) == min(dif1, dif2))
  if(Ind ==2){flip = TRUE}
  return(list(IMSE = min(dif1, dif2), flip = flip))
}


#################################################
#function for evaluating accuracy of estimation
###############################################
eval.accuracy <- function(phi.hat, phi,numLongiPoints, numFunctPoints,
                          Y.hat, Y.star,
                          Y.mean, Y.mean.hat,Y.mean.init, 
                          beta.hat, beta.init, beta){
 
  #create holders
  phi.acc <- rep(0,2)
  beta.acc <- rep(0,3)
  beta.init.acc <- rep(0,3)
  #xi.acc <- rep(0,2)
  #psi.acc <- matrix(0, nrow =2, ncol =2)
  #rownames(psi.acc) <- c("k1", "k2") #two eigenfunctions 
  #colnames(psi.acc) <- c("f1", "f2") #same k
  #zeta.acc <- matrix(0, nrow =2, ncol =2)
  #rownames(zeta.acc) <- c("k1", "k2") #two eigenfunctions scores for the same k
  #colnames(zeta.acc) <- c("f1", "f2")
  #sigma.acc <- matrix(0, nrow =2, ncol =2)
  #rownames(sigma.acc) <- c("k1", "k2") #two eigenfunctions 
  #colnames(sigma.acc) <- c("Subj", "SubjVis") #same k
  #x <- NULL
  #x.pred <- NULL
  y_in_smp <- NULL
  #y_out_smp <- NULL
  mu_in_smp <- NULL
  mu_init_smp <- NULL
  #mu_out_smp <- NULL
  #accuracy of estimating eigenfunctions
  if(!is.null(phi.hat) & !is.null(phi)){
    res <- IMSE_fn(fn_hat = phi.hat$phi_1, fn = phi$phi_1, len = numFunctPoints)    
    #flip sign of eigenfunctions for plots and scores calculation if res$flip == TRUE
    if(res$flip == TRUE){phi.hat$phi_1 <- -phi.hat$phi_1}
    #save results
    phi.acc[1] <- res$IMSE
    #second eigenfunction, k=2
    res <- IMSE_fn(fn_hat = phi.hat$phi_2, fn=phi$phi_2, len = numFunctPoints)    
    #flip sign of eigenfunctions for plots and scores calculation if  res$flip == TRUE
    if(res$flip == TRUE){phi.hat$phi_2 <- -phi.hat$phi_2}
    #save results
    phi.acc[2] <- res$IMSE
  }#end if !null phi
  
  #evaluating accuracy for estimating scores eigenfunctions psi for k=1
  #if(!is.null(psi.hat) & !is.null(psi)){
  #  res <- IMSE_fn(fn_hat = psi.hat$psi_11, fn=psi$psi_11, len = numLongiEvalPoints)    
  #  #flip sign of eigenfunctions for plots and scores calculation if  res$flip == TRUE
  #  if(res$flip == TRUE){psi.hat$psi_11 <- -psi.hat$psi_11
  #                       if(!is.null(zeta.hat)){zeta.hat$zeta_11 <- -zeta.hat$zeta_11}}
  #  #save results
  #  psi.acc[1,1] <- res$IMSE
  #  #second psi for k=1
  #  if(!is.null(psi.hat$psi_12)){
  #    res <- IMSE_fn(fn_hat = psi.hat$psi_12, fn=psi$psi_12, len = numLongiEvalPoints)    
  #    #flip sign of eigenfunctions for plots and scores calculation if  res$flip == TRUE
  #    if(res$flip == TRUE){psi.hat$psi_12 <- -psi.hat$psi_12
  #                         if(!is.null(zeta.hat)){zeta.hat$zeta_12 <- -zeta.hat$zeta_12}}
  #    #save results
  #    psi.acc[1,2] <- res$IMSE
  #  }#end if !is.null(psi.hat$psi_1_2)
  #  #evaluating accuracy for estimating visits covariance eigenfunctions psi for k=2
  #  res <- IMSE_fn(fn_hat = psi.hat$psi_21, fn=psi$psi_21, len = numLongiEvalPoints)    
  #  #flip sign of eigenfunctions for plots and scores calculation if  res$flip == TRUE
  #  if(res$flip == TRUE){psi.hat$psi_21 <- -psi.hat$psi_21
  #                       if(!is.null(zeta.hat)){zeta.hat$zeta_21 <- -zeta.hat$zeta_21}}
  #  #save results
  # psi.acc[2,1] <- res$IMSE
    
  #  if(!is.null(psi.hat$psi_22)){
  #     res <- IMSE_fn(fn_hat = psi.hat$psi_22, fn=psi$psi_22, len = numLongiEvalPoints)    
  #  #flip sign of eigenfunctions for plots and scores calculation if  res$flip == TRUE
  #     if(res$flip == TRUE){psi.hat$psi_22 <- -psi.hat$psi_22
  #                         if(!is.null(zeta.hat)){zeta.hat$zeta_22 <- -zeta.hat$zeta_22}}
  #    #save results
  #    psi.acc[2,2] <- res$IMSE
  #   }
  # }#end evaluation for accuracy of psi
  
  ###############################
  #add accuracy for beta functions
  ###############################
  
  if(!is.null(beta.hat) & !is.null(beta) ){
    #beta0(s)
    beta.acc[1] <- sum((beta.hat[,1]-beta[,1])^2)/numFunctPoints
    #beta1(s)
    beta.acc[2] <- sum((beta.hat[,2]-beta[,2])^2)/numFunctPoints
    #beta2(s)
    beta.acc[3] <- sum((beta.hat[,3]-beta[,3])^2)/numFunctPoints
  }
  #out of sample prediction
  if(!is.null(beta.init) & !is.null(beta) ){
    #beta0(s)
    beta.init.acc[1] <- sum((beta.init[,1]-beta[,1])^2)/numFunctPoints  
    #beta1(s)
    beta.init.acc[2] <- sum((beta.init[,2]-beta[,2])^2)/numFunctPoints  
    #beta2(s)
    beta.init.acc[3] <- sum((beta.init[,3]-beta[,3])^2)/numFunctPoints  
  }
  
  ##################################################
  #evaluating accuracy for estimating scores xi
  ##################################################  
  # 
  # if(!is.null(xi.hat) & !is.null(xi) &!is.null(subjID)){
  #   if (model =="unspecified"){
  #     #accuracy in estimating xi_hat
  #     #k=1
  #     dif <- (xi.hat$xi_1 - xi$xi_1)^2
  #     #average over visits
  #     df <- data.frame(dif, subjID)
  #     names(df) <- c("dif", "subj")
  #     int_dif <- aggregate(dif ~ subj, df, mean)$dif
  #     #save results
  #     xi.acc[1] <- sum(int_dif)/n #k=1
  #     #k=2
  #     dif <- (xi.hat$xi_2 - xi$xi_2)^2
  #     #average over visits
  #     df <- data.frame(dif, subjID)
  #     names(df) <- c("dif", "subj")
  #     int_dif <- aggregate(dif ~ subj, df, mean)$dif
  #     #save results
  #     xi.acc[2] <- sum(int_dif)/n #k=2
  #   }#end if model = unspecified
  #   
  #   if (model =="exchangable"){
  #     xi.acc <- matrix(0, nrow =2, ncol =2)
  #     rownames(xi.acc) <- c("k1", "k2") #two eigenfunctions 
  #     colnames(xi.acc) <- c("Subj", "SubjVis") #same k
  #     #accuracy in estimating xi_hat
  #     #k=1
  #     #subject-specific scores
  #     dif <- (xi.hat$xi_i1 - xi$xi_i1)^2
  #     #save results
  #     xi.acc[1,1] <- sum(dif)/n #k=1
  #     
  #     #subject-visit specific scores
  #     dif <- (xi.hat$xi_ij1 - xi$xi_ij1)^2
  #     #average over visits
  #     df <- data.frame(dif, subjID)
  #     names(df) <- c("dif", "subj")
  #     int_dif <- aggregate(dif ~ subj, df, mean)$dif
  #     #save results
  #     xi.acc[1,2] <- sum(int_dif)/n #k=1
  #     
  #     #k=2
  #     #subject-specific scores
  #     dif <- (xi.hat$xi_i2 - xi$xi_i2)^2
  #     #save results
  #     xi.acc[2,1] <- sum(dif)/n #k=2
  #     
  #     #subject-visit specific scores
  #     dif <- (xi.hat$xi_ij2 - xi$xi_ij2)^2
  #     #average over visits
  #     df <- data.frame(dif, subjID)
  #     names(df) <- c("dif", "subj")
  #     int_dif <- aggregate(dif ~ subj, df, mean)$dif
  #     #save results
  #     xi.acc[2,2] <- sum(int_dif)/n #k=2
  #   }#end if model = exchangable
  # }#end accuracy of estimation of xi
  
  ###################################################
  #evaluating accuracy of sigma for exchangable model
  ###################################################
  # if(!is.null(sigma.hat) & !is.null(sigma)){
  #   
  #   sigma.acc[1,1] <- abs(sigma.hat$k1[1] - sigma$k1[1])
  #   sigma.acc[1,2] <- abs(sigma.hat$k1[2] - sigma$k1[2])
  #   sigma.acc[2,1] <- abs(sigma.hat$k2[1] - sigma$k2[1])
  #   sigma.acc[2,2] <- abs(sigma.hat$k2[2] - sigma$k2[2])
  # }
    
  ########################################################################
  #evaluating accuracy for estimating scores of eigenfunction scores zeta
  ######################################################################## 
  # if(!is.null(zeta.hat) & !is.null(zeta)){
  #   #first eigenfunction k=1  
  #   dif <- (zeta.hat$zeta_11 - zeta$zeta_11)^2
  #   #save results
  #   zeta.acc[1,1] <- sum(dif)/n
  #   
  #   if(!is.null(zeta.hat$zeta_12)){
  #     dif <-  (zeta.hat$zeta_12 - zeta$zeta_12)^2
  #     #save results
  #     zeta.acc[1,2] <- sum(dif)/n}
  #   
  #   #second eigenfunction k=2 
  #   dif <- (zeta.hat$zeta_21 - zeta$zeta_21)^2
  #   #save results
  #   zeta.acc[2,1] <- sum(dif)/n
  #   
  #   if(!is.null(zeta.hat$zeta_22)){
  #     dif <-  (zeta.hat$zeta_22 - zeta$zeta_22)^2
  #     #save results
  #     zeta.acc[2,2] <- sum(dif)/n}
  # }#end accuracy of zetas
  
  #####################################
  #accuracy of mean estimation
  #####################################
  #in sample initial mean estimate
  if(!is.null(Y.mean) & !is.null(Y.mean.init) ){
    dif <- (Y.mean - Y.mean.init)^2
    #intergate over time s
    int_dif <- apply(dif, 1, sum)/numFunctPoints
    #save results
    mu_init_smp <- sum(int_dif)/dim(Y.mean)[1] #averaging over visits
  }
  
  #in sample final mean estimate
  if(!is.null(Y.mean) & !is.null(Y.mean.hat) ){
    dif <- (Y.mean - Y.mean.hat)^2
    #intergate over time s
    int_dif <- apply(dif, 1, sum)/numFunctPoints
    #save results
    mu_in_smp <- sum(int_dif)/dim(Y.mean)[1] #averaging over visits
  }
  
  # #out of sample mean
  # if(!is.null(Y.pred.mean) & !is.null(Y.pred.mean.hat) ){
  #   dif <- (Y.pred.mean - Y.pred.mean.hat)^2
  #   #intergate over time s
  #   int_dif <- apply(dif, 1, sum)/NumFunctPoints
  #   #save results
  #   mu_out_smp <- sum(int_dif)/dim(Y.pred.mean)[1] #averaging over visits
  # }
  
  #######################################################
  #accuracy of subject-visit specific curves estimation
  #######################################################
  # #in sample subject-visit specific mean
  # if(!is.null(X) & !is.null(X.hat) ){
  #   dif <- (X - X.hat)^2
  #   #intergate over time s
  #   int_dif <- apply(dif, 1, sum)/NumFunctPoints
  #   #save results
  #   x <- sum(int_dif)/dim(X)[1] #averaging over visits
  # }
  # 
  # #out of  sample subject-visit specific mean
  # if(!is.null(X.pred) & !is.null(X.pred.hat) ){
  #   dif <- (X.pred - X.pred.hat)^2
  #   #intergate over time s
  #   int_dif <- apply(dif, 1, sum)/NumFunctPoints
  #   #save results
  #   x.pred <- sum(int_dif)/dim(X.pred)[1] #averaging over number of curves
  # }
  
  ########################################################
  #add accuracy for in sample and out of sample prediction
  ########################################################
  
  if(!is.null(Y.hat) & !is.null(Y.star) ){
    dif <- (Y.hat - Y.star)^2
    #intergate over time s
    int_dif <- apply(dif, 1, sum)/numFunctPoints
    #save results
    y_in_smp <- sum(int_dif)/dim(Y.star)[1] #averaging over visits
  }
  #out of sample prediction
  # if(!is.null(Y.pred.hat) & !is.null(Y.pred.star) ){
  #   dif <- (Y.pred.hat - Y.pred.star)^2
  #   #intergate over time s
  #   int_dif <- apply(dif, 1, sum)/NumFunctPoints
  #   #save results
  #   y_out_smp <- sum(int_dif)/dim(Y.pred.hat)[1] #averaging over points out of sample
  # }
  #return results
  return(list(phi.acc = phi.acc,
              mean.insmpl.acc = mu_in_smp, mean.insmpl.init.acc = mu_init_smp, 
              insmpl.acc = y_in_smp, beta.acc = beta.acc, beta.init.acc = beta.init.acc))
}


#####################################################################################
#function for evaluating accuracy of beta estimation (confidence interval & acc)
#####################################################################################
eval.beta.ci <- function(beta.hat, beta.se, beta.true, numFunctPoints){
  #fixed effects functions  beta(s)
  beta0 <- beta.true[,1]
  beta1 <- beta.true[,2]
  beta2 <- beta.true[,3]
  
  
  #length of 95% confidence interval 
  length.ave <- colSums(2*qnorm(0.975)*beta.se)/numFunctPoints
    
  #beta0(s)
  #confidence band 
  lower <- beta.hat[,1]-qnorm(0.975)*beta.se[,1]
  upper <- beta.hat[,1]+qnorm(0.975)*beta.se[,1]
  #count number of points 
  k1 <- 0 
  for(i in 1:numFunctPoints){
    if(beta0[i]>=lower[i] & beta0[i]<=upper[i]) 
      k1 <- k1+1
  }
  coverage.rate.beta0 <- k1/numFunctPoints
  
  #beta1(s)
  #confidence band 
  lower <- beta.hat[,2]-qnorm(0.975)*beta.se[,2]
  upper <- beta.hat[,2]+qnorm(0.975)*beta.se[,2]
  #count number of points 
  k2 <- 0 
  for(i in 1:numFunctPoints){
    if(beta1[i]>=lower[i] & beta1[i]<=upper[i]) 
      k2 <- k2+1
  }
  coverage.rate.beta1 <- k2/numFunctPoints
  
  #beta2(s)
  #confidence band 
  lower <- beta.hat[,3]-qnorm(0.975)*beta.se[,3]
  upper <- beta.hat[,3]+qnorm(0.975)*beta.se[,3]
  #count number of points 
  k3 <- 0 
  for(i in 1:numFunctPoints){
    if(beta2[i]>=lower[i] & beta2[i]<=upper[i]) 
      k3 <- k3+1
  }
  coverage.rate.beta2 <- k3/numFunctPoints
  
  coverage.rate <- c(coverage.rate.beta0,coverage.rate.beta1,coverage.rate.beta2)
  
  return(list(length.ave=length.ave,coverage.rate=coverage.rate))
}


#beta.acc
beta.acc <- function(beta.hat, beta.true, numFunctPoints){
  if(!is.null(beta.hat) & !is.null(beta.true)){
    #beta(s)
    beta.acc <- c(sum((beta.hat[,1]-beta.true[,1])^2)/numFunctPoints,
                  sum((beta.hat[,2]-beta.true[,2])^2)/numFunctPoints,
                  sum((beta.hat[,3]-beta.true[,3])^2)/numFunctPoints)
    return(beta.acc)
  }
}

#####################################################
#function for evaluating accuracy of mean estimation
#####################################################
mean.ise <- function(Y.mean.hat, Y.mean.true, numFunctPoints){
  dif <- (Y.mean.hat - Y.mean.true)^2
  #intergate over time s
  int_dif <- apply(dif, 1, sum)/numFunctPoints
  #save results
  mu_smp <- sum(int_dif)/dim(Y.mean.hat)[1] #averaging over visits
  return(mu_smp)
}


#######################################################################
#function for evaluating accuracy of bootstrap beta confidence interval
#######################################################################
eval.betaBp.ci <- function(beta.lower, beta.upper, beta.true, numFunctPoints){
  #fixed effects functions  beta(s)
  beta0 <- beta.true[,1]
  beta1 <- beta.true[,2]
  beta2 <- beta.true[,3]
  
  beta0.band <- cbind(beta.lower[,1], beta.upper[,1])
  beta1.band <- cbind(beta.lower[,2], beta.upper[,2])
  beta2.band <- cbind(beta.lower[,3], beta.upper[,3])
  
  #length of 95% confidence interval 
  length.ave <- c(sum(beta0.band[,2]-beta0.band[,1]),sum(beta1.band[,2]-beta1.band[,1]),sum(beta2.band[,2]-beta2.band[,1]))/numFunctPoints
  
  #beta0(s)
  #count number of points 
  k1 <- 0 
  for(i in 1:numFunctPoints){
    if(beta0[i]>=beta0.band[i,1] & beta0[i]<=beta0.band[i,2]) 
      k1 <- k1+1
  }
  coverage.rate.beta0 <- k1/numFunctPoints
  
  #beta1(s)
  k2 <- 0 
  for(i in 1:numFunctPoints){
    if(beta1[i]>=beta1.band[i,1] & beta1[i]<=beta1.band[i,2]) 
      k2 <- k2+1
  }
  coverage.rate.beta1 <- k2/numFunctPoints
  
  #beta2(s)
  k3 <- 0 
  for(i in 1:numFunctPoints){
    if(beta2[i]>=beta2.band[i,1] & beta2[i]<=beta2.band[i,2]) 
      k3 <- k3+1
  }
  coverage.rate.beta2 <- k3/numFunctPoints
  
  coverage.rate <- c(coverage.rate.beta0,coverage.rate.beta1,coverage.rate.beta2)
  
  return(list(length.ave=length.ave,coverage.rate=coverage.rate))
}



