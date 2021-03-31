library(MASS)
############################################################################
GenData <- function(Nsubj=100, numFunctPoints = 101, min_visit=8, max_visit=12,
                    numLongiPoints = 41, sigma_sq = 1.5, delta = 0,
                         sigma_z11 = 3, sigma_z12 = 1.5, 
                         sigma_z21 = 2, sigma_z22 = 1,
                         model = c("exchangable", "iid")){
  #Arguments:
  #Nsubj -- number of subjects, default to 100
  #numFunctPoints - number of points at which functions are observed
  #min_visit = 8 #(min # visits is 8) - minimum number of visits for each subject
  #max_visit = 12 #(max # visits is 12) - maximum number of visits for each subject
  #numLongiPoints - total number of possible visits. Subjects have uniform visits (from  min_visit to max_visit) 
  #anywhere in the range between [0,1] of length J 
  #sigma_sq -- white noise variance defaults to 1.5 for SNR = 5 (see Ana-Maria paper)
  #numLongiEvalPoints -- number of points to evaluate inner eigenfunctions (for visits) on 
  #inner (visits) eigenfunction scores
  #sigma_z11 <- 3
  #sigma_z12 <- 1.5
  #sigma_z21 <- 2
  #sigma_z22 <- 1
  #delta <- 0 -- the scalar controling the magnitude of deviation from the null model
  #model = c("exchangable", "iid") -- type of covariance structure to generate data from
  ############################################################
  #time points at which functional responses are collected
  s <- seq(0, 1, length = numFunctPoints)
  
  ########################################
  #Select time points for visits j
  ########################################
  #vector of visit argvalues T of length J 
  T <- seq(0, 1, len = numLongiPoints) 
  
  #select number of visits per subject from uniform distribution
  runifdisc <- function(n, min=0, max=1) {sample(min:max, n, replace=TRUE)}
  m_i <- runifdisc(n=Nsubj, min = min_visit, max = max_visit)
  
  #vector of subject visits
  subjID <- rep(1:Nsubj, m_i)
  n <- length(subjID)
  
  #select for each subject i from all possible argvalues T, subject visits of length m_i
  Tij <- lapply(m_i, function(x){sort(sample(T, size = x, replace = FALSE))})
  #Keep indicators for the visit number selected
  visitID <- lapply(Tij, function(x) which(T %in% x))
  #convert into a vector
  Tij <- unlist(Tij)
  visitID <- unlist(visitID)
  
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
  colnames(Cov) <- c("x1", "x2")
  
  ########################################
  #Define mean response functions
  ########################################
  beta0 <- sqrt(2)*cos(3*pi*s) + 1
  beta1 <- 2 + cos(2*pi*s)
  beta2 <- 2 + sin(pi*s)
  mean <- Cov%*%t(cbind(beta1,beta2)) + matrix(rep(beta0,n),nr=n,byrow = TRUE) 
  
  ########################################
  #Define deviation from the null model
  #z_i(T) = zeta_i1 psi_1(T) + zeta_i2 psi_2(T)
  ########################################
  #define phi_1 and phi_2
  #psi_1 <- sin(2*pi*T)
  #psi_2 <- sin(4*pi*T)
  psi_11 <- sqrt(2)*cos(2*pi*T)
  psi_12 <- sqrt(2)*sin(2*pi*T)
  psi_21 <- sqrt(2)*cos(4*pi*T)
  psi_22 <- sqrt(2)*sin(4*pi*T)
  #generate xis with n = #subjects from normal distrubution with mean 0 and variances defined at the beginning of the script
  #define zeta_1 and zeta_2
  zeta_11 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11))
  zeta_12 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z12))
  zeta_21 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21))
  zeta_22 <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z22))
  
  #create z_i(T_ij)
  z1 <- rep(zeta_11, m_i)*psi_11[visitID] + rep(zeta_12, m_i)*psi_12[visitID]
  z2 <- rep(zeta_21, m_i)*psi_21[visitID] + rep(zeta_22, m_i)*psi_22[visitID]
  ################################################
  #Define mean function at each time (s) and visit (T) point 
  ################################################
  
  if (model == "exchangable"){
    #k=1
    #subject specific random effects
    xi_i1_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z11))
    xi_i1 <- rep(xi_i1_coef, m_i) 
    #subject-visit specific effects
    xi_ij1 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z12))
    #k=2
    #define subject-visit specific effects
    xi_i2_coef <- rnorm(n = Nsubj, mean = 0, sd = sqrt(sigma_z21))
    xi_i2 <- rep(xi_i2_coef, m_i)
    #subject-visit specific effects
    xi_ij2 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z22))  
    
    #define X_i(s, T_ij)
    X = t(matrix(phi_1(s))%*%(xi_i1 + xi_ij1 + delta*z1)) + 
      t(matrix(phi_2(s))%*%(xi_i2 + xi_ij2 + delta*z2))
  }#end if model == exchangable
  
  if (model == "iid"){
    #k=1
    #subject-visit specific effects
    xi_ij1 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z11+sigma_z12))
    #k=2
    #subject-visit specific effects
    xi_ij2 <- rnorm(sum(m_i), mean =0, sd = sqrt(sigma_z21+sigma_z22))  
    
    #define X_i(s, T_ij)
    X <-  t(matrix(phi_1(s))%*%(xi_ij1+delta*z1)) + t(matrix(phi_2(s))%*%(xi_ij2+delta*z2))
  
  }#end if model == iid
  
  ########################################
  #Generate random error terms
  #there are (i=length(s))*(j=sum(m_i)) = #time points * #visits for all subjects
  ########################################
  epsilon <- matrix(rnorm(n*numFunctPoints, mean = 0, sd = sqrt(sigma_sq)), 
                    nrow = n, ncol = numFunctPoints)
  
  ########################################
  #combine to get data Y_ij(s) = mu(s,T_ij) + X_i(s, T_ij) + epsilon_ij(s)
  ########################################
  Y = mean + X + epsilon #response function used for evaluation (matrix form)
  Y.star = mean + X #used for accuracy evaluation
  
  return(list(data=list(subjID = subjID, Tij = Tij, visitID= visitID, funcArg = s,Y = Y, Cov = Cov),
              X = X, Y.star = Y.star, mean = mean, beta = cbind(beta0,beta1,beta2), phi = phi))
  
}#end of function GenerateData
