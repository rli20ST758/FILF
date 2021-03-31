library(face)
library(mgcv)

fgee <- function(formula, Y, Cov, s, subjID, Tij, numLongiPoints = 41, 
                 control = list(FPCA1.pve=0.95, FPCA2.pve=0.95,  
                                FPCA1.knots=20, FPCA2.knots=12),
                 corstr = c("unspecified", "exchangeable","independent")){
  
  # Arguments:
  # 1.fomular 
  # parameter fomular: a formula expression which indicates a linear relationship between Y and Cov
  #
  # 2.data 
  # parameter Y -- a matrix of which each row corresponds to one curve observed on a regular and dense grid 
  # (dimension of N by m; N = total number of observed functions; m = number of grid points) 
  # parameter Cov -- including all covariates of interest.      
  # parameter s -- argvals a string indicating the functional domain variable name in data   
  # parameter subjID -- subject id; vector of length N with each element corresponding a row of Y
  # parameter Tij -- actual time of visits at which a function is observed
  # parameter numLongiPoints -- total number of evaluation time points for visits; defaults to 41
  # 
  # 3.parameter control: a list with the named components: 
  # FPCA1.pve -- proportion of variance explained used to choose the number of principal components; 
  # defaults to 0.95; see fpca.face
  # FPCA2.pve -- parameter needed for corstr of type "unspecified", 
  # proportion of variance explained used to choose the number of principal components; 
  # defaults to 0.95; see fpca.sc
  # FPCA1.knots -- number of knots to use or the vectors of knots used for obtain 
  # a smooth estimate of a covarance function; defaults to 20; see fpca.face
  # FPCA2.knots -- parameter needed for corstr of type "unspecified", 
  # number of B-spline basis functions for estimation of the mean function and 
  # bivariate smoothing of the covariance surface; defaults to 12; see fpca.sc
  #
  # 4.model:
  # parameter corstr -- type of covariance structure considering in model
  ############################################################   
  
  # formula for initial estimate
  model.formula <- as.character(formula)
  stopifnot(model.formula[1] == "~" & length(model.formula) == 3)
  # compile the formula again
  formula <- paste("Y", model.formula[1], model.formula[3])
  # select covariates are used in formula
  Cov.fm <- c() 
  for(i in 1:ncol(Cov)){
    if(grepl(colnames(Cov)[i], model.formula[3])){
      Cov.fm <- c(Cov.fm, i)
    }
  }
  
  TT <- seq(min(Tij), max(Tij), length.out=numLongiPoints) #Tij close to TT[v]
  visitID <- sapply(Tij, function(a) which(abs(TT-a) == min(abs(TT-a))))
  psi.MF <- NULL #psi at Tij
  # holder for estimation results
  psi <- NULL
  # a simple vector to matrix conversion function
  vec2mat <- function(x){matrix(x,ncol=length(s),byrow=TRUE)}
  
  # additional setups
  k0 <- 10 #number of basis functions for fixed effects functions
  bs <- "ps"
  # vectorize data
  y <- c(t(Y))
  t <- rep(s,times=nrow(Y))
  for(i in 1:length(Cov.fm)){
    assign(colnames(Cov)[Cov.fm[i]], Cov[,Cov.fm[i]])
    assign(paste("x", i, sep=""), rep(Cov[,Cov.fm[i]], each=length(s))) 
  }
  ###########################################
  # Step 1: Initial mean estimate
  ##########################################

  # initial mean estimates
  fit_init <- pffr(as.formula(formula), yind=s, algorithm="bam",
                   bs.yindex = list(bs = bs, k = k0, m = c(2, 2))) #initial fit
  Y.mean_init <- predict(fit_init, type="terms")
  Y.mean.init <- fitted(fit_init) #fitted means
  Y.resid <- as.matrix(Y - Y.mean.init)

  # initial estimation of beta(s) and their stardard error
  plot.data <- {
    pdf(NULL)
    plot.init <- plot(fit_init,n=length(s), seWithMean=TRUE, pages=1)
    invisible(dev.off())
  }
  beta.init.0 <- plot.init[[1]]$fit + fit_init$coefficients[1]
  beta.init.cov <- matrix(unlist(lapply(1:length(Cov.fm), function(x) plot.init[[x+1]]$fit)),ncol=length(Cov.fm),byrow=FALSE)
  beta.init = cbind(beta.init.0, beta.init.cov)
  # to get stardard error of initial beta(s)
  beta.init.se <- matrix(unlist(lapply(1:(length(Cov.fm)+1), function(x) plot.init[[x]]$se/2)),ncol=(length(Cov.fm)+1),byrow=FALSE)
  ########################################
  # Step 2: Estimate marginal eigenfunctions
  ########################################
  fpca_margin <- fpca.face(Y.resid, center = FALSE, argvals=s, knots=control$FPCA1.knots, pve=control$FPCA1.pve)
  phi <- fpca_margin$efunctions*sqrt(length(s))
  lambda.phi <- fpca_margin$evalues/length(s)
  
  ########################################
  #Step 3: Re-estimate setup
  ########################################
  
  if(corstr == "unspecified"){
    ########################################
    # estimate conditional eigenfunctions using faces
    ########################################
    scores <- fpca_margin$scores/sqrt(length(s))
    # convert scores into data frame
    data_scores <- list()
    for(k in 1:ncol(scores))
      data_scores[[k]] <- data.frame(y = scores[,k], argvals = Tij, subj = subjID)
    
    # convert data into faces structure
    T_dense <- seq(0,1, length.out=numLongiPoints)
    # fit faces
    fpca_condition <- lapply(data_scores,function(x){
      face.sparse(x, argvals.new = T_dense, center = FALSE,
                  newdata = x, calculate.scores = TRUE, pve = control$FPCA2.pve,
                  knots = control$FPCA2.knots)
    })
    
    # add sigma as eigenvalue to select new eigenvalues and eigenfunctions
    psi <- lapply(fpca_condition,function(x){as.matrix(x$eigenfunctions)})
    lambda.psi <- lapply(fpca_condition,function(x){x$eigenvalues})
    sigma2 <- lapply(fpca_condition,function(x){x$sigma2}) 
    # calculate the sum of selected lambdas and sigma2 which will be used to get the proportion of the lambda later
    # get the sum of all lambda and sigma2
    lambda.sum <- sum(unlist(lambda.psi)) + sum(unlist(sigma2)) 
    # total pve explained by sigma2
    sigma2.pve <- sum(unlist(sigma2)/lambda.sum)
    # the proportion of pve explained by lambda
    lambda.pve <- sort(unlist(lambda.psi)/lambda.sum, decreasing=TRUE)
    # find the rank of the smallest lambda.psi we select
    rank <- which(cumsum(lambda.pve)+sigma2.pve > control$FPCA2.pve)[1]
    # the value of the smallest lambda.psi we select
    threshold <- sort(unlist(lambda.psi), decreasing=TRUE)[rank]
    
    # calculate psi's at each point of Tij
    psi.MF <- lapply(fpca_condition,function(x){as.matrix(as.matrix(x$eigenfunctions)[visitID,])})
    
    ########################################
    # Assign global eigenfunctions
    ########################################
    efuncs <- list()
    index <- 0
    for(k in 1:ncol(phi)){
      efuncs[[k]] <- psi.MF[[k]]%x%phi[,k]
      for(j in 1:ncol(psi.MF[[k]])){
        # if the eigenvalue is large enough, then we would select its eigenfunctions as basis
        if(lambda.psi[[k]][j] >= max(threshold, 0.05*lambda.sum)){ 
          index <- index + 1
          assign(paste("efunc",index,sep=""), efuncs[[k]][,j], envir = parent.frame()) # make it global variable
        }
      }
    }
    
    # select all error terms(sigma2) to build up independent part. In this way, exchangeable and independent model are nested in unspcified model
    efuncs2 <- list()
    index2 <- 0  
    for(k in 1:ncol(phi)){
      efuncs2[[k]] <- rep(phi[,k], times = nrow(Y))
      index2 <- index2 + 1
      assign(paste("efunc2",index2,sep=""), efuncs2[[k]], envir = parent.frame()) ##make it global variable
    }
    
  }#end if corstr == unspecified
  
  if(corstr == "exchangeable"){
    
    ########################################
    # Assign global eigenfunctions
    ########################################
    efuncs <- list()
    index <- 0
    for(k in 1:ncol(phi)){
      efuncs[[k]] <- rep(phi[,k], times = nrow(Y))
      index <- index + 1
      assign(paste("efunc",index,sep=""), efuncs[[k]], envir = parent.frame()) ##make it global variable
    } 

  }#end if corstr == exchangeable

  if(corstr == "independent"){
    
    ########################################
    # Assign global eigenfunctions
    ########################################
    efuncs <- list()
    index <- 0
    for(k in 1:ncol(phi)){
      efuncs[[k]] <- rep(phi[,k], times = nrow(Y))
      index <- index + 1
      assign(paste("efunc", index, sep=""), efuncs[[k]], envir = parent.frame()) ##make it global variable
    } 
  }#end if corstr == independent
  
  ########################################
  #Step 4: Re-estimate model
  ########################################
  
  g <- rep(as.factor(subjID),each=length(s))
  g_vis <- rep(as.factor(paste(subjID,"_",visitID,sep="")),each=length(s)) #indicator for each subjects's visit
  
  model0 <- paste("s(t, by=x",1:length(Cov.fm),", bs=bs,k=k0) ", collapse="+", sep="")
  model <- paste("y ~ s(t, bs=bs,k=k0) +", model0)
  phiForm <- paste("s(g, by = efunc",1:index,", bs='re') ", collapse="+",sep="")
  
  if(corstr == "unspecified"){
    if(index > 0){
      phiForm_vis <- paste("s(g_vis, by = efunc2",1:index2,", bs='re') ", collapse="+",sep="")
      newformula <- paste(model, "+", phiForm, "+", phiForm_vis)
    }else{
      phiForm_vis <- paste("s(g_vis, by = efunc2",1:index2,", bs='re') ", collapse="+",sep="")
      newformula <- paste(model, "+", phiForm_vis)
    }
  }#end if corstr == unspecified
  
  if(corstr == "exchangeable"){
    phiForm_vis <- paste("s(g_vis, by = efunc",1:index,", bs='re') ", collapse="+",sep="")
    newformula <- paste(model, "+", phiForm, "+", phiForm_vis)
  }#end if corstr == exchangeable
  
  if(corstr == "independent"){
    phiForm_vis <- paste("s(g_vis, by = efunc",1:index,", bs='re') ", collapse="+",sep="")
    newformula <- paste(model, "+", phiForm_vis)
  }#end if corstr == independent

  fit <- bam(as.formula(newformula), discrete = TRUE, nthreads = 2)

  ########################################
  #Step 5: Model fit
  ########################################
  
  Y.pred <- predict(fit, type="terms")
  Y.mean <- rowSums(Y.pred[,1:(length(Cov.fm)+1)]) + fit$coefficients[1]
  Y.mean <- vec2mat(c(Y.mean))
  sigma2 <- fit$sig2
  
  # estimation of beta(s) and their stardard error
  plot.data <- {
    pdf(NULL)
    plot <- plot(fit,n=length(s),seWithMean=TRUE, pages=1)
    invisible(dev.off())
  }
  beta.0 <- plot[[1]]$fit + fit$coefficients[1]
  beta.cov <- matrix(unlist(lapply(1:length(Cov.fm), function(x) plot[[x+1]]$fit)),ncol=length(Cov.fm),byrow=FALSE)
  beta = cbind(beta.0, beta.cov)
  # to get stardard error of initial beta(s)
  beta.se <- matrix(unlist(lapply(1:(length(Cov.fm)+1), function(x) plot[[x]]$se/2)),ncol=(length(Cov.fm)+1),byrow=FALSE)
  
  return(list(phi = phi, psi=psi,
              Y.mean = Y.mean, 
              Y.mean.init = Y.mean.init,
              beta = beta,
              beta.init = beta.init,
              beta.se = beta.se,
              beta.init.se = beta.init.se,
              sigma2 = sigma2,
              final.fit = fit,
              corstr = corstr))
}  


