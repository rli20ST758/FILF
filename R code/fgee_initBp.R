fgee_initBp <- function(formula, Y, Cov, s, subjID, Tij,
                        numLongiPoints = 41, n.Bp=300){
  
  #Arguments:
  # The bootstrap algorithm is from Park et al. (2017)
  #1.fomular 
  #parameter fomular: a formula expression which indicates a linear relationship between Y and Cov
  #
  #2.data 
  #parameter Y -- a matrix of which each row corresponds to one curve observed on a regular and dense grid 
  #(dimension of N by m; N = total number of observed functions; m = number of grid points) 
  #parameter Cov -- including all covariates of interest.      
  #parameter s -- argvals a string indicating the functional domain variable name in data   
  #parameter subjID -- subject id; vector of length N with each element corresponding a row of Y
  #parameter Tij -- actual time of visits at which a function is observed
  #n.Bp: the times of resampling
  
  #formula for initial estimate
  model.formula <- as.character(formula)
  stopifnot(model.formula[1] == "~" & length(model.formula) == 3)
  #complie the formula again
  formula <- paste("Y", model.formula[1], model.formula[3])
  #select covariates are used in formula
  Cov.fm <- c() 
  for(i in 1:ncol(Cov)){
    if(grepl(colnames(Cov)[i], model.formula[3])){
      Cov.fm <- c(Cov.fm, i)
    }
  }
  n.cov <- length(Cov.fm)
  subjID.uniq <- unique(subjID)
  
  # additional setups
  k0 <- 10 #number of basis functions for fixed effects functions
  bs <- "ps"
  for(i in 1:n.cov){
    assign(colnames(Cov)[Cov.fm[i]], Cov[,Cov.fm[i]])
  }
  ###########################################
  #Step 1: Initial mean estimate
  ##########################################
  
  #initial mean estimates
  fit_init <- pffr(as.formula(formula), yind=s, algorithm = "bam",
                   bs.yindex = list(bs = bs, k = k0, m = c(2, 2))) ## initial fit
  Y.mean_init <- predict(fit_init, type="terms")
  Y.mean.init <- fitted(fit_init) ## fitted means
  
  #initial estimation of beta(s) 
  plot.data <- {
    pdf(NULL)
    plot.init <- plot(fit_init,n=length(s), seWithMean=TRUE, pages = 1)
    invisible(dev.off())
  }
  beta.init.0 <- plot.init[[1]]$fit + fit_init$coefficients[1]
  beta.init.cov <- matrix(unlist(lapply(1:n.cov, function(x) plot.init[[x+1]]$fit)),ncol=n.cov,byrow=FALSE)
  beta.init = cbind(beta.init.0, beta.init.cov)

  #####################################################################
  #Step 2: To get confidence bands of initial beta(s) by bootstroap
  ####################################################################

  for(i in 1:n.Bp){
    #resample the subject indexa from the index set
    index <- sample(subjID.uniq, size=length(subjID.uniq), replace=TRUE) 
    #find correspending rows to original data for subjects from the resample set 
    Bp.row <- unlist(lapply(1:length(index), function(i){which(subjID==index[i])}))
    Y.Bp <- Y[Bp.row, ]
    Cov.Bp <- Cov[Bp.row, Cov.fm]
    colnames(Cov.Bp) <- paste("X.Bp", 1:n.cov, sep="")
    ## formula for bootstrap estimate
    formula.Bp <- paste("Y.Bp", "~", paste("X.Bp", 1:n.cov, collapse="+", sep=""))
    # bootstrap estimates for beta(s)
    fit.Bp <- pffr(as.formula(formula.Bp), yind=s, data=data.frame(Y.Bp=Y.Bp, Cov.Bp),
                   algorithm = "bam", bs.yindex = list(bs=bs, k=k0, m=c(2, 2))) 
    plot.data <- {
      pdf(NULL)
      plot.Bp <- plot(fit.Bp,n=length(s),seWithMean=TRUE,pages = 1)
      invisible(dev.off())
    }
    
    if(i==1){
      beta.Bp <- lapply(1:(n.cov+1), function(x) plot.Bp[[x]]$fit)
      beta.Bp[[1]] <- beta.Bp[[1]] + fit.Bp$coefficients[1]
    } else {
      for(kk in 1:(n.cov+1)){
        beta.Bp[[kk]] = cbind(beta.Bp[[kk]], plot.Bp[[kk]]$fit)
      }
      beta.Bp[[1]][,i] <- beta.Bp[[1]][,i] + fit.Bp$coefficients[1]
    }
  }
  
  #To get confidence bands of beta(s)
  beta.band.lower <- matrix(unlist(lapply(1:(n.cov+1), function(x){
    apply(beta.Bp[[x]], 1, function(row)  mean(row)-qnorm(0.975)*sd(row) ) })), ncol=(n.cov+1), byrow=FALSE)
  beta.band.upper <- matrix(unlist(lapply(1:(n.cov+1), function(x){
    apply(beta.Bp[[x]], 1, function(row)  mean(row)+qnorm(0.975)*sd(row) ) })), ncol=(n.cov+1), byrow=FALSE)
  
  #To get confidence bands of beta(s)
  #beta.band.lower <- matrix(unlist(lapply(1:(length(Cov.fm)+1), function(x){
  #  apply(beta.Bp[[x]], 1, function(row) quantile(row, 0.025)) })), ncol=(length(Cov.fm)+1), byrow=FALSE)
  #beta.band.upper <- matrix(unlist(lapply(1:(length(Cov.fm)+1), function(x){
  #  apply(beta.Bp[[x]], 1, function(row) quantile(row, 0.975)) })), ncol=(length(Cov.fm)+1), byrow=FALSE)
  
  return(list(Y.mean.init = Y.mean.init,
              beta.init = beta.init,
              beta.Bp = beta.Bp,
              beta.band.lower = beta.band.lower,
              beta.band.upper = beta.band.upper))
}  



