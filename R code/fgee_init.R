fgee_init <- function(formula, Y, Cov, s){
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
  ## additional setups
  k0 <- 10 ##number of basis functions for fixed effects functions
  bs <- "ps"
  for(i in 1:length(Cov.fm)){
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
  Y.resid <- as.matrix(Y - Y.mean.init)
  
  #initial estimation of beta(s) and their stardard error
  plot.data <- {
    pdf(NULL)
    plot.init <- plot(fit_init, n=length(s), seWithMean=TRUE, pages=1)
    invisible(dev.off())
  }
  
  beta.init.0 <- plot.init[[1]]$fit + fit_init$coefficients[1]
  beta.init.cov <- matrix(unlist(lapply(1:length(Cov.fm), function(x) plot.init[[x+1]]$fit)),ncol=length(Cov.fm),byrow=FALSE)
  beta.init = cbind(beta.init.0, beta.init.cov)
  #To get stardard error of initial beta(s)
  beta.init.se <- matrix(unlist(lapply(1:(length(Cov.fm)+1), function(x) plot.init[[x]]$se/2)),ncol=(length(Cov.fm)+1),byrow=FALSE)
  
  return(list(Y.mean.init = Y.mean.init,
              beta.init = beta.init,
              beta.init.se = beta.init.se))
}
