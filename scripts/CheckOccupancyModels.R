################################################################
# CHECK OCCUPANCY MODELS
# from Suraci et al. 2021. Global Change Biology
#
# Model checking and Bayesian p-value calculation function 
# for beta binomial occupancy models 
################################################################
library(tidyverse)
library(rstan)
library(loo)
library(bayesplot)
library(boot)
library(updog)

# Some additional pval functions
cvar<-function(x) sd(x)/mean(x)
skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}

# modCheck produces several of the key model checking outputs from the StanModelChecking_BetabBinomial.R file as a list
modCheck<-function(species, modType){
  
  # Definte model to use
  mod1 <- get(paste(species,modType, sep = ""))
  # Define accompanying data set
  dat1 <- get(paste(species,".dat", sep = ""))$modDat  
  dat2 <- get(paste(species,".dat", sep = ""))$projDat
  
  #-----------------------------------------------------------------
  #-----------------------------------------------------------------
  # EXTRACT MODEL POSTERIORS AND SAVE ALL DATA
  mo<-rstan::extract(mod1)
  mo$A_mean <- apply(mo$A, c(1,3), mean)
  mo$B_mean <- apply(mo$B, c(1,3), mean)
  
  #-------------------------
  #-------------------------
  # CONVERGENCE
  #####
  # Extract model summary and check for Rhat values over 1.1
  msum<-summary(mod1)$summary
  rh<-max(msum[,"Rhat"])
  # Check out the trace plots for higher order coefficients obs-level coefficients
  trace<-rstan::traceplot(mod1, pars = c('a','b')) + ggtitle(paste(spname, "Max Rhat =", rh, sep = " "))
  #####
  
  #-------------------------
  #-------------------------
  # POSTERIOR PREDICTIVE DIST AND CHECKS
  #####
  # Simulate new data from model and calculate bayesian p-values
  reps = dim(mo$A)[1]
  N = dat1$N
  J = dat1$J
  y_new <- array(numeric(),c(reps, N))
  z <- array(numeric(),c(reps, N))
  Tobs <- vector('double', length = reps)
  Tnew <- vector('double', length = reps)
  TOproj <- array(numeric(),c(reps, J))
  TNproj <- array(numeric(),c(reps, J))
  CHobs <- vector('double', length = reps)
  CHnew <- vector('double', length = reps)
  COproj <- array(numeric(),c(reps, J))
  CNproj <- array(numeric(),c(reps, J))
  for(i in 1:reps){
    
    # Recover z matrix
    # Likelihood of no detections if psi == 1 divided by the total likelihood (i.e., psi == 1 + likelihood of no detections if psi == 1)
    
    # Calculate parameters of beta-binom dist from up and rho
    alpha = mo$mu[i,] * (1 - mo$rho[i]) / mo$rho[i]
    beta = (1 - mo$mu[i,]) * (1 - mo$rho[i]) / mo$rho[i]
    
    # Get probability of z = 1 when y = 0
    z_prob <- (mo$psi[i,] * dbetabinom(x = rep(0, N), size = dat1$S, mu = mo$mu[i,], rho = mo$rho[i], log = F)) / ((1-mo$psi[i,]) + (mo$psi[i,] * dbetabinom(x = rep(0, N), size = dat1$S, mu = mo$mu[i,], rho = mo$rho[i], log = F)))
    
    z_rand <- rbinom(n = N, size = 1, prob = z_prob)
    z[i,] <- ifelse(dat1$Q == 1, 1, z_rand)
    
    # Create new data
    mu.lam <- mo$mu[i,]*z[i,]
    y_new[i,] <- rbetabinom(n = N, size = dat1$S, mu = mu.lam, rho = mo$rho[i])
    
    
    # BAYESIAN P-VALUES
    #-----------------------
    # Freeman-Tukey Residual
    terr = (sqrt(dat1$y) - sqrt(mu.lam * dat1$S))^2
    terrnew = (sqrt(y_new[i,]) - sqrt(mu.lam * dat1$S))^2
    # Overall model sum
    Tobs[i] = sum(terr)
    Tnew[i] = sum(terrnew)
    # Sum by project
    TOproj[i,] = tapply(terr, dat1$projO, sum)
    TNproj[i,] = tapply(terrnew, dat1$projO, sum)
    
    #-----------------------
    # Chi-squared stat
    cherr = ((dat1$y - (mu.lam * dat1$S))^2) / sqrt(mu.lam * dat1$S + 0.01)
    cherrnew = ((y_new[i,] - (mu.lam * dat1$S))^2) / sqrt(mu.lam * dat1$S + 0.01)
    # Overall model sum
    CHobs[i] = sum(cherr)
    CHnew[i] = sum(cherrnew)
    # Sum by project
    COproj[i,] = tapply(cherr, dat1$projO, sum)
    CNproj[i,] = tapply(cherrnew, dat1$projO, sum)
  }
  
  # Get data.frame of p-values for model specified above
  Tp<-length(which(Tnew >= Tobs))/length(Tnew)
  Cp<-length(which(CHnew >= CHobs))/length(CHnew)
  pv<-data.frame(Data = "FullMod", Tukey = Tp, ChiSq = Cp)
  
  # P-value by project
  for(j in 1:J){
    proj<-dat2$Project[j]
    Tpval<-length(which(TNproj[,j] >= TOproj[,j]))/length(TNproj[,j])
    Cpval<-length(which(CNproj[,j] >= COproj[,j]))/length(CNproj[,j])
    # print(paste(proj, ": p = ", round(Tpval,4), sep = ""))
    
    pv.temp<-data.frame(Data = proj, Tukey = Tpval, ChiSq = Cpval)
    pv<-rbind(pv, pv.temp)
  }
  
  # Rename p-value data frame
  pv[,2:3]<-round(pv[,2:3],4)
  #####
  
  #-------------------------
  #-------------------------
  # Model Checking Visualizations (following Gabry et al 2019 J Roy Stat Soc)
  #####
  
  # Investigating divergences
  theme_update(axis.text = element_text(size = 20))
  hmc_diagnostics <- nuts_params(mod1)
  color_scheme_set("darkgray")
  
  div_style <- scatter_style_np(div_color = "green", div_size = 2.5, div_alpha = 0.75)
  scat <- 
    mcmc_scatter(
      mod1,
      size = 1.5,
      alpha = 2/3,
      pars = c("a[1]", "a[2]"),
      # transform = list(tau1 = "log"), 
      np = hmc_diagnostics,
      np_style = div_style
    )
  scatter<-scat + ggtitle(spname)
  #####
  
  mchecks<-list(trace, scatter, pv)
  return(mchecks)
}

# RUN FUNCTION
# Read in all occupancy model output files 
# ***** Created in RunOccupancyModels.R Script *****
path = "output/BetaBinom_AllSp_Final/"
modFiles<-list.files(path = path)
modFiles<-modFiles[grep('.rda', modFiles)]
tPlots<-list()
sPlots<-list()
pvList<-list()

for(i in 1:length(modFiles)){
  # Load mod file
  load(paste(path,modFiles[i], sep = "/"), verbose = T)
  # Get species and model type
  fn<-strsplit(modFiles[i], split = c("_20"))[[1]][1]
  spname <-strsplit(fn, split = "BB_")[[1]][1]
  mtype<-strsplit(fn, split = spname)[[1]][2]
  
  # Run function
  temp<-modCheck(species = spname, modType = mtype)
  
  tPlots[[i]]<-temp[[1]]
  sPlots[[i]]<-temp[[2]]
  pvList[[i]]<-temp[[3]]
  names(pvList)[[i]]<-spname
  
}

# Print trace plots to file
pdf(file = "output/AllSp_traceplots.pdf")
invisible(lapply(tPlots, print))
dev.off()

# Print divergence plots to file
pdf(file = "output/AllSp_DivergencePlots.pdf")
invisible(lapply(sPlots, print))
dev.off()

# Save p-values
save(pvList, file = "output/AllSp_Pvalues.rda")
