#####################################################################################
# RUN SINGLE SPECIES OCCUPANCY MODELS
# from Suraci et al. 2021. Global Change Biology
#
# This script prepares the data for and fits single species occupancy models
# in which observations at each camera site are separated into seasonal sampling
# periods of at most 6 months. Human presence and footprint variables are centered 
# and scaled prior to species level dataset preparation such that values are standardized 
# across all species. Project-level random intercepts and slopes from human disturbance, 
# and human disturbance polynomial terms are fit for both occupancy and detection probability.  
####################################################################################
library(tidyverse)
library(rstan) 
library(loo)

# Get occupancy model data
covs <- read.csv(file = "Suraci_etal_2021_GCB_CovariateData.csv", header = T) # Model covariates
spProj <- read.csv(file = "Suraci_etal_2021_GCB_SpeciesByProject.csv", header = T) # project x species matrix
detections <- read.csv(file = "Suraci_etal_2021_GCB_DetectionData.csv", header = T) # Detections
# Remove camera sites with missing HFP values
covs<-covs[which(is.na(covs$HFP_1k)==F),]

# Read in model settings csv
mset<-read.csv("SpeciesModSettings.csv", header = T, colClasses = c('character', 'numeric','numeric',rep("character",3),rep("numeric",3),'character'))

# Get date (for file naming)
date<-unlist(str_split(str_sub(as.character(Sys.time()), start = 1, end = 10), "-"))
date<-paste(date[1],date[2],date[3], sep = "")

# Create scale2 helper function
scale2<-function(x, Log = NA, SD = 1){
  if(is.na(Log)==F) x = log(x + Log)
  scale.x = (x - mean(x)) / (SD * sd(x))
  return(scale.x)
}

# ** SCALE DISTURBANCE COVARIATES **
# Scale human presence and human footprint OUTSIDE of the individual species models 
# so that the range of values is comparable across species
covs$human_scale<-scale2(covs$dr_human, Log = 0.001, SD = 2)
covs$FP_scale<-scale2(covs$HFP_1k, Log = 1, SD = 2)

# Define function packaging occupancy model data for use in Stan
# randSub argument allows for taking a random subset (of specified size) 
# of data for model testing. Ignored if not specified. 
# Returns a named list of length 3: 'modDat' = data for Stan model,
# 'camDat' = camera-level covariates, 'projDat' = project-level covariates.
occData<-function(species, order, randSub = NA){
  
  #______________________
  # Define data set
  # Determine which rows to use based on projects in which species was detected
  pkeep<-spProj$Project[which(spProj[,species]==1)] # Projects where species was detected on at least one camera
  rkeep<-which(covs$Project %in% pkeep & is.na(covs$Baited)==F) # Rows to keep for all relevant data sets
  
  # Take random subset, if specified
  if(is.na(randSub)==F) rkeep <- sample(rkeep, size = randSub, replace = FALSE)
  
  # Get covariate data set truncated to applicable rows
  covData<-covs[rkeep,]; covData<-droplevels(covData) 
  
  #________________
  # Organize y data (observations)
  obData<-as.data.frame(detections[rkeep,]) # get observations data set truncated to applicable rows
  survey.dat<-obData$Weeks # Get number of trials (weeks operational)
  y.dat<-obData[,which(names(obData)==species)] #Extract observations for species (defined in function call)
  
  #________________
  # Organize covariates
  # Get covariates for Observation-level models
  hum.dat<-covData$human_scale # Human detections for each sampling period - SCALED OUTSIDE OF MODEL DATASET
  hum.dat2<-hum.dat^2
  hfp.dat<-covData$FP_scale # Footprint within 1km - SCALED OUTSIDE OF MODEL DATASET
  hfp.dat2<-hfp.dat^2
  
  # drprey<-pull(covData, paste("dr_",prey,"Prey", sep = "")) # Prey detections for each sampling period
  forest100<-covData$Forest_100
  forest1k<-covData$Forest_1k
  NPP<-covData$NPP_1km
  Lat<-covData$Lat
  Long<-covData$Long
  Baited<-covData$Baited %>% 
    scale(., center = T, scale = F)
  isSummer<-ifelse(covData$Season=='Summer',1,0) %>% 
    scale(., center = T, scale = F)
  # "No" and "Partial" hunting combined
  if(order == 'Carnivora') hunt <- covData$Carn_hunted else hunt <- covData$Ung_hunted
  Hunt.dat<- hunt %>% 
    scale(., center = T, scale = F)
  
  # Data transformations
  forest100.dat<-scale2(forest100, SD = 2)
  forest1k.dat<-scale2(forest1k, SD = 2)
  NPP.dat<-scale2(NPP, SD = 2)
  Lat.dat<-scale2(Lat, SD = 2)
  Long.dat<-scale2(Long, SD = 2)
  
  # Prepare covariate matrix for project-level varying coefficients (n x Kp matrix). 
  #X[,1] = 1 (for intercept); X[,2] = Human Pres; X[,3] = Human Footprint 
  Xp<-cbind(rep(1, length(hum.dat)), hum.dat, hfp.dat, hum.dat2, hfp.dat2)
  
  # Prepare covariate matrix for fixed coefficients. Differs for occupancy (n x Kfa) and detection (n x Kfb)
  Xfa<-cbind(Hunt.dat,forest1k.dat, NPP.dat, isSummer, Long.dat, Lat.dat)
  Xfb<-cbind(Hunt.dat, forest100.dat, Baited, isSummer)
  colnames(Xfa) <- c('Hunt','Forest_1k','NPP','Summer','Long','Lat')
  colnames(Xfb) <- c('Hunt','Forest_100m','Baited','Summer')
  
  #________________
  # Create indexing variables
  # projC = project index number for each unique camera site
  # projO = project index number for each unique observation
  # camO = camera site index number for each unique observation
  J<-length(pkeep) # Total number of projects
  C<-length(unique(covData$CamID)) # Total number of camera sites
  projO<-as.numeric(as.factor(covData$Project))
  camO<-as.numeric(as.factor(covData$CamID))
  proj_cams<-tapply(covData$CamID, covData$Project, function(x)length(unique(x))) # Number of cameras per project
  projC<-rep(1, proj_cams[1]) # Create group index variable
  for(i in 2:J){
    projC<-append(projC, rep(i, proj_cams[i]))
  }
  
  #________________
  # Get data for Project-level variables
  # Used in plotting, but not analysis
  
  # Extract data for appropriate projects
  pvars<-proj.vars[proj.vars$Project %in% pkeep,]
  
  # Prep covariates
  hum.p<-pvars$HPD_scale
  hfp.p<-pvars$HFP_scale
  
  # Combine all Project-level data into a J x L matrix. 
  Ua<-cbind(rep(1, dim(pvars)[1]), hum.p, hfp.p)
  Ub<-cbind(rep(1, dim(pvars)[1]), hum.p, hfp.p)
  
  # Specify other key values
  Q<-ifelse(y.dat>0,1,0) # Naive occupancy - for determining likelihood incrementation
  N<-length(y.dat) # Number of observations (each sampling period at each camera site)
  Kp<-dim(Xp)[2] # Number of project-level coefficients on observation-level data
  Kfa<-dim(Xfa)[2] # Number of fixed coefficients on observation-level data for occupancy
  Kfb<-dim(Xfb)[2] # Number of fixed coefficients on observation-level data for detection
  La<-dim(Ua)[2] # Number of covariates on project-level coefficients for occupancy
  Lb<-dim(Ub)[2] # Number of covariates on project-level coefficients for detection
  
  #________________
  # Compile data for Stan
  data<-list(y = y.dat, S = survey.dat, Xp = Xp, Xfa = Xfa, Xfb = Xfb, Ua = Ua, Ub = Ub, J = J, Kp = Kp, Kfa = Kfa, Kfb = Kfb, N = N, C = C, La = La, Lb = Lb, projO = projO, projC = projC, camO = camO, Q=Q)
  
  # Bundle data and results to return
  ret.list<-list(data, covData, pvars)
  ret.list<-setNames(ret.list, c('modDat','camDat','projDat'))
  return(ret.list)
}

#-----------------------------------------------------------------------
# Set model run parameters
chains = 3
cores = chains

# Model and parameters to use
mod.file = "OccMod_BetaBinom.stan"
pars = c('A','B','A_mu','A_sigma','B_mu','B_sigma','a0','b0','a','b','log_lik','mu','rho','psi')
#-----------------------------------------------------------------------
# Fit Beta Binomial Random Effects model for all species

# Loop through species and fit models
for(i in 1:nrow(mset)){
  # Define species data set to use
  sp.dat<-occData(species = mset$SpeciesID[i], order = mset$Order[i])
  modData <- sp.dat$modDat
  
  # Model settings
  iter = mset$iter[i]
  warm = mset$warm[i]
  init = mset$init[i]
  delta = mset$delta[i]
  
  #--------------------------
  #--------------------------
  # Fit model with Camera Intercept
  
  # Monitor start time for each species
  t1<-Sys.time(); print(t1)
  print(paste(mset$SpeciesID[i], "2 START", sep = " "))
  
  # Fit model
  mod<-stan(file = mod.file,
            data = modData,
            pars = pars,
            iter = iter,
            warmup = warm,
            init = init,
            chains = chains,
            cores = cores,
            control = list(adapt_delta = delta))
  
  # Monitor end time for each species
  print(paste(mset$SpeciesID[i], "2 END", sep = " "))
  t2<-Sys.time(); print(difftime(t2,t1))
  
  # Rename model and data objects and save
  modName = paste(mset$SpeciesID[i], "BB_REsq", sep = "")
  datName = paste(mset$SpeciesID[i], ".dat", sep = "")
  fileName = paste(mset$SpeciesID[i], "BB_REsq_", date, ".rda", sep = "")
  assign(modName, mod)
  assign(datName, sp.dat)
  save(list = c(modName, datName), file = paste('ouptut/',fileName))
}


#-----------------------------------------------------------------------
# EXTRACT AND SAVE STAN MODEL OUTPUT
# Extracts and saves Stan model output, along with all associated data sets
# used to fit model, in a convenient form for additional analysis and visualization

modExtract<-function(species, modType, dataPath, plotPath){
  
  # Define model to use
  mod1 <- get(paste(species,modType, sep = ""))
  # Define accompanying data set
  dat<-get(paste(species,".dat", sep = ""))
  
  #-----------------------------------------------------------------
  #-----------------------------------------------------------------
  # EXTRACT MODEL POSTERIORS AND SAVE ALL DATA
  mo<-rstan::extract(mod1)
  mo$A_mean <- apply(mo$A, c(1,3), mean)
  mo$B_mean <- apply(mo$B, c(1,3), mean)
  
  # Add posteriors to data sets (as "mout") and save
  dat$mout<-mo # add model posteriors to data sets list
  assign(x = paste(species,"Mod",sep = ""), dat) # rename datasets list with species name
  # Save
  date<-str_sub(as.character(Sys.time()), start = 1, end = 10)
  dFile = paste("output/", species, modType, "_", date, ".rda", sep = "")
  save(list = c(paste(species,"Mod",sep = "")), file = dFile)

}

# Run for all models
path = "output/" # Set path to retrieve saved model objects
modFiles<-list.files(path = path)
modFiles<-modFiles[grep('BB_', modFiles)]

for(i in 1:length(modFiles)){
  # Load mod file
  load(paste(path,modFiles[i], sep = "/"), verbose = T)
  # Get species and model type
  fn<-strsplit(modFiles[i], split = c("_20"))[[1]][1]
  spname <-strsplit(fn, split = "BB_")[[1]][1]
  mtype<-strsplit(fn, split = spname)[[1]][2]
  
  # Run modExtract
  modExtract(species = spname, modType = mtype, dataPath = "output/BetaBinom_AllSp_Final/")
  
  # Delete model
  rm(list = c(paste(spname,mtype, sep = ""), paste(spname,".dat", sep = "")))
}