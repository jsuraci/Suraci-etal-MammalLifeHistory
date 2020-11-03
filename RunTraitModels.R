#####################################################################################
#
#                                  NA Cam Trap
#                           Trait Response Model
#                                  2020-07-17
#
# Project level responses by each species to human presence (A1, B1) and human
# footprint (A2, B2) are modeled as functions of species traits.  
# Observation error (estimated as sd of all project level response
# coefficient estimates for a given species) is incorporated for each observed
# response estimate.
# NOTE: Model selection is inconclusive. USING ONLY MODEL PC2  FOR ALL RESPONSES
####################################################################################
library(tidyverse)
library(rstan)
library(boot)
library(lme4)
library(loo)
library(bayesplot)
library(vegan)

# Load data
load(file = "Data_Sets/NACam_OccModDatasets_AllSp_BP_500m_2020721.rda", verbose = T)
# Get mean lat and lon for each project
projLL<- dCovs %>% 
  dplyr::select(Project, CamID, Lat, Long) %>% 
  mutate(Project = as.character(Project)) %>% 
  unique() %>% 
  group_by(Project) %>% 
  summarise(Lat_mean = mean(Lat), Long_mean = mean(Long)) %>% 
  arrange(Project) %>% as.data.frame()



# Read in all model output files 
# ***** Created in RunOccupancyModels.R Script *****
path = "Model_output/BetaBinom_AllSp_Final/"
modFiles<-list.files(path = path, full.names = T)
for(i in 1:length(modFiles)) load(modFiles[i], verbose = T)

# Define species to use (Artiodactyla and Carnivora)
spn<-list.files(path = path)
spnames<-str_split(spn, pattern = "BB_", simplify = T)[,1]
nsp<-length(spnames)

# Helper functions
#####
# Preform some data transformations
scale2<-function(x, Log = NA, SD = 1){
  if(is.na(Log)==F) x = log(x + Log)
  scale.x = (x - mean(x)) / (SD * sd(x))
  return(scale.x)
}
cri<-function(x){
  m<-mean(x)
  ul<-quantile(x, probs = c(0.025,0.975))
  return(c(ul[1], m, ul[2]))
}
cri95<-function(x){
  m<-mean(x)
  ul<-quantile(x, probs = c(0.025,0.975))
  return(c(m, ul[1], ul[2]))
} # As above, but with mean first (mean, lower CI, upper CI)
cri90<-function(x){
  m<-mean(x)
  ul<-quantile(x, probs = c(0.05,0.95))
  return(c(m, ul[1], ul[2]))
} # As above, but with mean first (mean, lower CI, upper CI)
cri80<-function(x){
  ul<-quantile(x, probs = c(0.1,0.9))
  return(c(ul[1], ul[2]))
} # Calculates mean and 80% cri of input
skew <- function(x) {
  xdev <- x - mean(x)
  n <- length(x)
  r <- sum(xdev^3) / sum(xdev^2)^1.5
  return(r * sqrt(n) * (1 - 1/n)^1.5)
}
SSE<-function(x) sum((x-mean(x))^2)
cvar<-function(x) sd(x)/mean(x)
#####

#________________________________________________________________________
#________________________________________________________________________
# RUN TRAIT PCA
#####
# Get species traits
traits<-read.csv("SpeciesTraits.csv", header = T)

# Recode Activity cycle as nocturnal or not
traits$Nocturnal <- ifelse(traits$ActivityCycle==1,1,0)

# Calculate diet diversity (Using shannon index on elton traits, a la Santini et al. 2019)
DietDiv<-traits %>% 
  dplyr::select(Diet.Inv:Diet.PlantO) %>% 
  diversity(index = 'shannon')
traits$DietDiv<-DietDiv

# Select which traits to join to occ mod data
traitsSelect <- traits %>% 
  filter(Species %in% spnames) %>% 
  arrange(Species) %>% 
  dplyr::select(Species, Order, Family, Genus, 
                AdultBodyMass_g, ActivityCycle, WeaningAge_d, 
                LitterSize, SexualMaturityAge_d, MaxLongevity_m, Nocturnal,
                HomeRange_km2, StrictCarn, VertDiet, DietDiv, TrophicLevel)

# Collect variables for PCA
spt<-traitsSelect %>% 
  dplyr::select(AdultBodyMass_g, HomeRange_km2, SexualMaturityAge_d, 
                LitterSize, MaxLongevity_m, WeaningAge_d, VertDiet, DietDiv) %>% 
  mutate(AdultBodyMass_g = log(AdultBodyMass_g), HomeRange_km2 = log(HomeRange_km2), WeaningAge_d = log(WeaningAge_d))

# Run PCA
spPCA<-prcomp(spt, center = T, scale. = T)


## Quick biplot
biplot(spPCA, col = c('white','red'))
text(spPCA$x[,1], spPCA$x[,2], spnames)

# Add first two pc's to traitsSelect
traitsSelect$PC1<-spPCA$x[,1]
traitsSelect$PC2<-spPCA$x[,2]
traitsSelect$PC3<-spPCA$x[,3]

#####

#________________________________________________________________________
#________________________________________________________________________
# PREP DATA 
#####

# Extract means and 95% CIs for each coefficient
respSp<-data.frame()
respProj<-data.frame()
for(i in 1:length(spnames)){
  # Get model coefficients for species i
  mod<-get(paste(spnames[i], "Mod", sep = ""))$mout
  # Get and simplify covariates
  covs<-get(paste(spnames[i], "Mod", sep = ""))$projDat
  covs<-dplyr::select(covs, Project)
  
  # Extract project level coefficient mean and sd
  am<-apply(mod$A,c(2,3),mean)
  as<-apply(mod$A,c(2,3),sd)
  bm<-apply(mod$B,c(2,3),mean)
  bs<-apply(mod$B,c(2,3),sd)
  
  # Extract average coefficients (across projects) mean and sd
  Am<-apply(mod$A_mean,2,mean)
  As<-apply(mod$A,3,sd) # Uses variation across all MCMC runs for all projects
  Asig<-apply(mod$A_sigma, 2, mean)
  Bm<-apply(mod$B_mean,2,mean)
  Bs<-apply(mod$B,3,sd) # Uses variation across all MCMC runs for all projects
  Bsig<-apply(mod$B_sigma, 2, mean)
  
  # Combine into a data frame, with species name and project
  tempP1<-data.frame(Species = rep(spnames[i], nrow(covs)), Project = covs$Project) # Species name as first column of each data frame
  tempP2<-cbind(tempP1, am[,1:3], as[,1:3], bm[,1:3], bs[,1:3])
  
  tempS1<-data.frame(Species = spnames[i]) # Species name as first column of each data frame
  tempS2<-cbind(tempS1, t(Am[1:3]), t(Asig[1:3]), t(As[1:3]), t(Bm[1:3]), t(Bsig[1:3]), t(Bs[1:3]))
  
  # Rename columns
  colnames(tempP2)[3:ncol(tempP2)] = c('A0_mu','A1_mu','A2_mu','A0_sd','A1_sd','A2_sd','B0_mu','B1_mu','B2_mu','B0_sd','B1_sd','B2_sd')
  colnames(tempS2)[2:ncol(tempS2)] = c('A0_mu','A1_mu','A2_mu','A0_sigma','A1_sigma','A2_sigma','A0_sd','A1_sd','A2_sd','B0_mu','B1_mu','B2_mu','B0_sigma','B1_sigma','B2_sigma','B0_sd','B1_sd','B2_sd')
  
  # Add to data frame for all species
  respProj<-rbind(respProj, tempP2)
  respSp<-rbind(respSp, tempS2)
}

# Add species traits (and project-level covariates for respProj)
respProj<-left_join(respProj, traitsSelect)
respProj<-left_join(respProj, projLL)
respSp<-left_join(respSp, traitsSelect)

# Add numeric code for each project, family, species and taxonomic nesting
respProj <- respProj %>% 
  arrange(Order, Family, Species) %>% 
  droplevels()
respProj$projO<-as.numeric(as.factor(respProj$Project))
respProj$spO<-as.numeric(respProj$Species)
respProj$famO<-as.numeric(respProj$Family)
respProj$ordO<-as.numeric(respProj$Order)
# Species nested in family
famSp<-rep(NA, length(unique(respProj$spO)))
for(i in 1:length(famSp)){
  famSp[i]<-unique(respProj$famO[respProj$spO==i])
}
# Family nested in order
ordFam<-rep(NA, length(unique(respProj$famO)))
for(i in 1:length(ordFam)){
  ordFam[i]<-unique(respProj$ordO[respProj$famO==i])
}
#####


#________________________________________________________________________
#________________________________________________________________________
# RUN CANDIDATE MODELS

# Prep model covariates
#####
# Define covariates
PC1 <- respProj$PC1
PC2 <- respProj$PC2
TL2 = ifelse(respProj$TrophicLevel == 2, 1, 0) %>% 
  scale(., center = T, scale = F) %>% .[,1] # + = Omnivore
TL3 = ifelse(respProj$TrophicLevel == 3, 1, 0) %>% 
  scale(., center = T, scale = F) %>% .[,1] # + = Carnivore
OR = scale(respProj$ordO, center = T, scale = F)[,1]
LAT = scale2(respProj$Lat_mean, SD = 2)
LON = scale2(respProj$Long_mean, SD = 2)

# Interactions
PC1xTL2 <- PC1 * TL2
PC1xTL3 <- PC1 * TL3
PC2xTL2 <- PC2 * TL2
PC2xTL3 <- PC2 * TL3
PC1xPC2 <- PC1 * PC2
PC1xOR <- PC1 * OR
PC2xOR <- PC2 * OR

# Polynomials
PC1sq <- PC1^2
PC2sq <- PC2^2

#####

#_____________________
# Data prep functions for each candidate model
trDat_PCM<-function(respPar = c('A1','A2','B1','B2'), obsErr = c('y','sp')){
  # Prep X data
  Xr = cbind(rep(1, nrow(respProj)))
  Xf = cbind(PC1, PC2, PC1xPC2, LAT, LON)
  
  # Prep obs error data (with appropriate alignment to each observation using spO)
  if(respPar == "A1"){
    mu_sp<-respSp$A1_mu[respProj$spO]
    sigma_sp<-respSp$A1_sd[respProj$spO]
    sigma_y<-respProj$A1_sd
    y = respProj$A1_mu
  }
  if(respPar == "A2"){
    mu_sp<-respSp$A2_mu[respProj$spO]
    sigma_sp<-respSp$A2_sd[respProj$spO]
    sigma_y<-respProj$A2_sd
    y = respProj$A2_mu
  }
  if(respPar == "B1"){
    mu_sp<-respSp$B1_mu[respProj$spO]
    sigma_sp<-respSp$B1_sd[respProj$spO]
    sigma_y<-respProj$B1_sd
    y = respProj$B1_mu
  }
  if(respPar == "B2"){
    mu_sp<-respSp$B2_mu[respProj$spO]
    sigma_sp<-respSp$B2_sd[respProj$spO]
    sigma_y<-respProj$B2_sd
    y = respProj$B2_mu
  }
  
  # Choose observation error type (species level or observation level)
  if(obsErr == 'sp') sigma_o <- sigma_sp else sigma_o <- sigma_y
  
  # Organize data
  modDat = list(y = y,
                mu_sp = mu_sp,
                sigma_o = sigma_o,
                Xr = Xr,
                Xf = Xf,
                spO = respProj$spO,
                projO = respProj$projO,
                famO = respProj$famO,
                ordO = respProj$ordO,
                famSp = famSp,
                ordFam = ordFam,
                N = nrow(respProj),
                J = max(respProj$projO),
                S = max(respProj$spO),
                Fa = max(respProj$famO),
                Or = max(respProj$ordO),
                Kr = ncol(Xr),
                Kf = ncol(Xf)
  )
  
  
  
  # Return
  return(modDat)
} 

#_____________________
# MODEL SETTINGS

# Model specs
iter = 4000
warm = 3000
chains = 3
cores = chains
delta = 0.99
tdepth = 15

# Model to use
modFile = "TraitResponse_FINAL_FamSpRE.stan"

# Parameters to save (Depends on model choice above)
pars = c('Bf_mu','Bf_sigma','Bf','Bs_sigma','Bs', 'b','sigma_p', 'z', 'y_new', 'cvar_new', 'log_lik')

#_____________________
# FIT MODELS
#####
ynames<-c('A1','A2','B1','B2')
mnames<-c('PCM')

saveList<-c()

# Loop through each response and model time
# Fit model and store model, output and data in a list
for(i in 1:length(ynames)){
  for(j in 1:length(mnames)){
    
    # Get model names
    modName<-paste(ynames[i], mnames[j], sep = "_")
    
    # Monitor loop progress
    print(modName)
    t1<-Sys.time()
    
    # Get applicable data generating function
    df<-get(paste("trDat_",mnames[j], sep = ""))
    # Generate data and fit model
    dat<-df(respPar = ynames[i], obsErr = 'y')
    mod<-stan(file = modFile,
              data = dat,
              pars = pars,
              iter = iter,
              warmup = warm,
              chains = chains,
              cores = cores,
              control = list(adapt_delta = delta,
                             max_treedepth = tdepth))
    (t2<-Sys.time()); print(difftime(t2,t1))
    # Gather into list for saving
    mout<-rstan::extract(mod)
    modList<-list(dat,mod,mout)
    modList<-setNames(modList, c('modDat','mod','mout'))
    # Assign correct model name
    assign(modName, modList)
    
    # Add model name to save list
    saveList<-c(saveList, modName)
    
  }
}

#####


# Save everything
save(list = c('A1_PCM','A2_PCM','B1_PCM','B2_PCM'), file = "TraitResponseModels.rda")
