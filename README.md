# Suraci-etal-MammalLifeHistory
This repo containts all code necessary for reproducing the analyses in Suraci et al. 2021 "Disturbance type and species life history predict mammal responses to humans," Global Change Biology.  See below for a description of each file included here. All data necessary for replicating the analyses in the above publication area available on figshare at https://doi.org/10.6084/m9.figshare.14444600.v1

__ANALYSIS SCRIPTS__
Located in the _scripts_ directory

_RunOccupancyModels.R_ - Prepares data and fits single species occupancy models (via call to Stan) for each species.  The modExtract function produces convenient lists of data sets from each fitted model object for use in further analyses and plotting.

_CheckOccupancyModels.R_ - Function to produce Bayesian p-values and run other model checking routines (e.g., traceplots, divergence plots) for all single-species occupancy models.

_RunTraitModels.R_ - Prepares data and runs (1) mammal trait PCA and (2) mammal response vs. trait linear mixed effects models (via call to Stan). Requires first running single species occupancy models


__STAN PROGRAMS__
Located in the _scripts_ directory

_OccMod_BetaBinom.stan_ - Stan program fit single species, beta-binomial occupancy model including population- and camarea site-level random effects.  Called in  "RunOccupancyModels.R" script.

_RunTraitModels.R_ - Stan program to fit mammal response vs. trait models (linear mixed effects models with observation error). Called in "RunTraitModels.R" script.


__HELPER FILES__

_SpeciesModSettings.csv_ - Species-specific Stan model settings used in "RunOccupancyModels.R" script.

_SpeciesTraits.csv_ - Mammal species trait values used in fitting the trait PCA. Used in "RunTraitModels.R" script.
