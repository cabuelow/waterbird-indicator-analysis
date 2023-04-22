# joint species distribution model (JSDM)

# libaries

library(Hmsc)
library(MASS)
library(tidyverse)
library(sf)
library(raster)
library(GGally)

# data

birds.env <- read.csv('data/waterbird-basin-threats.csv')

# check multi-collinearity in predictor variables

colnames(birds.env)
predictors.cor1 <- cor(dplyr::select(birds.env, Consumptive_Water_Loss,
                                    Pesticide_Loading,
                                    Phosphorus_Loading, Nitrogen_Loading,
                                    Aquaculture_Pressure, Mercury_Deposition,
                                    ire_pc_sse, pre_mm_syr))
symnum(predictors.cor1)

predictors.cor2 <- cor(dplyr::select(birds.env, Consumptive_Water_Loss,
                                     Cropland,
                                     Phosphorus_Loading, Nitrogen_Loading,
                                     Aquaculture_Pressure, Mercury_Deposition, 
                                     ire_pc_sse, pre_mm_syr))
symnum(predictors.cor2)

# plot individual corr coefficients

ggcorr(dplyr::select(birds.env, Consumptive_Water_Loss,
                     Pesticide_Loading,
                     Phosphorus_Loading, Nitrogen_Loading,
                     Aquaculture_Pressure, Mercury_Deposition, 
                     ire_pc_sse, pre_mm_syr), label = TRUE)

# pesticide loading and cropland are highly correlated, but important threats
# run two models, one with each, and choose model with highest probability of effect

pred1 <- birds.env %>% dplyr::select(Consumptive_Water_Loss,
                                     Pesticide_Loading, 
                                     Phosphorus_Loading, Nitrogen_Loading,
                                     Aquaculture_Pressure, Mercury_Deposition, 
                                     ire_pc_sse, pre_mm_syr)

pred2 <- birds.env %>% dplyr::select(Consumptive_Water_Loss,
                                     Cropland,
                                     Phosphorus_Loading, Nitrogen_Loading,
                                     Aquaculture_Pressure, Mercury_Deposition,
                                     ire_pc_sse, pre_mm_syr)
# specify hmsc model

xycoords <- as.matrix(birds.env %>% dplyr::select(X,Y))
colnames(xycoords) <- c('x-coordinate', 'y-coordinate')
row.names(xycoords) <- 1:nrow(birds.env)
y <- as.matrix(birds.env %>% dplyr::select(Australasian.Darter:Australasian.Bittern))
write.csv(y, 'outputs/mod-site-basin-raneff_y.csv', row.names = F) # save to use for model diagnostic/plotting scripts
x <- pred1
write.csv(x, 'Outputs/mod-site-basin-raneff_x_pesticide.csv', row.names = F) # save to use for model diagnostic/plotting scripts
#x <- pred2
#write.csv(x, 'Outputs/mod-site-basin-raneff_x_cropland.csv', row.names = F) # save to use for model diagnostic/plotting scripts
formula <- paste0(colnames(x), collapse = '+')

# random effect at sampling unit to use latent factors to model spp associations
#studyDesign = data.frame(sample = as.factor(1:nrow(birds.env)))
#rL = HmscRandomLevel(units = studyDesign$sample)

# if setting up a spatial random effect
studyDesign <- data.frame(sample = as.factor(1:nrow(birds.env)))
rL.spatial <- HmscRandomLevel(sData = xycoords)
rL.spatial <- setPriors(rL.spatial, nfMin=3, nfMax=5) # here set number of LVs (3-5)

m <- Hmsc(Y = y, XData = x, 
          XFormula = ~Consumptive_Water_Loss+Pesticide_Loading+Phosphorus_Loading+Nitrogen_Loading+Aquaculture_Pressure+Mercury_Deposition+ire_pc_sse+pre_mm_syr,
          studyDesign = studyDesign, 
          ranLevels = list(sample = rL.spatial), 
          distr = 'probit')

# set mcmc sampling parameters
# takes 4 hours to run!!

nChains <- 4
test.run <- FALSE
if (test.run){
  # with this option, the vignette runs fast but results are not reliable
  thin = 1
  samples = 10
  transient = 5
  verbose = 0
} else {
  # with this option, the vignette evaluates slow but it reproduces the results of
  # the .pdf version
  thin = 100
  samples = 1000
  transient = ceiling(0.5*samples*thin)
  verbose = 1 # set to 0 to suppress output reporting how MCMC sampling is proceeding
}

# run model
set.seed(123)

system.time( # takes 20 mins to run for 50 spp, 21 covars (1000 samples, 4 chains, 100 thinning, transient = burnin)
m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient, adaptNf = rep(ceiling(0.4*samples*thin),1),
                       nChains = nChains, nParallel = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
)
# NOTE: adaptNf controls the number of interations during which the number of latent factors are adapted

# save model

saveRDS(m, 'outputs/models/mod-pesticide_spatialRF.rds')
#saveRDS(m, 'outputs/models/mod-non-spatial_cropland.rds')


