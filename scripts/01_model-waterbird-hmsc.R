# joint species distribution model (JSDM)

# libaries

library(Hmsc)
library(tidyverse)
library(GGally)
library(scales)

# data

birds.env <- read.csv('data/waterbird-basin-threats-100_bioregion_notshore.csv') 

# there was multicollinearity in several of the predictors,
# so we have chosen a set that aren't multicollinear, chosen based on relevance to management

preds <- birds.env %>%
  dplyr::select(Consumptive_Water_Loss, # select predictors only
                                     Pesticide_Loading, 
                                     Phosphorus_Loading, Nitrogen_Loading,
                                     Aquaculture_Pressure,
                                     pre_mm_syr) 
# rescale precipitation from 0 to 1 so on same scale as other predictors
preds$pre_mm_syr <- scales::rescale(preds$pre_mm_syr, to = c(0,1))

# check multi-collinearity in predictor variables (greater than or less than 0.7)

ggcorr(preds, label = TRUE)
ggsave('outputs/corr-predictors.png', width = 10, height = 10)

# use ggpairs to explore data

species <- colnames(dplyr::select(birds.env, Australasian.Darter:Musk.Duck)) # get species names

# loop through and make a ggpairs plot for each species
# to explore relationships between each species presence-absence and predictors

#my_fn <- function(data, mapping, method="loess", ...){
 # p <- ggplot(data = data, mapping = mapping) + 
  #  geom_point() + 
   # geom_smooth(method=method, ...)
  #p
#}

#for(i in 1:length(species)){
#ggpairs(birds.env[,c(i+1, 92:109)], lower = list(continuous = my_fn))
#ggsave(paste0('outputs/spp-plots/', species[i], '_ggpairs-plot.png'), width = 12, height = 10)
#}

# specify hmsc model

xycoords <- as.matrix(birds.env %>% dplyr::select(X,Y))
colnames(xycoords) <- c('x-coordinate', 'y-coordinate')
row.names(xycoords) <- 1:nrow(birds.env)
y <- as.matrix(birds.env %>% dplyr::select(Australasian.Darter:Musk.Duck))
x <- preds

# spatial random effect for each sampling location (basin)

studyDesign <- data.frame(sample = as.factor(1:nrow(birds.env)), bioregion = as.factor(birds.env$Q_REG))
rL.spatial <- HmscRandomLevel(sData = xycoords)
rL.spatial <- setPriors(rL.spatial, nfMin=3, nfMax=5) # here set number of LVs (3-5)
rL <- HmscRandomLevel(units = unique(studyDesign$bioregion))

m <- Hmsc(Y = y, XData = x, 
          XFormula = as.formula(paste('~', paste(colnames(x), collapse = '+'))),
          studyDesign = studyDesign, 
          ranLevels = list('sample' = rL.spatial, 'bioregion' = rL), 
          distr = 'probit')

# set mcmc sampling parameters 

nChains <- 4
test.run <- FALSE

if (test.run == TRUE){
  thin = 1
  samples = 10
  transient = 5
  verbose = 0
} else {
  thin = 10
  samples = 1000
  #transient = ceiling(0.5*samples*thin)
  transient = 1000
  verbose = 1 # set to 0 to suppress output reporting how MCMC sampling is proceeding
}

# run model

set.seed(123)

system.time( # takes 15 mins to run (1000 samples, 4 chains, 10 thinning, 1000 burnin)
  m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient, 
                  nChains = nChains, verbose = verbose, updater=list(GammaEta=FALSE))
)
# save model

saveRDS(m, 'outputs/models/mod-spatialRF_final_notshore.rds')




