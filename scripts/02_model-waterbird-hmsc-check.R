# check model fit

# libaries

library(Hmsc)
library(tidyverse)
source('scripts/functions/calculate-residuals.R')

# load model

m <- readRDS('outputs/models/mod-spatialRF_final_notshore.rds')
m

# wait to do this for the final models
# check mcmc convergence

mpost <- convertToCodaObject(m)

# plot effective sample size of regression parameters and potential scale reduction factors
# the effective sample size is the sample size controlled for autocorrelation among sequential posterior samples
# psrf = gelman-rubin convergence diagnostic - compares the different chains to each other - close to one - chains are giving similar results
# psrf should be close to one, effective sample sizes high

# beta regression parameters for each species and covariate (species niches)

hist(effectiveSize(mpost$Beta), main="ess(beta)")
effectiveSize(mpost$Beta) # presence-absence models can be notoriously hard to get high eff. sample size with, so opt to visually inspect trace-plots (below)
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
gelman.diag(mpost$Beta, multivariate = F)$psrf

# get trace plots

pdf('outputs/beta-trace-plots_notshore.pdf')
plot(mpost$Beta)
dev.off()

# omega species-species residual covariance matrix parameters (estimated as spatial random effects)
# note conv.diagnostics for omega are usually not as good as for beta parameters
# easier to estimate fixed rather then random effects
# use option 1 for <10 species
# option 2 for >10 species - subsample of 100 randomly selected pairs of spp (avoid excessive computation)

# option 1 (<10 species)

#hist(effectiveSize(mpost$Omega[[1]]), main="ess(omega)")
#hist(gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf, main="psrf(omega)")

# option 2 (>10 species)

ns <- ncol(m$Y)
sppairs <- matrix(sample(x = 1:ns^2, size = 100))
tmp <- mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega <- effectiveSize(tmp)
psrf.omega <- gelman.diag(tmp, multivariate=FALSE)$psrf
hist(ess.omega)
hist(psrf.omega)

# alpha parameter = estimated spatial scale of the random effect 

plot(mpost$Alpha[[1]])
round(summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5, 0.975))[[2]], 2)
# the spatial scale of the random effect - 50% = median value
# confidence intervals cross zero, so not strong evidence for an effect, 
# but mcmc convergence not good- so estimates of this parameter not reliable
# but see Ovaskainen's comments: still best to keep the spatial random effect to account for autocorrelation
# https://github.com/hmsc-r/HMSC/issues/129

# now check structural model assumptions

preds <- computePredictedValues(m)
dim(preds) # large array storing data matrices of predictions from each posterior sample

# get mean predicted value across all matrices for each species, then calculate residual values

preds.mean <- apply(preds, FUN = mean, MARGIN = c(1,2))

# check structural model assumptions with DS residuals

png('outputs/residuals_notshore.png', width = 500, height = 300)
resid(m$Y, preds.mean, plotds = T, family = 'binomial')
dev.off()

# evaluate model fit/explanatory power for each species

spp.fit <- evaluateModelFit(hM = m, predY = preds)
names <- colnames(m$Y)

mfit <- data.frame(spp = names, AUC = spp.fit$AUC, TjurR2 = spp.fit$TjurR2)
mfit

write.csv(mfit, 'outputs/model-fit_notshore.csv', row.names = F)

# evaluate model performance based on cross validation, i.e. predictive power
# increasing the number of folds means that more data is availabel for fitting the model
# which can lead to greater predictive performance (but computationally intensive to do, so don't do LOO-cross validation)

partition <- createPartition(m, nfolds = 5)
predsCV1 <- computePredictedValues(m, partition=partition)
MF <- evaluateModelFit(hM=m, predY=predsCV1)
MF

cv <- data.frame(spp = names, AUC = MF$AUC, TjurR2 = MF$TjurR2)

write.csv(cv, 'outputs/cross-validation_pred-power.csv', row.names = F)

# conditional cross validation

partition <- createPartition(m, nfolds = 5)
predsCV1 <- computePredictedValues(m, partition=partition, partition.sp(c(1:90)), mcmcStep = 100)
MF <- evaluateModelFit(hM=m, predY=predsCV1)
MF

cv <- data.frame(spp = names, AUC = MF$AUC, TjurR2 = MF$TjurR2)

write.csv(cv, 'outputs/cond-cross-validation_pred-power.csv', row.names = F)
