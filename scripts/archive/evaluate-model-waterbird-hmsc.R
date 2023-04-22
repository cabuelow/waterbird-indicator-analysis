# Evaluate and plot JSDM

# libaries

library(Hmsc)
library(MASS)
library(tidyverse)
library(tmap)
library(sf)
library(raster)
library(corrplot)
library(RColorBrewer)
library(randomcoloR)
library(ggcorrplot)
library(abind)
library(cluster)
library(factoextra)

# load model and data

m.spatial <- readRDS('outputs/models/mod-spatial-3-10LVs.rds')
m.spatial

dat <- read.csv('data/waterbird-basin-threats-comp-cases.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')

# extract posterior distribution from model object

mpost <- convertToCodaObject(m.spatial) # extract the posterior distribution

# check structural model assumptions

preds <- computePredictedValues(m.spatial)
dim(preds) # large array storing 2000 data matrices of predictions from each posterior sample
# get mean predicted value across all 2000 matrices for each species, then calulate residual values
preds.mean <- apply(preds, FUN = mean, MARGIN = c(1,2))
# calculate residuals for each species at each site
resids <- m.spatial$Y-preds.mean
# TODO: calculating and checking randomised quantile residuals

# evaluate model fit/explanatory power for each species

evaluateModelFit(hM = m.spatial, predY = preds)

# partition explained variance among fixed and random effects

head(m.spatial$X)
VP <- computeVariancePartitioning(m.spatial, group = c(rep(1,13), rep(2,3)), groupnames = c("river-threat","env-covariate"))
plotVariancePartitioning(m.spatial, VP = VP)

vardf <- data.frame(VP$vals)
vardf$cat <- c('Threats', 'Environmental', 'Spatial')
vardf.long <- vardf %>%
  pivot_longer(Australasian.Darter:Australian.Shelduck, names_to = 'Species', values_to = 'Variance_prop')

ggplot(vardf.long) +
  aes(x = Species, y = Variance_prop, fill = cat) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  ylab('Proportion variance explained') +
  xlab('') +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(3, 'Set2')) +
  theme(legend.title = element_blank())

ggsave('outputs/spp-var-explained.png', width = 8, height = 6)

# evaluate predictive power with cross-validation

partition <- createPartition(m.spatial, nfolds = 2, column = "sample")
cvpreds.spatial <- computePredictedValues(m.spatial, partition=partition,updater=list(GammaEta=FALSE))
cvMF.spatial <- evaluateModelFit(hM=m.spatial, predY=cvpreds.spatial)
cvMF.spatial

# plot estimats of Beta parameters for each spp as heatmap
# red = posterior probability for the parameter being positive is greater than 0.95 
# blue = posterior probability for the parameter being negative is greater than 0.95 
# white = no strong statistical support for parameter being positive or negative

postBeta <- getPostEstimate(m.spatial, parName = "Beta")
#pdf('outputs/spp-regression-params.pdf')
plotBeta(m.spatial, post = postBeta, param = "Support", 
         supportLevel = 0.95, mgp = c(1.5,1.5,0), 
         mar = c(5,0,5,0), cex = c(0.4,0.4,0.5))
#dev.off()

beta.mean <- data.frame(postBeta$mean)
beta.mean$variable <- c('intercept', as.character(colnames(x)))
beta.mean.long <- beta.mean %>%
  pivot_longer(Australasian.Darter:Australian.Shelduck, names_to = 'Species', values_to = 'Beta_mean') %>% 
  filter(variable != 'intercept')

beta.support <- data.frame(postBeta$support)
beta.support$variable <- c('intercept', as.character(colnames(x)))
beta.support.long <- beta.support %>%
  pivot_longer(Australasian.Darter:Australian.Shelduck, names_to = 'Species', values_to = 'Beta_support') %>% 
  filter(variable != 'intercept')

beta.supportNeg <- data.frame(postBeta$supportNeg)
beta.supportNeg$variable <- c('intercept', as.character(colnames(x)))
beta.supportNeg.long <- beta.supportNeg %>%
  pivot_longer(Australasian.Darter:Australian.Shelduck, names_to = 'Species', values_to = 'Beta_supportNeg') %>% 
  filter(variable != 'intercept')

beta.mean.long$support <- beta.support.long$Beta_support
beta.mean.long$supportNeg <- beta.supportNeg.long$Beta_supportNeg

# filter for spp-ind relationships with greater than 75% support

beta.pos <- filter(beta.mean.long, support > 0.75)
beta.neg <- filter(beta.mean.long, supportNeg > 0.75)

beta.posneg <- rbind(beta.pos, beta.neg)

ggplot(beta.mean.long) +
  aes(x = Species, y = Beta_mean, col = variable) +
  geom_point() +
  theme_classic() +
  ylab('Coefficient') +
  xlab('') +
  coord_flip() +
  scale_fill_manual(values = randomColor(16)) +
  theme(legend.title = element_blank())

#ggsave('outputs/spp-beta-coefs2.png', width = 8, height = 6)

# summarise spp coefficents from posterior distribution 
# mean, median, and credible intervals

beta.coefs1 <- data.frame(summary(mpost$Beta[[1]])$statistics) # first chain
beta.quant1 <- data.frame(summary(mpost$Beta[[1]])$quantiles) # first chain
variable <- rep(c('intercept', colnames(x)), nrow(beta.coefs1)/(length(colnames(x))+1))
names <- colnames(y)
spp <- c()
for(i in 1:length(colnames(y))){
  sp <- rep(names[i], nrow(beta.coefs1)/(length(colnames(y))-1))
  spp <- c(spp, sp)
}

beta.coefs1 <- beta.coefs1 %>% 
  mutate(CI.2.5 = beta.quant1$X2.5.,
         CI.97.5 = beta.quant1$X97.5.,
         variable = variable,
         species = spp) %>% 
  filter(variable != 'intercept')

for(i in 1:length(colnames(x))){
  pl <- filter(beta.coefs1, variable == paste(colnames(x)[i]))
  ggplot(pl) +
    aes(x = species, y = Mean) +
    geom_point() +
    geom_errorbar(aes(ymin=CI.2.5, ymax=CI.97.5), width=.1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_fill_manual(values = distinctColorPalette(16)) +
    theme_classic() +
    ylab('Coefficient') +
    xlab('') +
    coord_flip() +
    theme(legend.title = element_blank()) 
  
  ggsave(paste0('outputs/spp-ind-coefs/', colnames(x)[i], '.png'), width = 8, height = 6)
}

# extract estimates of species-species associations
# computeAssociations function converts the covariances to correlations

OmegaCor <- computeAssociations(m.spatial)
corr.mat <- OmegaCor[[1]]$mean
ggcorrplot(corr.mat, tl.cex = 7)
ggsave('outputs/spp-correlation-matrix.png', width = 10, height = 5.5)

# extract species and site loadings on LVs

etaPost <- getPostEstimate(m.spatial, "Eta") # site loadings 
lambdaPost <- getPostEstimate(m.spatial, "Lambda") # spp loadings

# plot variation in spp richness over enviro gradients
# Note: no non-focal variables are listed, so  predictions for the net effect of
# the focal variable, that is the total effect of the focal variable and 
# any non-focal variables covarying with the focal variable

for(i in 1:ncol(x)){
  Gradient <- constructGradient(m.spatial, focalVariable = paste(colnames(x)[i]))
  Gradient$XDataNew
  predY <- predict(m.spatial, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                   ranLevels=Gradient$rLNew, expected=TRUE)
  pdf(paste0('outputs/spp-env-gradients/spp-richness_', colnames(x)[i], '.pdf'), width =9, height =7) 
  plotGradient(m.spatial, Gradient, pred=predY, measure="S", showData = TRUE)
  dev.off()
}

# plot variation over enviro gradients for indicator spp 

spp <- colnames(m.spatial$Y) # index of spp. in model
iter <- c(1, 4, 7, 6, 8, 15, 16, 19, 21, 24, 33, 35, 40)

for(j in iter){
  sp <- spp[j]
for(i in 1:ncol(x)){
  Gradient <- constructGradient(m.spatial, focalVariable = paste(colnames(x)[i]))
  Gradient$XDataNew
  predY <- predict(m.spatial, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                   ranLevels=Gradient$rLNew, expected=TRUE)
  pdf(paste0('outputs/spp-env-gradients/', sp, '_', colnames(x)[i], '.pdf'), width =9, height =7) 
  plotGradient(m.spatial, Gradient, pred=predY, measure="Y", index = j, showData = TRUE)
  dev.off()
}}

# determine regions of similar species composition
#**TODO here look at finnish bird vignette for predicting to new basins without spp data

predY <- predict(m.spatial, XData=m.spatial$XData, studyDesign=m.spatial$studyDesign,
                 ranLevels=m.spatial$rL, expected=TRUE, predictEtaMean = TRUE)
predY2 <- apply(abind(predY,along=3),c(1,2),mean)

fviz_nbclust(predY2, kmeans, method='wss')
spp.comp <- kmeans(predY2, 8)

# dat for plotting

dat2 <- dat %>% 
  select(HYBAS_ID, X, Y)
dat2$cluster <- as.factor(spp.comp$cluster)

# conditional prob occurrence ind species

ind.spp <- predY2[,iter]
dat2 <- cbind(dat2, ind.spp)
dat2$sppRich <- rowSums(predY2)

# convert to basins for plotting

basins <- basins.qld %>% 
  inner_join(dat2, by = 'HYBAS_ID')

# plot

tmap_mode('view')

tm_shape(basins) + 
  tm_polygons(col = 'cluster')

tm_shape(basins) + 
  tm_polygons(col = 'sppRich')

s <- colnames(ind.spp)
tm_shape(basins) + 
  tm_polygons(col = paste0(s[5]))

# plot predicted indicator species occurrence and scenarios

