# plot covariate effect, species interactions, etc.

# libaries

library(Hmsc)
library(MASS)
library(tidyverse)
library(RColorBrewer)
library(randomcoloR)
library(ggcorrplot)

dat <- read.csv('data/waterbird-basin-threats_final.csv')
y <- as.matrix(dat %>% dplyr::select(Australasian.Darter:Musk.Duck))
x <- dat %>% dplyr::select(Aquaculture_Pressure:pre_mm_syr)

# load model

m <- readRDS('outputs/models/mod-non-spatial.rds')
m

# partition explained variance among fixed and random effects

head(m$X) # check to make sure groupnames match column names here 
VP <- computeVariancePartitioning(m, group = c(rep(1,12), rep(2,3)), groupnames = c("river-threat", "env-covariate"))
plotVariancePartitioning(m, VP = VP)

vardf <- data.frame(VP$vals)
vardf$cat <- c('River threat', 'Environmental', 'Species interactions')
vardf.long <- vardf %>%
  pivot_longer(Australasian.Darter:Musk.Duck, names_to = 'Species', values_to = 'Variance_prop')

ggplot(vardf.long) +
  aes(x = Species, y = Variance_prop, fill = cat) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  ylab('Proportion variance explained') +
  xlab('') +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(3, 'Set2')) +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')

ggsave('outputs/spp-var-explained.png', width = 7, height = 6)

# plot estimats of Beta parameters for each spp as heatmap
# red = posterior probability for the parameter being positive is greater than 0.95 
# blue = posterior probability for the parameter being negative is greater than 0.95 
# white = no strong statistical support for parameter being positive or negative

postBeta <- getPostEstimate(m, parName = "Beta")

pdf('outputs/spp-regression-params.pdf')
plotBeta(m, post = postBeta, param = "Support", 
         supportLevel = 0.95, mgp = c(1.5,1.5,0), 
         mar = c(5,0,5,0), cex = c(0.4,0.4,0.5))
dev.off()

beta.mean <- data.frame(postBeta$mean)
beta.mean$variable <- c(as.character(colnames(m$X)))
beta.mean.long <- beta.mean %>%
  pivot_longer(Australasian.Darter:Musk.Duck, names_to = 'Species', values_to = 'Beta_mean') %>% 
  filter(variable != '(Intercept)')

beta.support <- data.frame(postBeta$support)
beta.support$variable <- c(as.character(colnames(m$X)))
beta.support.long <- beta.support %>%
  pivot_longer(Australasian.Darter:Musk.Duck, names_to = 'Species', values_to = 'Beta_support') %>% 
  filter(variable != '(Intercept)')

beta.supportNeg <- data.frame(postBeta$supportNeg)
beta.supportNeg$variable <- c(as.character(colnames(m$X)))
beta.supportNeg.long <- beta.supportNeg %>%
  pivot_longer(Australasian.Darter:Musk.Duck, names_to = 'Species', values_to = 'Beta_supportNeg') %>% 
  filter(variable != '(Intercept)')

beta.mean.long$support <- beta.support.long$Beta_support
beta.mean.long$supportNeg <- beta.supportNeg.long$Beta_supportNeg

# filter for spp-ind relationships with greater than 95% support

beta.pos <- filter(beta.mean.long, support > 0.95)
beta.pos$direction <- rep('positive', nrow(beta.pos))
beta.neg <- filter(beta.mean.long, supportNeg > 0.95)
beta.neg$direction <- rep('negative', nrow(beta.neg))

beta.posneg <- rbind(beta.pos, beta.neg)

ggplot(filter(beta.posneg, !variable %in% c("tmp_dc_syr", "pre_mm_syr", "ele_mt_sav"))) +
  aes(x = Species, y = Beta_mean, col = direction, shape = variable) +
  geom_point() +
  scale_shape_manual(values = c(7, 8,9, 10, 11, 12, 14, 15, 16))+
  theme_classic() +
  ylab('Coefficient') +
  xlab('') +
  geom_vline(xintercept = 0) +
  coord_flip() +
  scale_fill_manual(values = randomColor(16)) +
  theme(legend.title = element_blank())

ggsave('outputs/spp-beta-coefs-95-support.png', width = 8, height = 6)

# summarise spp coefficents from posterior distribution 
# mean, median, and credible intervals

mpost <- convertToCodaObject(m)

beta.coefs1 <- data.frame(summary(mpost$Beta[[1]])$statistics) # first chain
beta.quant1 <- data.frame(summary(mpost$Beta[[1]])$quantiles) # first chain
variable <- rep(as.character(colnames(m$X)), nrow(beta.coefs1)/(length(colnames(m$X))))
names <- colnames(beta.mean[,-length(beta.mean)])
spp <- c()
for(i in 1:length(names)){
  sp <- rep(names[i], nrow(beta.coefs1)/(length(names)))
  spp <- c(spp, sp)
}

beta.coefs1 <- beta.coefs1 %>% 
  mutate(CI.2.5 = beta.quant1$X2.5.,
         CI.97.5 = beta.quant1$X97.5.,
         variable = variable,
         species = spp) %>% 
  filter(variable != 'intercept')

for(i in 2:length(colnames(m$X))){
  pl <- filter(beta.coefs1, variable == paste(colnames(m$X)[i]))
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
  
  ggsave(paste0('outputs/', colnames(m$X)[i], '.png'), width = 8, height = 6)
}

# extract estimates of species-species associations
# computeAssociations function converts the covariances to correlations

OmegaCor <- computeAssociations(m)
corr.mat <- OmegaCor[[1]]$mean
ggcorrplot(corr.mat, tl.cex = 7)
ggsave('outputs/spp-correlation-matrix.png', width = 10, height = 5.5)

# extract species and site loadings on LVs

etaPost <- getPostEstimate(m, "Eta") # site loadings 
lambdaPost <- getPostEstimate(m, "Lambda") # spp loadings

# plot variation in spp richness over enviro gradients
# Note: no non-focal variables are listed, so  predictions for the net effect of
# the focal variable, that is the total effect of the focal variable and 
# any non-focal variables covarying with the focal variable
##TODO: see finish bird vignette for making predictions under concrete scenarios

for(i in 2:length(colnames(m$X))){
  Gradient <- constructGradient(m, focalVariable = paste(colnames(m$X)[i]))
  Gradient$XDataNew
  predY <- predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                   ranLevels=Gradient$rLNew, expected=TRUE)
  pdf(paste0('outputs/spp-richness_', colnames(m$X)[i], '.pdf'), width =9, height =7) 
  plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE)
  dev.off()
}

##TODO adapt from Renee's study if want to make predictions of 
# prob. of spp occurrence under different threat levels


# save cpue predictions under different scenarios for plotting later
# scenario 1: increase mean patch size of mangroves by 20 ha in each basin

Xdata <- m.spatial$XData
Xdata[,'mang_area_mn'] <- Xdata['mang_area_mn'] + 20 # double check units on this stuff
predY <- predict(m.spatial, XData=Xdata, studyDesign=m.spatial$studyDesign,
                 ranLevels=m.spatial$rL, expected=TRUE)
predY.array <- simplify2array(predY)
predY.m <- apply(predY.array, 1:2, mean) # double-check if should use median instead mean, might do some sensitivity for final analysis

# plot counterfactural vs restoration scenario
# calc total cpue in each basin

cpue.bas <- data.frame(BASIN = x$BASIN_NAME, exp(m.spatial$Y)) %>% 
  group_by(BASIN) %>% 
  summarise_all(sum) %>% 
  mutate(scenario = rep('BAU', nrow(.)))

cpue.bas.scenario <- data.frame(BASIN = x$BASIN_NAME, exp(predY.m)) %>% 
  group_by(BASIN) %>% 
  summarise_all(sum) %>% 
  mutate(scenario = rep('mang20', nrow(.)))

df <- rbind(cpue.bas, cpue.bas.scenario)

write.csv(df, 'Outputs/hmsc-scenario-pred.csv', row.names = F)


