# plot covariate effect, species interactions, etc.

# libaries

library(Hmsc)
library(MASS)
library(tidyverse)
library(RColorBrewer)
library(randomcoloR)
library(ggcorrplot)
library(factoextra)
library(FactoMineR)
library(cluster)
library(scales)
library(ggrepel)
library(abind)
source('scripts/plot-helpers.R')

# load model

m <- readRDS('outputs/models/mod-spatialRF_final.rds')

# partition explained variance among fixed and random effects

VP <- computeVariancePartitioning(m)
plotVariancePartitioning(m, VP = VP)

vardf <- data.frame(VP$vals)
vardf$cat <- c(rep('Threats', 5), rep('Environmental', 1), 'Species interactions', 'Bioregion')
vardf.long <- vardf %>%
  pivot_longer(Australasian.Darter:Musk.Duck, names_to = 'Species', values_to = 'Variance_prop')

ggplot(vardf.long) +
  aes(x = Species, y = Variance_prop, fill = cat) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  ylab('Proportion variance explained') +
  xlab('') +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(4, 'Set2')) +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')

ggsave('outputs/spp-var-explained.png', width = 7, height = 10)

# community-level effects - species richness vs. threats
# here we set non.focal variables to 1, so just get the marginal effects of each focal variable
# this is the effect of the focal variable, independent of the others (fixes non-focal variables to most likely values)

mpost <- convertToCodaObject(m)

prednames <- colnames(m$XData)
spp <- list()
trend <- list()
for(i in seq_along(prednames)){
  Gradient <- constructGradient(m, focalVariable= prednames[i],
                               non.focalVariables=1)
  predY <- predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                  ranLevels=Gradient$rLNew, expected=TRUE)
  plotdat <- get_spprich(m, Gradient, pred=predY, measure="S", las=1,
                         showData = TRUE, showPosteriorSupport = FALSE)
  spp[[i]] <- data.frame(Threat = prednames[i], plotdat[[2]])
  trend[[i]] <- data.frame(Threat = prednames[i], plotdat[[1]])
}
sppdat <- do.call(rbind, spp) %>% filter(Threat != 'pre_mm_syr')
trendat <- do.call(rbind, trend) %>% filter(Threat != 'pre_mm_syr')

ggplot() +
  geom_point(data = sppdat, aes(x = x, y = y), col = 'lightgrey') +
  geom_line(data = trendat, aes(x = x, y = richness_50), size = 1) +
  geom_ribbon(data = trendat, aes(x = x, ymin = richness_2.5, ymax = richness_97.5), fill = "grey", alpha = 0.3) +
  geom_ribbon(data = trendat, aes(x = x, ymin = richness_25, ymax = richness_75), fill = "darkgrey", alpha = 0.5) +
  facet_wrap(~Threat, scales = 'free_x') +
  xlab('') +
  ylab('Species Richness') +
  theme_classic()

ggsave('outputs/species-trends.png', width = 6, height = 4)

# plot individual species marginal effects 

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support",
         plotTree = FALSE, supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))

# summarise beta coefficents from posterior distribution for more enhanced control over plotting etc.
# extract mean coefs and credible intervals from each mcmc chain

beta.coefs <- list()
for(i in seq_along(mpost$Beta)){ # loop through mcmc chains to extract
  b.stats <- data.frame(summary(mpost$Beta[[i]])$statistics)
  b.quant <- data.frame(summary(mpost$Beta[[i]])$quantiles)
  beta.coefs[[i]] <- b.stats %>%
    mutate(CI.2.5 = b.quant$X2.5.,
           CI.97.5 = b.quant$X97.5.,
           CI.25 = b.quant$X25.,
           CI.75 = b.quant$X75.,
           variable = rep(as.character(colnames(m$X)), length(colnames(m$Y))),
           species = rep(as.character(colnames(m$Y)), each = length(colnames(m$X)))) %>%
    filter(variable != '(Intercept)') %>%
    pivot_longer(Mean:CI.75, names_to = 'stat', values_to = 'val')
}

# bind and calculate mean aross mcmc chains

beta <- do.call(rbind, beta.coefs) %>%
  group_by(variable, species, stat) %>%
  summarise(val = median(val)) %>%
  pivot_wider(names_from = 'stat', values_from = 'val') %>% 
  mutate(probable_pos = ifelse(CI.75>0 & CI.25>0, 1, 0),
         probable_neg = ifelse(CI.75<0 & CI.25<0, 1, 0))

# plot beta coefs with either weak effect (50% CIs) or strong effect (95% CIS)

beta2 <- beta %>% 
  filter(probable_pos == 1 | probable_neg == 1) %>% 
  mutate(probable_pos = ifelse(probable_pos == 1, 'Positive', 'Negative')) %>% 
  filter(variable != 'pre_mm_syr')

ggplot(beta2) +
  aes(x = species, y = Mean, col = variable) +
  geom_point() +
  geom_linerange(aes(ymin=CI.2.5, ymax=CI.97.5), alpha = 0.5) + # 95% CiS
  geom_linerange(aes(ymin=CI.25, ymax=CI.75), size = 2, alpha = 0.5) + # 50% CIs
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = distinctColorPalette(16)) +
  facet_wrap(~variable, ncol = 5) +
  theme_classic() +
  ylab('Coefficient') +
  xlab('') +
  coord_flip() +
  theme(legend.position = 'none') 

ggsave('outputs/beta-coefficients_100.png', width = 10, height = 8.5)

# turn into heatmap for easier visualisation - +ve or -ve responses

ggplot(beta2) +
  aes(x = variable, y = species, fill = factor(probable_pos)) +
  geom_tile() +
  xlab('') +
  ylab('') +
  #facet_wrap(~, ncol=5) +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1))

ggsave('outputs/indicator-response-heatmap.png', width = 5, height = 10)

# cluster analysis to identify species that are respoding simililarly to threats

spp.response <- beta %>% 
  filter(variable != 'pre_mm_syr') %>% 
  mutate(response = probable_pos) %>% 
  mutate(response = ifelse(probable_neg == 1, 2, response)) %>% 
  pivot_wider(id_cols = species, names_from = variable, values_from = response)

set.seed(123)

dat2 <- data.frame(spp.response[,-1])
row.names(dat2) <- spp.response$species

# try kmeans

fviz_nbclust(dat2, kmeans, method = "wss")
km <- kmeans(dat2, 7, nstart = 50)

# visualise

fviz_cluster(km, data = dat2,
             axes = c(1,3),
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07", "black", 'pink', 'purple', 'gold'),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot"
)

# interpretation plots

dat.clust <- data.frame(spp.response, cluster = km$cluster)
write.csv(dat.clust, 'data/cluster-dat.csv', row.names = F)

# heatmap of scores for each indicator in each cluster

dat_sum <- dat.clust %>% 
  pivot_longer(cols = Aquaculture_Pressure:Phosphorus_Loading,
               names_to = 'indicator',
               values_to = 'response') %>% 
  mutate(positive = ifelse(response == 1, 1, 0),
         negative = ifelse(response == 2, 1, 0),
         neutral = ifelse(response == 0, 1, 0)) %>% 
  group_by(cluster, indicator) %>% 
  summarise(Positive = sum(positive)/n(),
            Negative = sum(negative)/n(),
            Neutral = sum(negative)/n()) %>% 
  pivot_longer(cols = Positive:Negative,
    values_to = 'Proportion',
               names_to = 'direction')

ggplot(dat_sum) +
  aes(x = factor(direction), y = indicator, fill = Proportion) +
  geom_tile() +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1) +
  xlab('') +
  ylab('') +
  facet_wrap(~ cluster, ncol=4) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.2))

ggsave('outputs/cluster-heatmap.png', width = 6.5, height = 3.5)

beta4 <- beta %>% 
  left_join(dplyr::select(dat.clust, species, cluster)) %>% 
  filter(probable_pos == 1 | probable_neg == 1) %>% 
  filter(variable != 'pre_mm_syr')

ggplot(beta4) +
  aes(x = species, y = Mean, col = factor(cluster)) +
  geom_point() +
  geom_linerange(aes(ymin=CI.2.5, ymax=CI.97.5), alpha = 0.5) + # 95% CiS
  geom_linerange(aes(ymin=CI.25, ymax=CI.75), size = 2, alpha = 0.5) + # 50% CIs
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = distinctColorPalette(16)) +
  facet_wrap(~variable, ncol = 5) +
  theme_classic() +
  ylab('Coefficient') +
  xlab('') +
  coord_flip() #+
 # theme(legend.position = 'none') 

ggsave('outputs/beta-coefficients_100_cluster.png', width = 10, height = 8.5)

# now plot spp responses grouped by their cluster

beta3 <- beta2 %>% 
  left_join(dat.clust, by = 'species')

ggplot(beta3) +
  aes(x = variable, y = species, fill = factor(probable_pos)) +
  geom_tile() +
  xlab('') +
  ylab('') +
  facet_wrap(~cluster, ncol=4, scales = 'free_y') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.06),
        axis.text.x = element_text(angle = 45, vjust = 0.98, hjust = 1))

ggsave('outputs/indicator-response-heatmap.png', width = 12, height = 7)

### extra for looking at residual correlations between species - species interactions or missing covariates

# extract estimates of species-species associations
# computeAssociations function converts the covariances to correlations

OmegaCor <- computeAssociations(m)
corr.mat <- OmegaCor[[1]]$mean # first random level
ggcorrplot(corr.mat, tl.cex = 7)
ggsave('outputs/spp-correlation-matrix.png', width = 10, height = 5.5)

# extract species and site loadings on LVs

etaPost <- getPostEstimate(m, "Eta") # site loadings 
lambdaPost <- getPostEstimate(m, "Lambda") # spp loadings

# predict species richness (S) and probability of occurrence (Y)
# marginal effects (prediction type 1)

# set all non-focal variables to prediction type '1' (ie. marginal effects)

spp <- which(colnames(m$Y) == 'Pacific.Golden.Plover')
focalv <- 'Aquaculture_Pressure'
#focalv <- 'Nitrogen_Loading'
nonfocal <- colnames(dplyr::select(m$XData, -focalv))

ls <- list()
for(i in 1:length(nonfocal)){
  xx <- list(1) # set type to 1 for marginal effects of focal variable
  names(xx) <- nonfocal[i]
  ls[[i]] <- xx
}
names(ls) <- nonfocal

Gradient <- constructGradient(m, focalVariable = focalv, non.focalVariables = ls)

predY <- predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                 ranLevels=Gradient$rLNew, expected=TRUE)

plotGradient(m, Gradient, pred=predY, measure="Y", index = spp, showData = TRUE) # single spp predictions
#plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE) # spp richness predictions

# set nitrogent to minimum and get effect of aquaculture when nitrogen is low

ls[[which(names(ls) == 'Nitrogen_Loading')]] <- list(3, max(m$XData$Nitrogen_Loading))

Gradient <- constructGradient(m, focalVariable = focalv, non.focalVariables = ls)

predY <- predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                 ranLevels=Gradient$rLNew, expected=TRUE)

plotGradient(m, Gradient, pred=predY, measure="Y", index = spp, showData = TRUE) # single spp predictions

# can I just  loop through the above, and find out when the indicator capacity is conditional another threat 
# (e.g., the effect of a threat changes if considering the min or the max of the gradient)

# for each species, loop through each 