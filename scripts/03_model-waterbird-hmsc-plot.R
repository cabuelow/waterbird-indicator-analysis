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

for(i in seq_along(prednames)){
  
  Gradient = constructGradient(m, focalVariable= prednames[i],
                               non.focalVariables=1)
  
  predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew,
                  ranLevels=Gradient$rLNew, expected=TRUE)
  
  png(paste0('outputs/community_', prednames[i], '.png'), height = 300, width = 400)
  plotGradient(m, Gradient, pred=predY, measure="S", las=1,
               showData = TRUE, showPosteriorSupport = FALSE)
  dev.off()
}

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


# TODO; spatially explicity predicitons


### extra dont use below!

# cluster analysis to identify species that are respoding simililarly to threats

spp.response <- beta %>% 
  filter(variable != 'pre_mm_syr') %>% 
  mutate(response = probable_pos) %>% 
  mutate(response = ifelse(probable_neg == 1, 2, response)) %>% 
  pivot_wider(id_cols = species, names_from = variable, values_from = response)

set.seed(123)

dat2 <- data.frame(spp.response[,-1])
row.names(dat2) <- spp.response$species

# visualise distance values for clustering

fviz_dist(get_dist(dat2, method = 'euclidian'), lab_size = 8)

# hierarchical cluster analysis

dat_hclust <- hcut(get_dist(dat2, method = 'euclidian'), k = 7, isdiss = T)

# plot dendrogram

fviz_dend(dat_hclust, rect = TRUE)
#ggsave('outputs/typologies/governance-typology-dendrogram.png', width = 11, height = 7)

# visualize the silhouette

fviz_silhouette(dat_hclust)

# try kmeans

km <- kmeans(dat2, 8)

# visualise

fviz_cluster(km, data = dat2,
             axes = c(1,3),
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07", "black", 'pink', 'purple'),
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
  facet_wrap(~ cluster, ncol=5) +
  theme_classic()

ggsave('outputs/cluster-heatmap.png', width = 8.3, height = 2)

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

# pca

set.seed(123)
res.pca <- PCA(dat2, graph = FALSE)

# visualize eigenvalues/variances

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

# extract the results for variables

var <- get_pca_var(res.pca)
var

# contributions of variables to PCs

fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 4, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 5, top = 10)

# plot

fviz_pca_var(res.pca, axes = c(1,2), col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 

fviz_pca_var(res.pca, axes = c(1,3), col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 

fviz_pca_var(res.pca, axes = c(1,4), col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 

# save scores for countries on each pca

scores <- data.frame(get_pca_ind(res.pca)$coord)
colnames(scores) <- c('Nitro_conswater', 'Aqua_phos', 'Pest_phos',
                      'Techniques', 'Services')
scores$species <- spp.response$species
#scores$Institution <- dat$Institution
scores <- scores %>% 
  left_join(dat.clust)

# plot

ggplot(scores) +
  aes(x =  Nitro_conswater, y = Aqua_phos, 
      #size = Communities_eia,
      col = factor(cluster), 
      #label = Country.y.y, 
      alpha = 0.9) +
  geom_point(size = 2) +
  geom_vline(aes(xintercept = 0),linetype = 'dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', alpha = 0.5) +
  #geom_text_repel(size = 2, show.legend = F, max.overlaps = 20) +
  scale_color_brewer(palette = 'Set2') +
  theme_classic() +
  guides(alpha = FALSE,
         color=guide_legend(title="Typology")) #+
#geom_text(show.legend = FALSE)

ggsave('outputs/pca_biodiverse-techniques-services-biodiversity.png', width = 5, height = 3)

###

ggplot(scores) +
  aes(x =  Techniques.services, y =  Enviro.impacts.socio.cultural, 
      #size = Communities_eia,
      col = factor(cluster), 
      #label = Country.y.y, 
      alpha = 0.9) +
  geom_point(size = 2) +
  geom_vline(aes(xintercept = 0),linetype = 'dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', alpha = 0.5) +
  #geom_text_repel(size = 2, show.legend = F, max.overlaps = 20) +
  scale_color_brewer(palette = 'Set2') +
  theme_classic() +
  guides(alpha = FALSE,
         color=guide_legend(title="Typology")) #+
#geom_text(show.legend = FALSE)

ggsave('outputs/pca_techniques-services-enviro-impacts-socio-cultural.png', width = 5, height = 3)



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