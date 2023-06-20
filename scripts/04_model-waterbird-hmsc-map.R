# map threats and species' distributions

library(sf)
library(tmap)
library(tidyverse)
library(raster)
library(abind)
library(scales)
library(Hmsc)
library(plotly)

# load data

dat <- read.csv('data/waterbird-basin-threats-100_bioregion_notshore.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
qld <- st_read('data/qld-shp/queensland-polygon.shp')
m <- readRDS('outputs/models/mod-spatialRF_final_notshore.rds')
m

# use model to predict probability of spp occurrence will increase or decrease depending on threat
# Black swan and increased phosphorus and consumptive waterloss
# set predictor of interest to max of observed environmental gradient at our basins

XDataNew <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Phosphorus_Loading = max(m$XData$Phosphorus_Loading))

XDataNew2 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Consumptive_Water_Loss = max(m$XData$Consumptive_Water_Loss))

# make predictions

# first make baseline predictions of prob. of occurrence for all species
# use the same covriate/predictor values as we built the model on at all of the hydrobasins
predY <- predict(m, XData=m$XData, studyDesign=m$studyDesign,
                 ranLevels=m$ranLevels, expected=TRUE)

# make predictions of spp. prob of occurrence where phosphorus is at max value, given known environmental gradient
predY2 <- predict(m, XData=XDataNew, studyDesign=m$studyDesign,
                 ranLevels=m$ranLevels, expected=TRUE)

# make predictions of spp. prob of occurrence where consumptive water loss is at max value, given known environmental gradient
predY3 <- predict(m, XData=XDataNew2, studyDesign=m$studyDesign,
                  ranLevels=m$ranLevels, expected=TRUE)

# extract the mean prob. of occurrence from our posterior density distribution of predictions 
EpredY <- data.frame(apply(abind(predY, along = 3), c(1,2), mean))
EpredY2 <- data.frame(apply(abind(predY2, along = 3), c(1,2), mean))
EpredY3 <- data.frame(apply(abind(predY3, along = 3), c(1,2), mean))

# add in Hydrobasin ID (for plotting)
EpredY$HYBAS_ID <- XDataNew[,1]

# calculating the predicted change in prob. of occurrence in a species given an increase threat
# subtract the baseline probability of occurrence from our scenario-based prob. occurrence where we increase threat levels
EpredY$Probability <- rescale(EpredY2$Black.Swan - EpredY$Black.Swan, to = c(0,1)) # rescale from 0-1 for easier visualisation
EpredY$Probability2 <- rescale(EpredY3$Black.Swan - EpredY$Black.Swan, to = c(0,1))

# group the changes in prob. of occurrence into categories based on whether increased occurrence
# is the same or different for different threats

plotdf <- EpredY %>% 
  mutate(group = ifelse(Probability <= quantile(EpredY$Probability, 0.5) & Probability2 <= quantile(EpredY$Probability2, 0.5), 1, NA)) %>% 
  mutate(group = ifelse(Probability >= quantile(EpredY$Probability, 0.5) & Probability2 <= quantile(EpredY$Probability2, 0.5), 2, group)) %>% 
  mutate(group = ifelse(Probability <= quantile(EpredY$Probability, 0.5) & Probability2 >= quantile(EpredY$Probability2, 0.5), 3, group)) %>% 
  mutate(group = ifelse(Probability >= quantile(EpredY$Probability, 0.5) & Probability2 >= quantile(EpredY$Probability2, 0.5), 4, group)) %>% 
  mutate(group  = factor(group))

# plot correlations between changes in prob. of occurrence given two threats

a <- ggplot(select(plotdf, group, HYBAS_ID, Probability, Probability2)) +
  geom_point(aes(x = Probability, 
                 y = Probability2,
                 label = HYBAS_ID,
                 colour = group)) +
  xlab('Increased probability of occurrence with high Phosphorus') +
  ylab('Increased probability of occurrence with high Consumptive Water Loss') +
  geom_vline(xintercept = quantile(EpredY$Probability, 0.5), alpha = 0.5) +
  geom_hline(yintercept = quantile(EpredY$Probability2, 0.5), alpha = 0.5) +
  scale_color_brewer(palette = 'Set2') +
  #geom_abline() +
  theme_classic() +
  theme(legend.position = 'none')


ggsave('outputs/black-swan-scenario-correlation_notshore.png', width = 5.5, height = 5.5)

# map predictions

predY.sf <- basins.qld %>% 
  inner_join(plotdf, by = 'HYBAS_ID')

m1 <- tm_shape(qld) +
  tm_fill() + # can change this to tm_polygons if you want borders
tm_shape(predY.sf) +
  tm_polygons('Probability')

m2 <- tm_shape(qld) +
  tm_fill() +
tm_shape(predY.sf) +
  tm_polygons('Probability2')

mm <- tmap_arrange(m1, m2)
mm

tmap_save(mm, 'outputs/map-scenario_black-swan_notshore.png', width = 6.666, height = 4)

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(predY.sf) +
  tm_fill('group', palette = 'Set2', legend.show = F, title = 'Indicator group') +# turn on legend by saying T
  tm_layout(frame = T,
            legend.position = c(0.4,0.8),
            legend.width = 2) +
  tm_add_legend(type = "fill", 
                col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
                labels = c("No response", "Indicates Phosphorous loading (PL)", "Indicates Consumptive Water Loss (CWL)", "Indicates CWL & PL"))
m3
tmap_save(m3, 'outputs/map-scenario_black-swan-categories_notshore.png',  height = 6, width = 7)

### End here

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(filter(predY.sf, HYBAS_ID == 5080068310)) +
  tm_polygons('group', palette = 'Set2', legend.show = T, title = 'Indicator group') +# turn on legend by saying T
  tm_layout(frame = F,
            legend.position = c(0.5,0.8))
m3
tmap_save(m3, 'outputs/map-scenario_black-swan-categories_notshore_1.png',  height = 4, width = 3)

