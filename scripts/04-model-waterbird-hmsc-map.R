# map threats and species' distributions

library(sf)
library(tmap)
library(tidyverse)
library(raster)
library(abind)
library(scales)

# load data

dat <- read.csv('data/waterbird-basin-threats-100_bioregion.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp') %>% 
  st_transform(crs = 4326)
qld <- st_read('data/qld-shp/queensland-polygon.shp') %>% 
  st_transform(crs = 4326)
m <- readRDS('outputs/models/mod-spatialRF_final.rds')
m

# use model to predict probability of spp occurrence will increase or decrease depending on threat
# Black swan and increased phosphorus

alph <- 0.5 # percent increase in threat level

XDataNew <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Phosphorus_Loading = 1)
 # mutate(Phosphorus_Loading = Phosphorus_Loading + (Phosphorus_Loading*alph)) %>% 
  #mutate(Phosphorus_Loading = ifelse(Phosphorus_Loading == 0, 
   #                                  min(m$XData$Phosphorus_Loading[m$XData$Phosphorus_Loading>0]),
    #                                 Phosphorus_Loading))
XDataNew2 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Consumptive_Water_Loss = 1)
  #mutate(Consumptive_Water_Loss = Consumptive_Water_Loss + Consumptive_Water_Loss*alph) %>% 
  #mutate(Consumptive_Water_Loss = ifelse(Consumptive_Water_Loss == 0, 
   #                                  min(m$XData$Consumptive_Water_Loss[m$XData$Consumptive_Water_Loss>0]),
    #                                 Consumptive_Water_Loss))

# make predictions

predY <- predict(m, XData=m$XData, studyDesign=m$studyDesign,
                 ranLevels=m$ranLevels, expected=TRUE)
predY2 <- predict(m, XData=XDataNew, studyDesign=m$studyDesign,
                 ranLevels=m$ranLevels, expected=TRUE)
predY3 <- predict(m, XData=XDataNew2, studyDesign=m$studyDesign,
                  ranLevels=m$ranLevels, expected=TRUE)
EpredY <- data.frame(apply(abind(predY, along = 3), c(1,2), mean))
EpredY2 <- data.frame(apply(abind(predY2, along = 3), c(1,2), mean))
EpredY3 <- data.frame(apply(abind(predY3, along = 3), c(1,2), mean))
EpredY$HYBAS_ID <- XDataNew[,1]
EpredY$Probability <- rescale(EpredY2$Black.Swan - EpredY$Black.Swan, to = c(0,1))
EpredY$Probability2 <- rescale(EpredY3$Black.Swan - EpredY$Black.Swan, to = c(0,1))

# compare predictions

EpredY <- EpredY %>% 
  mutate(group = ifelse(Probability < 0.5 & Probability2 < 0.5, 1, NA)) %>% 
  mutate(group = ifelse(Probability > 0.5 & Probability2 < 0.5, 2, group)) %>% 
  mutate(group = ifelse(Probability < 0.5 & Probability2 > 0.5, 3, group)) %>% 
  mutate(group = ifelse(Probability > 0.5 & Probability2 > 0.5, 4, group)) %>% 
  mutate(group  = factor(group))

ggplot(EpredY) +
  geom_point(aes(x = Probability, 
                 y = Probability2,
                 colour = group)) +
  xlab('Probability of detection with increased Phosphorus') +
  ylab('Probability of detection with increased Consumptive Water Loss') +
  geom_vline(xintercept = 0.5) +
  geom_hline(yintercept = 0.5) +
  scale_color_brewer(palette = 'Set2') +
  #geom_abline() +
  theme_classic() +
  theme(legend.position = 'none')

ggsave('outputs/black-swan-scenario-correlation.png', width = 5.5, height = 5)

# map predictions

predY.sf <- basins.qld %>% 
  inner_join(EpredY)

m1 <- tm_shape(qld) +
  tm_fill() +
tm_shape(predY.sf) +
  tm_polygons('Probability')

m2 <- tm_shape(qld) +
  tm_fill() +
tm_shape(predY.sf) +
  tm_polygons('Probability2')

mm <- tmap_arrange(m1, m2)
mm

tmap_save(mm, 'outputs/map-scenario_black-swan.png', width = 6.666, height = 4)

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(predY.sf) +
  tm_polygons('group', palette = 'Set2')
m3
tmap_save(m3, 'outputs/map-scenario_black-swan-categories.png', width = 6, height = 6)

### End here


set.seed(123)
EpredY$cluster <- as.factor(kmeans(EpredY, 7)$cluster) # regions of common profile
EpredY.scen$cluster <- as.factor(kmeans(EpredY.scen, 7)$cluster) # regions of common profile

EpredY$cumul.threat <- rowSums(dplyr::select(XDataNew, Aquaculture_Pressure:Wetland_Disconnectivity))
EpredY.scen$cumul.threat <- rowSums(dplyr::select(Xdata.scen, Aquaculture_Pressure:Wetland_Disconnectivity))

EpredY <- cbind(EpredY, XDataNew[,-1])
EpredY.scen <- cbind(EpredY.scen, Xdata.scen[,-1])


# maps of threat, spp richness and community composition

m1 <- tm_shape(predY.sf) +
  tm_polygons(col = 'cumul.threat', title = 'Threat')
m1

m2 <- tm_shape(predY.sf) +
  tm_polygons('spp.rich', title = 'Species richness')
m2

m3 <- tm_shape(predY.sf) +
  tm_polygons('cluster', title = 'Common profile')
m3

mm <- tmap_arrange(m1, m2, m3)
mm

tmap_save(mm, 'outputs/map-thrt-richness-profile.png', width = 10, height = 4)

# maps of change in prob. of occurrence of focal spp and change in threat level
# under scenario

m4 <- tm_shape(predY.scen.sf) +
  tm_polygons('change.threat', title = 'Pesticide Loading')
m4

m5 <- tm_shape(predY.scen.sf) +
  tm_polygons('change.spp', title = 'Pied Cormorant')
m5

mm <- tmap_arrange(m4, m5)
mm

tmap_save(mm, 'outputs/map-scenario_pesticide-pied-cormorant.png', width = 6.666, height = 4)
