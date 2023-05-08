# just find basins where a species can detect different threats, and plot that species

library(sf)
library(tmap)
library(tidyverse)
library(raster)
library(abind)
library(scales)
library(Hmsc)

# load data

dat <- read.csv('data/waterbird-basin-threats-100_bioregion.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
qld <- st_read('data/qld-shp/queensland-polygon.shp')
m <- readRDS('outputs/models/mod-spatialRF_final.rds')
coefs <- read.csv('outputs/spp-beta-coefs.csv')

# Australasian bittern responds positiveley to consumptive waterloss, nitrogen, and phosphorus

XDataNew <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Phosphorus_Loading = max(m$XData$Phosphorus_Loading))

XDataNew2 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Consumptive_Water_Loss = max(m$XData$Consumptive_Water_Loss))

XDataNew3 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Nitrogen_Loading = max(m$XData$Nitrogen_Loading))

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

# make predictions of spp. prob of occurrence where nitrogen loading is at max value, given known environmental gradient
predY4 <- predict(m, XData=XDataNew3, studyDesign=m$studyDesign,
                  ranLevels=m$ranLevels, expected=TRUE)

# extract the mean prob. of occurrence from our posterior density distribution of predictions 
EpredY <- data.frame(apply(abind(predY, along = 3), c(1,2), mean))
EpredY2 <- data.frame(apply(abind(predY2, along = 3), c(1,2), mean))
EpredY3 <- data.frame(apply(abind(predY3, along = 3), c(1,2), mean))
EpredY4 <- data.frame(apply(abind(predY4, along = 3), c(1,2), mean))

# add in Hydrobasin ID (for plotting)
EpredY$HYBAS_ID <- XDataNew[,1]

# calculating the predicted change in prob. of occurrence in a species given an increase threat
# subtract the baseline probability of occurrence from our scenario-based prob. occurrence where we increase threat levels
EpredY$Probability <- rescale(EpredY2$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1)) # rescale from 0-1 for easier visualisation
EpredY$Probability2 <- rescale(EpredY3$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1))
EpredY$Probability3 <- rescale(EpredY4$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1))

# classify basins so that if has detection probability 50% above baseline, can detect it

plotdf <- EpredY %>% 
  dplyr::select(-c(Australasian.Darter:Musk.Duck)) %>% 
  mutate(Phosphorus_Loading = ifelse(Probability > 0.5, 1, 0)) %>% 
  mutate(Consumptive_Water_Loss = ifelse(Probability2 > 0.5, 1, 0)) %>% 
  mutate(Nitrogen_Loading = ifelse(Probability > 0.5, 1, 0)) %>% 
  mutate(total = Phosphorus_Loading + Consumptive_Water_Loss + Nitrogen_Loading) %>% 
  filter(total == 1) %>% 
  pivot_longer(cols = Phosphorus_Loading:Nitrogen_Loading, names_to = 'Threat', values_to = 'Detect')

# map

predY.sf <- basins.qld %>% 
  inner_join(plotdf, by = 'HYBAS_ID')

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(predY.sf) +
  tm_polygons('Threat', palette = 'Set2', legend.show = T) + # turn on legend by saying T
  tm_layout(legend.position = c(0.5, 0.85),
            legend.title.size = 0.8,
            legend.text.size = 0.45)
m3
tmap_save(m3, 'outputs/map-scenario_Australasian.Bittern.png', width = 6, height = 6)

### End here
