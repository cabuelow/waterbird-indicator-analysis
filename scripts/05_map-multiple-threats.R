# just find basins where a species can detect different threats, and plot that species

library(sf)
library(tmap)
library(tidyverse)
library(raster)
library(abind)
library(scales)
library(Hmsc)
# NOTE** Below is all coded terribly inefficiently, need to fix when have time...

# load data

dat <- read.csv('data/waterbird-basin-threats-100_bioregion_notshore.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
qld <- st_read('data/qld-shp/queensland-polygon.shp')
m <- readRDS('outputs/models/mod-spatialRF_final_notshore.rds')
coefs <- read.csv('outputs/spp-beta-coefs_notshore.csv')

# Australasian bittern responds positiveley to consumptive waterloss, nitrogen, and phosphorus

XDataNew <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Phosphorus_Loading = max(m$XData$Phosphorus_Loading))

XDataNew2 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Consumptive_Water_Loss = max(m$XData$Consumptive_Water_Loss))

XDataNew3 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Nitrogen_Loading = max(m$XData$Nitrogen_Loading))

XDataNew4 <- data.frame(HYBAS_ID = dat$HYBAS_ID, m$XData) %>% 
  mutate(Nitrogen_Loading = max(m$XData$Pesticide_Loading))

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

# make predictions of spp. prob of occurrence where pesitcide loading is at max value, given known environmental gradient
predY5 <- predict(m, XData=XDataNew4, studyDesign=m$studyDesign,
                  ranLevels=m$ranLevels, expected=TRUE)

# extract the mean prob. of occurrence from our posterior density distribution of predictions 
EpredY <- data.frame(apply(abind(predY, along = 3), c(1,2), mean))#baseline
EpredY2 <- data.frame(apply(abind(predY2, along = 3), c(1,2), mean)) #phosphorus
EpredY3 <- data.frame(apply(abind(predY3, along = 3), c(1,2), mean)) #cons water loss
EpredY4 <- data.frame(apply(abind(predY4, along = 3), c(1,2), mean)) # nitrogen
EpredY5 <- data.frame(apply(abind(predY5, along = 3), c(1,2), mean)) #pesticide

# add in Hydrobasin ID (for plotting)
EpredY$HYBAS_ID <- XDataNew[,1]

# calculating the predicted change in prob. of occurrence in a species given an increase threat
# subtract the baseline probability of occurrence from our scenario-based prob. occurrence where we increase threat levels
EpredY$Australasian.Bittern_phosphorus_loading <- rescale(EpredY2$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1)) # rescale from 0-1 for easier visualisation
EpredY$Australasian.Bittern_cons_water_loss <- rescale(EpredY3$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1))
EpredY$Australasian.Bittern_nitrogen_loading <- rescale(EpredY4$Australasian.Bittern - EpredY$Australasian.Bittern, to = c(0,1))
EpredY$Musk.Duck_phosphorus_loading <- rescale(EpredY2$Musk.Duck - EpredY$Musk.Duck, to = c(0,1)) # rescale from 0-1 for easier visualisation
EpredY$Musk.Duck_cons_water_loss <- rescale(EpredY3$Musk.Duck - EpredY$Musk.Duck, to = c(0,1))
EpredY$Musk.Duck_nitrogen_loading <- rescale(EpredY4$Musk.Duck - EpredY$Musk.Duck, to = c(0,1))
EpredY$Mangrove.Honeyeater_nitrogen_loading <- rescale(EpredY$Mangrove.Honeyeater - EpredY4$Mangrove.Honeyeater, to = c(0,1)) # switch baseline becuase association is negative instead of positive
EpredY$Mangrove.Honeyeater_cons_water_loss <- rescale(EpredY3$Mangrove.Honeyeater - EpredY$Mangrove.Honeyeater, to = c(0,1)) 
EpredY$Chestnut.Teal_pesticide_loading <- rescale(EpredY5$Chestnut.Teal - EpredY$Chestnut.Teal, to = c(0,1)) 
EpredY$Chestnut.Teal_phosphorus_loading <- rescale(EpredY2$Chestnut.Teal - EpredY$Chestnut.Teal, to = c(0,1)) 
EpredY$Black.Swan_cons_water_loss <- rescale(EpredY3$Black.Swan - EpredY$Black.Swan, to = c(0,1)) 
EpredY$Black.Swan_phosphorus_loading <- rescale(EpredY2$Black.Swan - EpredY$Black.Swan, to = c(0,1)) 

# classify basins so that if has increased detection probability above 50th percentile, can detect it

aus.bittern <- EpredY %>% 
  dplyr::select(HYBAS_ID, Australasian.Bittern_phosphorus_loading:Australasian.Bittern_nitrogen_loading) %>% 
  mutate(Consumptive_Water_Loss = ifelse(Australasian.Bittern_cons_water_loss > quantile(EpredY$Australasian.Bittern_cons_water_loss, 0.5), 1, 0)) %>% 
  mutate(Phosphorus_Loading = ifelse(Australasian.Bittern_phosphorus_loading > quantile(EpredY$Australasian.Bittern_phosphorus_loading, 0.5), 1, 0)) %>% 
  mutate(Nitrogen_Loading = ifelse(Australasian.Bittern_nitrogen_loading > quantile(EpredY$Australasian.Bittern_nitrogen_loading, 0.5), 1, 0)) %>% 
  mutate(total = Phosphorus_Loading + Consumptive_Water_Loss + Nitrogen_Loading) %>% 
  filter(total == 1) %>% # filter where can only indicate one threat
  pivot_longer(cols = Phosphorus_Loading:Nitrogen_Loading, names_to = 'Threat', values_to = 'Detect') %>% 
  mutate(species = 'Australasian Bittern') %>% 
  dplyr::select(HYBAS_ID, species, Threat, Detect)

musk.duck <- EpredY %>% 
  dplyr::select(HYBAS_ID, Musk.Duck_phosphorus_loading, Musk.Duck_cons_water_loss, Musk.Duck_nitrogen_loading) %>% 
  mutate(Consumptive_Water_Loss = ifelse(Musk.Duck_cons_water_loss > quantile(EpredY$Musk.Duck_cons_water_loss, 0.5), 1, 0)) %>% 
  mutate(Phosphorus_Loading = ifelse(Musk.Duck_phosphorus_loading> quantile(EpredY$Musk.Duck_phosphorus_loading, 0.5), 1, 0)) %>% 
  mutate(Nitrogen_Loading = ifelse(Musk.Duck_nitrogen_loading > quantile(EpredY$Musk.Duck_nitrogen_loading, 0.5), 1, 0)) %>% 
  mutate(total = Phosphorus_Loading + Consumptive_Water_Loss + Nitrogen_Loading) %>% 
  filter(total == 1) %>% # filter where can only indicate one threat
  pivot_longer(cols = Phosphorus_Loading:Nitrogen_Loading, names_to = 'Threat', values_to = 'Detect') %>% 
  mutate(species = 'Musk Duck') %>% 
  dplyr::select(HYBAS_ID, species, Threat, Detect)

Mangrove.Honeyeater <- EpredY %>% 
  dplyr::select(HYBAS_ID, Mangrove.Honeyeater_cons_water_loss, Mangrove.Honeyeater_nitrogen_loading) %>% 
  mutate(Consumptive_Water_Loss = ifelse(Mangrove.Honeyeater_cons_water_loss > quantile(EpredY$Mangrove.Honeyeater_cons_water_loss, 0.5), 1, 0)) %>% 
  mutate(Nitrogen_Loading = ifelse(Mangrove.Honeyeater_nitrogen_loading > quantile(EpredY$Mangrove.Honeyeater_cons_water_loss, 0.5), 1, 0)) %>% 
  mutate(total = Consumptive_Water_Loss + Nitrogen_Loading) %>% 
  filter(total == 1) %>% # filter where can only indicate one threat
  pivot_longer(cols = Consumptive_Water_Loss:Nitrogen_Loading, names_to = 'Threat', values_to = 'Detect') %>% 
  mutate(species = 'Mangrove Honeyeater') %>% 
  dplyr::select(HYBAS_ID, species, Threat, Detect)

Chestnut.Teal <- EpredY %>% 
  dplyr::select(HYBAS_ID, Chestnut.Teal_pesticide_loading, Chestnut.Teal_phosphorus_loading) %>% 
  mutate(Pesticide_Loading = ifelse(Chestnut.Teal_pesticide_loading > quantile(EpredY$Chestnut.Teal_pesticide_loading, 0.5), 1, 0)) %>% 
  mutate(Phosphorus_Loading = ifelse(Chestnut.Teal_phosphorus_loading > quantile(EpredY$Chestnut.Teal_phosphorus_loading, 0.5), 1, 0)) %>% 
  mutate(total = Pesticide_Loading + Phosphorus_Loading) %>% 
  filter(total == 1) %>% # filter where can only indicate one threat
  pivot_longer(cols = Pesticide_Loading:Phosphorus_Loading, names_to = 'Threat', values_to = 'Detect') %>% 
  mutate(species = 'Chestnut Teal') %>% 
  dplyr::select(HYBAS_ID, species, Threat, Detect)

Black.Swan <- EpredY %>% 
  dplyr::select(HYBAS_ID, Black.Swan_cons_water_loss, Black.Swan_phosphorus_loading) %>% 
  mutate(Consumptive_Water_Loss = ifelse(Black.Swan_cons_water_loss > quantile(EpredY$Black.Swan_cons_water_loss, 0.5), 1, 0)) %>% 
  mutate(Phosphorus_Loading = ifelse(Black.Swan_phosphorus_loading > quantile(EpredY$Black.Swan_phosphorus_loading, 0.5), 1, 0)) %>% 
  mutate(total = Consumptive_Water_Loss + Phosphorus_Loading) %>% 
  filter(total == 1) %>% # filter where can only indicate one threat
  pivot_longer(cols = Consumptive_Water_Loss:Phosphorus_Loading, names_to = 'Threat', values_to = 'Detect') %>% 
  mutate(species = 'Black Swan') %>% 
  dplyr::select(HYBAS_ID, species, Threat, Detect)

# bind all

plot.df <- rbind(aus.bittern, musk.duck, Mangrove.Honeyeater, Chestnut.Teal, Black.Swan)

# map total indicator species in a basin

plot.df2 <- plot.df %>% 
  pivot_wider(names_from = 'species', values_from = 'Detect') 
plot.df2[is.na(plot.df2)] <- 0
plot.df2 <- plot.df2 %>% 
  mutate(num_indicator_species = `Australasian Bittern`+`Musk Duck`+`Mangrove Honeyeater`+`Chestnut Teal`+`Black Swan`) %>% 
  mutate(num_indicator_species = factor(num_indicator_species)) %>% 
  filter(num_indicator_species != 0)

predY.sf <- basins.qld %>% 
  inner_join(plot.df2, by = 'HYBAS_ID')

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(predY.sf) +
  tm_polygons('num_indicator_species', palette = 'Set2', legend.show = T) + # turn on legend by saying T
  tm_facets(by = 'Threat', ncol = 2, free.coords = F) +
  tm_layout(legend.position = c(0.01, 0.9),
            legend.title.size = 0.8,
            legend.text.size = 0.45)
m3
tmap_save(m3, 'outputs/map-scenario_number_indicator_species.png', width = 6, height = 6)

# map

predY.sf <- basins.qld %>% 
  inner_join(plot.df, by = 'HYBAS_ID') %>% 
  filter(Detect == 1)

m3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(predY.sf) +
  tm_polygons('Threat', palette = 'Set2', legend.show = T) + # turn on legend by saying T
  tm_facets(by = 'species', ncol = 3, free.coords = F) +
  tm_layout(legend.position = c(0.01, 0.9),
            legend.title.size = 0.8,
            legend.text.size = 0.45)
m3
tmap_save(m3, 'outputs/map-scenario_indicator_species.png', width = 6, height = 6)

### End here
