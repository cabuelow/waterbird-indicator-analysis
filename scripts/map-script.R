# map of study region and bird and threat variables

library(tmap)
library(tidyverse)
library(sf)

# read in data 

basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
birds.env <- read.csv('data/waterbird-basin-threats-100_bioregion.csv') 
bio <- st_read('data/Biogeographic_Regions.shp') %>% st_transform(st_crs(basins.qld))
qld <- st_read('data/qld-shp/queensland-polygon.shp')

# filter for basins with birds

basins.sub <- basins.qld %>% 
  inner_join(birds.env, by = 'HYBAS_ID')
  
# take a quick look

qtm(basins.sub, 'Aquaculture_Pressure')

# make a column of birds species richness and cumulative threats

basins.sub$Species.richness <- rowSums(dplyr::select(st_drop_geometry(basins.sub),
                                              Australasian.Darter:Musk.Duck))
# update threats selected below to only include those used in final model
basins.sub$Cumulative.threat <- rowSums(dplyr::select(st_drop_geometry(basins.sub),
                                               Agricultural_Water_Stress:Wetland_Disconnectivity))

# select bioregions for our study

bio.sub <- bio %>% 
  filter(Q_REG_NAME %in% birds.env$Q_REG_NAME) %>% 
  st_simplify(1000, preserveTopology = T)

# map bioregions

map <- tm_shape(qld) +
  tm_fill() +
  tm_shape(bio.sub) +
  tm_compass() +
  tm_scale_bar() +
  tm_polygons('Q_REG_NAME') +
  tm_layout(#legend.position = c('left', 'bottom'),
            legend.title.size = 0.1,
            legend.text.size = 0.3,
            bg.color = 'white')

tmap_save(map, 'outputs/study-area-map.png',  height = 4, width = 3)

# mapping threats and species

map1 <- tm_shape(qld) +
  tm_fill() +
tm_shape(basins.sub) +
  tm_polygons('Species.richness', title = 'Species Richness') +
  #tm_compass() +
  #tm_scale_bar() +
  tm_layout(legend.position = c('left', 'bottom'),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white')
map1

map2 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(basins.sub) +
  tm_polygons('Cumulative.threat', title = 'Cumulative Threat') +
  tm_compass(position = c(0.7, 0.12)) +
  tm_scale_bar(position = c(0.5, 0.06))  +
  tm_layout(legend.position = c('left', 'bottom'),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white')
map2

# save

map3 <- tmap_arrange(map1, map2)
map3

tmap_save(map3, 'outputs/study-area-map2.png',  height = 4, width = 6)
