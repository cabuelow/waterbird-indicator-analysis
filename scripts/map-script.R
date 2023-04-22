# map of study region and bird and threat variables

library(tmap)
library(tidyverse)
library(sf)
library(rnaturalearth)

# read in data 

basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
birds.env <- read.csv('data/waterbird-basin-threats_100.csv') 
aus <- countries110 %>% st_as_sf() %>% filter(ADMIN == 'Australia') 
qtm(aus)

# filter for basins with birds

basins.sub <- basins.qld %>% 
  inner_join(birds.env, by = 'HYBAS_ID')
  
# take a quick look

qtm(basins.sub, 'Aquaculture_Pressure')
#qtm(basins.sub, factor('Musk.Duck')) # need to figure why doesn't like this  binary variable

# crop aus for study area

qld <- aus %>% st_crop(st_bbox(basins.sub))
qtm(qld)

# make a column of birds species richness and cumulative threats

basins.sub$Species.richness <- rowSums(select(st_drop_geometry(basins.sub),
                                              Australasian.Darter:Musk.Duck))
# update threats selected below to only include those used in final model
basins.sub$Cumulative.threat <- rowSums(select(st_drop_geometry(basins.sub),
                                               Agricultural_Water_Stress:Wetland_Disconnectivity))

# publication quality map

map1 <- tm_shape(qld) +
  tm_polygons() +
tm_shape(basins.sub) +
  tm_polygons('Species.richness') +
  #tm_compass() +
  #tm_scale_bar() +
  tm_layout(legend.position = c('left', 'bottom'),
            bg.color = 'lightblue')
map1

map2 <- tm_shape(qld) +
  tm_polygons() +
  tm_shape(basins.sub) +
  tm_polygons('Cumulative.threat') +
  tm_compass() +
  tm_scale_bar() +
  tm_layout(legend.position = c('left', 'bottom'),
            bg.color = 'lightblue')
map2

# save

map3 <- tmap_arrange(map1, map2)
map3

tmap_save(map3, 'study-area-map.png')
