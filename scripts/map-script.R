# map of study region and bird and threat variables

library(tmap)
library(tidyverse)
library(sf)
library(terra)
library(tmaptools)
library(cowplot)
library(ggplot2)

# read in data 

data(World)
aus <- World %>% filter(name == 'Australia')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
birds.env <- read.csv('data/waterbird-basin-threats-100_bioregion_notshore.csv') 
bio <- st_read('data/Biogeographic_Regions.shp') %>% st_transform(st_crs(basins.qld))
qld <- st_read('data/qld-shp/queensland-polygon.shp')
fils <- list.files('data/Enviro_dat_QLD/river-threats', pattern = '.tif', full.names = T)
threats <- lapply(fils, rast)

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
                                                      Aquaculture_Pressure, Consumptive_Water_Loss,
                                                      Nitrogen_Loading, Pesticide_Loading, Phosphorus_Loading))

# select bioregions for our study

bio.sub <- bio %>% 
  filter(Q_REG_NAME %in% birds.env$Q_REG_NAME) %>% 
  st_simplify(1000, preserveTopology = T)

# map bioregions

map <- tm_shape(qld) +
  tm_fill() +
  tm_shape(bio.sub) +
  tm_fill('Q_REG_NAME', title = 'Bioregion', legend.show = T) +
  tm_shape(basins.sub) +
  tm_polygons(alpha = 0) +
  tm_compass(position = c(0.01, 0.08)) +
  tm_scale_bar(position = c(0.01, 0.02)) +
  tm_layout(legend.position = c(0.5, 0.8),
            #legend.title.size = 0.1,
            legend.text.size = 0.6,
            legend.width = 2,
            bg.color = 'white',
            frame = F)
map

tmap_save(map, 'outputs/study-area-map_bioregions-basins.png',  height = 4, width = 3)

# make an inset with australia

qld.rect <- bb_poly(qld)

ausinset <- tm_shape(aus) +
  tm_fill(col = 'lightgrey') +
  tm_shape(qld.rect) +
  tm_borders(lw = 2) +
  tm_layout(frame = T)
ausinset

# add to inset to map
sitemap_g <- tmap_grob(map)
ausinset_g <- tmap_grob(ausinset)

finalmap <- ggdraw() +
  draw_plot(sitemap_g) +
  draw_plot(ausinset_g,
            width = 0.2, height = 0.2,
            x = 0.5, y = 0.04)
finalmap

ggsave('outputs/study-area-map_bioregions-basins_inset.png', width = 5, height = 5)

# mapping threats and species

map1 <- tm_shape(qld) +
  tm_fill() +
tm_shape(basins.sub) +
  tm_polygons('Species.richness', title = 'Species Richness') +
  tm_compass(position = c(0.01, 0.08)) +
  tm_scale_bar(position = c(0.01, 0.02)) +
  tm_layout(legend.position = c(0.6, 0.80),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map1

tmap_save(map1, 'outputs/study-area-spp-richness.png',  height = 4, width = 3)

map2 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(basins.sub) +
  tm_fill('Cumulative.threat', title = 'Cumulative Threat') +
  tm_compass(position = c(0.01, 0.08)) +
  tm_scale_bar(position = c(0.01, 0.02)) +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map2
tmap_save(map2, 'outputs/study-area-cumulative-threat_v2.png',  height = 4, width = 3)

# plot threats

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[2]]) +
  tm_raster(title = 'Aquaculture Pressure') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/aqua-pressure-map.png',  height = 4, width = 3)

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[2]]) +
  tm_raster(title = 'Aquaculture Pressure') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/aqua-pressure-map.png',  height = 4, width = 3)

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[3]]) +
  tm_raster(title = 'Consumptive water loss') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/consumptive-water-loss-map.png',  height = 4, width = 3)

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[12]]) +
  tm_raster(title = 'Nitrogen loading') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/nitrogen-map.png',  height = 4, width = 3)

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[16]]) +
  tm_raster(title = 'Pesticide loading') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/pesticide-map.png',  height = 4, width = 3)

map3 <- tm_shape(qld) +
  tm_fill() +
  tm_shape(threats[[17]]) +
  tm_raster(title = 'Phosphorus loading') +
  tm_layout(legend.position = c(0.6, 0.70),
            legend.title.size = 0.8,
            legend.text.size = 0.6,
            bg.color = 'white',
            frame = F)
map3
tmap_save(map3, 'outputs/phosphorus-map.png',  height = 4, width = 3)

