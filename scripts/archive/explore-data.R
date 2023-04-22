library(tidyverse)
library(sf)
library(tmap)

# data

dat <- read.csv('data/waterbird-basin-threats-comp-cases.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')

# convert to basins for plotting

basins <- basins.qld %>% 
  left_join(dat, by = 'HYBAS_ID')

# plot

tmap_mode('view')

tm_shape(basins) + 
  tm_polygons(col = 'Agricultural_Water_Stress')

tm_shape(basins) + 
  tm_polygons(col = 'Australasian.Darter')
