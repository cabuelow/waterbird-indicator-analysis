# crop enviro data to QLD

library(sf)
library(mapview)

qld_basins <- st_read('data/Drainage Basins/drainage-basins.shp')
crs.geo <- "+proj=longlat +datum=WGS84"
qld_basins <- st_transform(qld_basins, crs = 4326)
#mapview(qld_basins)
#leafem::addMouseCoordinates(qld_basins)

# qld bounding box for cropping

bbox <- st_bbox(qld_basins)

bbox <- st_bbox(c(xmin = 490406, xmax = 1000008, ymax = 8264279, ymin = 7667647), crs = st_crs(32760))

bbox <- st_bbox(st_transform(st_as_sfc(bbox), crs = 4326))

b <- st_as_sfc(bbox)

fils <- list.files('data/Basins', full.names = T, pattern = '.shp')
fils.names <- list.files('data/Basins', full.names = F, pattern = '.shp')

for(i in 1:length(fils)){
hydro_bas.12 <- st_read(fils[i]) %>%  
  st_transform(crs = 32760) %>% 
  st_make_valid() %>% st_crop(bbox)
st_write(hydro_bas.12.crop, paste0('data/', fils.names[i]))
}

hydro_lak <- st_read('data/HydroAtlas_QLD/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp')
hydro_lak.crop <- st_crop(hydro_lak, bbox)
st_write(hydro_lak.crop, 'data/Enviro_dat_QLD/hydro_lakes_qld.shp', overwrite = T)

hydro_riv <- st_read('data/HydroAtlas_QLD/HydroRIVERS_v10_shp/HydroRIVERS_v10.shp')
hydro_riv.crop <- st_crop(hydro_riv, bbox)
st_write(hydro_riv.crop, 'data/Enviro_dat_QLD/hydro_rivers_qld.shp', overwrite = T)

qld_wetl <- st_read('data/DP_QLD_WETLAND_SYSTEM_100K_A/SWS_QLD_WETLAND_AREAS.shp')
qld_wetl.crop <- st_crop(qld_wetl, bbox)
st_write(qld_wetl.crop, 'data/Enviro_dat_QLD/wetland_areas_qld.shp', overwrite = T)

