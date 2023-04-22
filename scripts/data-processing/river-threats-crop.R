# explore river threats
# CABuelow
# 17-08-2020

# libraries

library(tidyverse)
library(ncdf4)
library(raster)
library(tmap)

# import data

river.fils <- list.files('data/RiverThreatDrivers/', full.names = T)

# extent of GBR bird data

max.lat <- -9
min.lat <- -32
max.lon <- 155
min.lon <- 137
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# open, crop, save nc files

for(i in 1:length(river.fils)){
  
nc <- nc_open(river.fils[i]) # use netcdf4 to expore netcdf attributes
brick <- brick(river.fils[i])

# crop

brick.crop <- crop(brick, geographic.extent)

# write raster

varname <- names(nc$var)
writeRaster(brick.crop, paste0('data/Enviro_dat_QLD/river-threats/', varname, '.tif'), overwrite = T)
}

