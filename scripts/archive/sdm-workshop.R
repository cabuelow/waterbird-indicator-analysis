# species distribution modelling
# 2020-08-09
# CABuelow

# libraries

library(dplyr) # data wrangling
library(dismo) # spp distribution modelling
library(sf) # spatial vector data (points, polygons, lines)
library(lubridate) # formatting dates and times
library(tmap) # making maps

# data

dat <- read.csv('data/GBR_Birdata.csv', fileEncoding="UTF-8-BOM")
basins <- st_read('data/Enviro_dat_QLD/hydro_basins_qld.shp')
#rivs <- st_read('data/Enviro_dat_QLD/hydro_rivers_qld.shp')
env <- stack('data/bioclim.data.crop.tif')

env.fils <- list.files('data/Enviro_dat_QLD/river-threats/', full.names = T)
riv.threat <- stack(env.fils)

# determine geographic extent of data

max.lat <- ceiling(max(dat$Lat))
min.lat <- floor(min(dat$Lat))
max.lon <- ceiling(max(dat$Lon))
min.lon <- floor(min(dat$Lon))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# get enviro data

#bioclim.data <- getData(name = "worldclim",
#                      var = "bio",
#                    res = 2.5,
#                  path = "data/")
#writeRaster(bioclim.data, 'data/bioclim.data.tif')

# crop bioclim data to geographic extent

#bioclim.data.crop <- crop(x = bioclim.data, y = geographic.extent)
#writeRaster(bioclim.data.crop, 'data/bioclim.data.crop.tif')

# format dates

dat$Start_Date <- mdy(dat$Start_Date)
dat$Start_Date_DateFormatted <- mdy(dat$Start_Date_DateFormatted)
dat$Start_Date_T2C <- dmy(dat$Start_Date_T2C)
dat$Start_Date_T2C_DateFormatted <- dmy(dat$Start_Date_T2C_DateFormatted)
dat$Finish_Date <- mdy(dat$Finish_Date)
dat$Finish_Date_DateFormatted <- mdy(dat$Finish_Date_DateFormatted)
dat$Finish_Date_T2C <- dmy(dat$Finish_Date_T2C)
dat$Finish_Date_T2C_DateFormatted <- dmy(dat$Finish_Date_T2C_DateFormatted)
str(dat)

# add a column for year

dat$Year <- year(dat$Start_Date)

# lets look at the magpie goose

target <- c(1995:2005) # target years between 1995 and 2005

mag <- dat %>% filter(Common.Name == 'Magpie Goose') %>% # filter for magpie goose between years 1995 and 2005
  filter(Year %in% target)

mag.sf <- st_as_sf(mag, coords = c('Lon', 'Lat'), crs = "+proj=longlat +datum=WGS84") # make into spatial dataframe for plotting 

# plot to check

tmap_mode(mode = 'view') # can use 'plot' for static maps or 'view' mode interactive maps

tm_shape(mag.sf) + # we layer on the shape and then add layer for how to plot it (e.g. dots vs. poygons vs. rasters)
  tm_dots()

# select only occurrences that intersect with basins

mag.sf.sub <- mag.sf %>% 
  st_join(basins, left = F) 

basins.sub <- basins %>% 
  st_join(mag.sf, left = F) %>% 
  group_by(HYBAS_ID) %>% 
  summarise()

# plot to check

tm_shape(mag.sf) +
  tm_dots() +
  tm_shape(mag.sf.sub) +
  tm_dots(col = 'red')

tm_shape(basins.sub) +
  tm_polygons() +
  tm_shape(mag.sf.sub) +
  tm_dots(col = 'red')

tm_shape(subset(env, 15)) + #can look at different env varibles by changing number
  tm_raster() +
  tm_shape(mag.sf.sub) +
  tm_dots(col = 'red')

tm_shape(subset(riv.threat, 1)) + # can look at different river threats by changing number
  tm_raster() +
  tm_shape(mag.sf.sub) +
  tm_dots(col = 'red')


# extract values of bioclim data at occurrence points

# check for colinearity of predictors

# use bioclim function to build sdm 
# p 46 of vignette

# predict

# model evaluation


