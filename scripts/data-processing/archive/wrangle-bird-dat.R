# split up bird data

library(tidyverse)
library(lubridate)
library(sf)

dat <- read.csv('data/GBR_Birdata.csv', fileEncoding="UTF-8-BOM")

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

# subset by year and export as .csv

dat.split <- split(dat, format(dat$Start_Date, "%Y"))

for(i in seq_along(dat.split)){
  dat.year <- dat.split[[i]]
  year <- year(dat.year[1,]$Start_Date)
  write.csv(dat.year, paste0('data/GBR_Birdata_Year/dat-', year, '.csv'), row.names = F)
}

# save as spatial points

for(i in seq_along(dat.split)){
  dat.year <- dat.split[[i]]
  sf <- st_as_sf(dat.year, coords = c('Lon', 'Lat'), crs = "+proj=longlat +datum=WGS84")
  year <- year(dat.year[1,]$Start_Date)
  st_write(sf, paste0('data/GBR_Birdata_Year_Spatial/dat-', year, '.shp'))
}

