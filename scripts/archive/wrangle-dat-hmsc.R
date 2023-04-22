# wrangle bird dat for hmsc

# libraries

library(tidyverse) # data wrangling
library(sf) # spatial vector data (points, polygons, lines)
library(lubridate) # formatting dates and times
library(tmap) # making maps
library(raster) # spatial raster data
library(GGally)

# data

dat <- read.csv('data/GBR_Birdata.csv', fileEncoding="UTF-8-BOM")
dat.sub <- read.csv('outputs/GBR_Birdata_sub.csv')
ind.spp <- read.csv('data/waterbird-indicator-species.csv')
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
env <- stack('data/bioclim.data.crop.tif') # load stack of biocimatic variables
riv.fils <- list.files('data/Enviro_dat_QLD/river-threats/', full.names = T) # list all river threat data
riv.threat <- stack(riv.fils) # load individual river threat rasters as a stack of river threat rasters, 'stack' function

# determine geographic extent of GBR bird data

max.lat <- ceiling(max(dat$Lat))
min.lat <- floor(min(dat$Lat))
max.lon <- ceiling(max(dat$Lon))
min.lon <- floor(min(dat$Lon))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

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

# select target years 1995-2005 and waterbird indicator species
# NOTE: hashed out becuase we just read in the subsetted bird data now that we've saved it 
# as 'GBR_Birdata_sub.csv'

#years <- c(1995:2005) # target years between 1995 and 2005
#species <- ind.spp$Scientific.name

#dat.sub <- dat %>% 
 # filter(Scientific.Name == 'Ardea alba') %>% 
  #filter(dat.sub, Year %in% c(1995:2005))

#dat.sub <- dat %>% 
 # filter(Scientific.Name %in% species) %>% 
#  filter(Year %in% years)

#write.csv(dat.sub, 'outputs/GBR_Birdata_sub.csv', row.names = F)

# see which indicator spp are missing

missing <- ind.spp %>% 
  filter(!Scientific.name %in% dat.sub$Scientific.Name)
missing

# data wrangling for exploratory plotting

dat.sub$obs <- c(rep(1, nrow(dat.sub))) # make a column called 'obs' with a '1' for each row

dat.sub2 <- dat.sub %>% # calculate total number of observations for each year and species
  group_by(Year, Common.Name) %>% 
  summarise(total.obs = sum(obs))

dat.sub3 <- dat.sub2 %>%  # calculate average, standard deviation, and total number of observations for each species
  group_by(Common.Name) %>% 
  summarise(avg.year = mean(total.obs), 
            sd = sd(total.obs), 
            total = sum(total.obs))
  
ggplot(dat.sub2) + # make plot - ggplot layers on different plot attributes with the '+' sign
  geom_bar(aes(x = Year, y = total.obs), stat = 'identity') +
  facet_wrap(~Common.Name, scales = 'free') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45))

ggsave('outputs/spp-obs-year.png', width = 11, height = 11) # save plot as a .png in outputs folder

# make into subsetted bird data into a spatial dataframe

dat.sub.sf <- st_as_sf(dat.sub, coords = c('Lon', 'Lat'), 
                       crs = "+proj=longlat +datum=WGS84")  # proj4 string for WGS84 geographic coordinate system

# join bird occurrence spatial dataframe with hydrobasin sampling units

birds <- dat.sub.sf %>% 
 st_join(basins.qld, left = F) # do a 'spatial join' with st_join, left = F so that we do an inner join

# save the bird-hydrobasin spatial dataframe as a shapefile - can explore in arcgis
st_write(birds, 'data/waterbirds_basins_1995-2005.shp', append = F, overwrite = T)

# add observations column to bird spatial dataframe so can do summaries and exploratory plots below

dat.sub.sf$obs <- c(rep(1, nrow(dat.sub.sf)))

# join hydrobasins with bird data, group by hydrobasin and species name
# then calcualte the number of ovservations in each hydrobasin for each species

basins <- basins.qld %>% 
  st_join(dat.sub.sf, left = F) %>% 
  group_by(HYBAS_ID, Common.Name) %>% 
  summarise(num.obs = sum(obs, na.rm = T)) # group_by and summarise aggregates all of the data to hydrobasin (i.e. gets rid of multiple basins that occur in spatial join, because there are multiple occurresnces in a basin)

# save as shapefile to explore in arcgis
#st_write(basins, 'data/basins_waterbirds.shp', append = F, overwrite = T)

# or plot with tmap in R

# set the tmap mode to view for interactive plotting
tmap_mode(mode = 'view') # can use 'plot' for static maps or 'view' mode interactive maps

# exploratory maps

# plot number of observations in each basin, and overlay all bird occurrences
# tmap is like ggplot that we layer plots on top of one another with '+' sign
# use 'tm_shape' to give the spatial data want to plot, and then tm_polygons or tm_dots
# to map as either polygons or points

tm_shape(basins) +
  tm_polygons(col = 'num.obs') +
  tm_shape(dat.sub.sf) +
  tm_dots()

# plot number of observations for individual species

tm_shape(filter(basins, Common.Name == 'Eastern Reef Egret')) + # can look at different species by updating species name
  tm_polygons(col = 'num.obs') +
  tm_shape(dat.sub.sf) +
  tm_dots()

# summarize river threats in basins
# Note: this is hashed out becuase takes 31 minutes to process all of the threat rasters
# but here we calculate the average threat value in each hydrobasin
# we save the output as a .csv below
# save the code so can reproduce

#basins.sp <- as(basins, 'Spatial')

#list <- list()

#system.time( # takes 31 mins
#for(j in 1:nlayers(riv.threat)){
 # rast <- subset(riv.threat, j)
  #df <- data.frame(Threat = NA, HYBAS_ID = NA, mean_threat = NA)
#for(i in 1:nrow(basins.sp)){
 # basin <- basins.sp[i,]
  #threat.mean <- raster::extract(rast, basin, fun = mean, na.rm = T)
  #df[i,1] <- names(rast)
  #df[i,2] <- basin$HYBAS_ID
  #df[i,3] <- threat.mean
#}
 # list[[j]] <- df
#})

#basin.riv.thrt <- do.call(rbind, list)

# pivot river threats to wide format

#basin.riv.thrt.wide <- basin.riv.thrt %>% 
 # pivot_wider(names_from = Threat, values_from = mean_threat)

#write.csv(basin.riv.thrt.wide, 'data/basins-mean-riv-threat.csv', row.names = F)

# read in the threat data calculated and saved above

basin.riv.thrt.wide <- read.csv('data/basins-mean-riv-threat.csv')

# convert birds dataframe to presence-absence in hydrobasins

birds.df <- data.frame(birds) # make the birds spatial dataframe into a normal dataframe

rows <- nrow(birds.df) # calcualte the number of rows in the birds dataframe (will use below)

# calculate the number of occurrences in each basin for each species

birds.occur <- birds.df %>% 
  mutate(Occurrence = rep(1, rows)) %>% 
  group_by(HYBAS_ID, Common.Name, Scientific.Name) %>% 
  summarise(Num_occur = sum(Occurrence, na.rm = T))

birds.occur.total <- birds.df %>% 
  mutate(Occurrence = rep(1, rows)) %>% 
  group_by(Common.Name, Scientific.Name) %>% 
  summarise(Num_occur = sum(Occurrence, na.rm = T))

# convert to wide format

birds.pa <- birds.occur %>% 
  dplyr::select(-Scientific.Name) %>% 
  pivot_wider(names_from = Common.Name, values_from = Num_occur)

# check percent occurrence for each species, remove rare species

tal <- data.frame(birds.pa %>% 
                    mutate(obs = rep(1, nrow(.))) %>% 
                    filter(cpue > 0) %>% 
                    group_by(spp) %>% 
                    summarise(perc.occur = (sum(obs)/nrow(.))*100))
tal

# convert to presence absence 

birds.pa[is.na(birds.pa)] <- 0
birds.pa[,-1][birds.pa[,-1] >1 ] <- 1
  
# join river threats with waterbird presence abcsence data by HYBAS_ID

birds.riv.thrt <- birds.pa %>% 
  left_join(basin.riv.thrt.wide, by = 'HYBAS_ID')

# join relevant covariates from basins and basin centroids

centroids <- basins.qld %>%
  st_transform(crs = 3112) %>% 
  st_centroid()

cent.coords <- data.frame(centroids %>% 
  st_coordinates())

cent.coords$HYBAS_ID <- centroids$HYBAS_ID

birds.env <- birds.riv.thrt %>% 
  left_join(dplyr::select(basins.qld, HYBAS_ID, inu_pc_slt, ria_ha_ssu, ele_mt_sav,
                   tmp_dc_syr, pre_mm_syr,
                   crp_pc_sse, ire_pc_sse, urb_pc_use, hft_ix_s93), by = 'HYBAS_ID') %>% 
  left_join(cent.coords, by = 'HYBAS_ID') %>% dplyr::select(-geometry)

#write.csv(birds.env, 'data/waterbird-basin-threats.csv', row.names = F)

birds.env <- read.csv('data/waterbird-basin-threats.csv')

birds.comp <- birds.env[complete.cases(birds.env),]

# look for correlations in predictor variables

corr_check <- function(Dataset, threshold){
  matriz_cor <- cor(Dataset)
  matriz_cor
  
  for (i in 1:nrow(matriz_cor)){
    correlations <-  which((abs(matriz_cor[i,i:ncol(matriz_cor)]) > threshold) & (matriz_cor[i,i:ncol(matriz_cor)] != 1))
    
    if(length(correlations)> 0){
      lapply(correlations,FUN =  function(x) (cat(paste(colnames(Dataset)[i], "with",colnames(Dataset)[x]), "\n")))
      
    }
  }
}

colnames(birds.env)
corr_check(birds.comp[,52:83], 0.85)

# get rid of flow disruption, invasive spp threats, and highly correlated variables (>0.85) then only select complete cases

birds.env2 <- birds.env %>% 
  dplyr::select(-Flow_Disruption, -Non.Native_Fishes_number, 
                -Non.Native_Fishes_percent, -River_Fragmentation,
                -Consumptive_Water_Loss, -Human_Water_Stress, 
                -Impervious_Surfaces, -Fishing_Pressure, -Cropland,
                -Mercury_Deposition, -Organic_Loading)

birds.comp2 <- birds.env2[complete.cases(birds.env2),]

#write.csv(birds.comp2, 'data/waterbird-basin-threats-comp-cases.csv', row.names = F)

colnames(birds.comp2)
corr_check(birds.comp2[,52:72], 0.85)

corr.plot <- ggpairs(birds.comp, columns = 52:72) # threats
corr.plot

# check basins with missing threat rasters

miss.basins <- birds.env2 %>% 
  filter(!HYBAS_ID %in% birds.comp2$HYBAS_ID)

miss.basins.sf <- basins.qld %>% 
  filter(HYBAS_ID %in% miss.basins$HYBAS_ID)

# plot to check

tmap_mode(mode = 'view') # can use 'plot' for static maps or 'view' mode interactive maps

# plot to check

tm_shape(subset(riv.threat, 19)) + # can look at different river threats by changing number
  tm_raster() +
  tm_shape(miss.basins.sf) +
  tm_polygons(alpha = 0.5)
