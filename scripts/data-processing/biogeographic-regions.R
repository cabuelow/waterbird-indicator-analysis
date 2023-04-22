library(sf)
library(tidyverse)
library(tmap)

birds.env <- read.csv('data/waterbird-basin-threats_100.csv') 
basins.qld <- st_read('data/Enviro_dat_QLD/hydro_basins_qld_lev08_valid.shp')
bas <- basins.qld %>% filter(HYBAS_ID %in% birds.env$HYBAS_ID)
bio <- st_read('data/Biogeographic_Regions.shp') %>% st_transform(st_crs(basins.qld))

basbio <- bas %>% st_join(bio, left = F, largest = T)
miss <- bas %>% filter(!HYBAS_ID %in% basbio$HYBAS_ID)
tmap_mode('view')
qtm(bio) + qtm(bas) + qtm(miss, symbols.col = 'red') + qtm(st_buffer(miss, 200000), symbols.col = 'red')

# missing hybas id is in SEQ so manually assign, 5080070600

df <- bas %>% st_join(bio, left = T, largest = T) %>% st_drop_geometry()

df2 <- df %>% 
  select(HYBAS_ID, Q_REG, Q_REG_NAME) %>% 
  mutate(Q_REG = ifelse(HYBAS_ID == 5080070600, 'SEQ', Q_REG),
         Q_REG_NAME = ifelse(HYBAS_ID == 5080070600, 'Southeast Queensland', Q_REG_NAME))

birds.env2 <- birds.env %>% 
  left_join(df2)

write.csv(birds.env2, 'data/waterbird-basin-threats-100_bioregion.csv', row.names = F)
