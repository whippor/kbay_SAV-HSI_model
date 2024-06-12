########################
### 91. HSI analyses ###
########################

# reproducibility
renv::restore()

# clear environment
rm(list = ls())

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, # is replacing the raster package
               viridis)

#####################################
#####################################

# load data
understorey <- terra::rast("data/d_suitability_data/understorey_HSI/understorey_HSI.tif")

canopy <- terra::rast("data/d_suitability_data/canopy_HSI/canopy_HSI.tif")

seagrass <- terra::rast("data/d_suitability_data/seagrass_HSI/seagrass_HSI.tif")

allmax <- terra::rast("data/d_suitability_data/all_HSI/allmax_HSI.tif")

bathymetry_under <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")

bathymetry_canop <- terra::rast("data/b_intermediate_data/canopy_bathymetry/bathymetry.grd")

bathymetry_seagr <- terra::rast("data/b_intermediate_data/seagrass_bathymetry/bathymetry.grd")


#####################################
#####################################

# understorey stats
# mask out all HSI below -30 m, above 3 m
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE)
understorey30 <- mask(understorey, 
                    bathmask)
plot(understorey30)

#understorey stats
terra::expanse(understorey) # 946.16 km^2 (total bay)
terra::hist(understorey)

terra::expanse(understorey30) # 373.91 km^2 (above -30)
terra::hist(understorey30)

## proportion suitability 
under_df <- data.frame(understorey30)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x=HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05))))

# overlap of HSI = 1 understorey with canopy HSI = 1
kelp_overlap <- c(canopy, understorey)
kelpDF <- data.frame(kelp_overlap)  
kelpDF %>%
  filter(HSI_value == 1) %>%
  tally(HSI_value.1 == 1) # 5186 cells
canopyDF <- data.frame(canopy)
canopyDF %>%
  filter(HSI_value == 1) %>%
  nrow() # 19430 cells
# 5186/19430  ~27%

#####################################
#####################################

# canopy stats
# mask out all HSI below -30 m, above 0 m
bathmask <- clamp(bathymetry_canop, lower = -30, upper = 0, values = FALSE)
canopy30 <- mask(canopy, 
                      bathmask)
plot(canopy30)

#canopy stats
terra::expanse(canopy) # 946.16 km^2 (total bay)
terra::hist(canopy)

terra::expanse(canopy30) # 352.60 km^2 (above -30)
terra::hist(canopy30)

## proportion suitability 
under_df <- data.frame(canopy30)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x=HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

# printed frequencies
as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05))))


#####################################
#####################################

#seagrass stats
# mask out all HSI below -10 m, above 1 m
bathmask <- clamp(bathymetry_seagr, lower = -10, upper = 1, values = FALSE)
seagrass10 <- mask(seagrass, 
                 bathmask)
plot(seagrass10)

#seagrass stats
terra::expanse(seagrass) # 946.16 km^2 (total bay)
terra::hist(seagrass)

terra::expanse(seagrass10) # 184.46 km^2 (above -10)
terra::hist(seagrass10)

## proportion suitability 
under_df <- data.frame(seagrass10)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x=HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05))))

4295/nrow(filter(test$HSI_value ==1))

#####################################
#####################################

# mask out all HSI below -30 m, above 3 m
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE)
allmax30 <- mask(allmax, 
                   bathmask)
plot(allmax30)

#combined stats
terra::expanse(allmax) # 946.16 km^2 (total bay)
terra::hist(allmax)

terra::expanse(allmax30) # 373.91 km^2 (above -30)
terra::hist(allmax30)

## proportion suitability 
under_df <- data.frame(allmax30)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x=HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()


as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05))))

max30DF <- data.frame(allmax30)
max30DF %>%
  filter(HSI_value > 0.5) %>%
  nrow() # 134436
# 134436/149562 0.8988647




plot(allmax, col = viridis(nrow(allmax), begin = 0.3))

plet(allmax,
     main = "SAV max HSI")


#####################################
#####################################

## NOTE!!!!
# must change all hist() to include maxcell = large number to get all values!

