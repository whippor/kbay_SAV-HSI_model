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
  
#####################################
#####################################

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

#####################################
#####################################

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










plot(allmax, col = viridis(nrow(allmax), begin = 0.3))

plet(allmax,
     main = "SAV max HSI")


#####################################
#####################################



