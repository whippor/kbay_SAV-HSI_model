#########################
### 888. HSI analyses ###
#########################

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
               viridis,
               leaflet)

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

presence <- terra::rast("data/c_submodel_data/seagrass_presence_HSI/presenceHSI.grd")

#####################################
#####################################

# canopy stats
# mask out all HSI below -30 m, above 0 m
bathmask <- clamp(bathymetry_canop, lower = -30, upper = 0, values = FALSE)
canopy30 <- mask(canopy, 
                      bathmask)
plot(canopy30)
plet(canopy30,
     #tiles = "Stadia.AlidadeSmooth",
     main = "Canopy Kelp HSI") 

#canopy stats
terra::expanse(canopy) # 946.16 km^2 (total bay)
terra::hist(canopy)

terra::expanse(canopy30) # 352.59 km^2 (above -30)
terra::hist(canopy30, maxcell = 1100000)

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
cantab <- as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05),
                                  include.lowest = TRUE)))
# values between 0.10 - 0.65
sum(cantab[3:18,2]) # 93420
# total cells: 141036
# proportion 
# 93420/141036 # 0.662
allones <- under_df %>%
  filter(HSI_value == 1)

# 16527 canopy HSI 1 from presence

#####################################
#####################################

# understorey stats
# mask out all HSI below -30 m, above 3 m
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE)
understorey30 <- mask(understorey, 
                      bathmask)
plot(understorey30)
plet(understorey30,
     #tiles = "Stadia.AlidadeSmooth",
     main = "Understorey Kelp HSI")

#understorey stats
terra::expanse(understorey) # 946.16 km^2 (total bay)
terra::hist(understorey, maxcell = 1100000)

terra::expanse(understorey30) # 373.91 km^2 (above -30)
terra::hist(understorey30, maxcell = 1100000)

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

# printed frequencies
undtab <- as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05),
                                  include.lowest = TRUE)))

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

#seagrass stats
# mask out all HSI below -7 m, above 3 m
bathmask <- clamp(bathymetry_seagr, lower = -10, upper = 3, values = FALSE)
seagrass3 <- mask(seagrass, 
                 bathmask)
plot(seagrass3)
plet(seagrass3,
     #tiles = "Stadia.AlidadeSmooth",
     main = "Seagrass HSI")

#seagrass stats
terra::expanse(seagrass) # 946.16 km^2 (total bay)
terra::hist(seagrass)

terra::expanse(seagrass3) # 194.47 km^2 (3 to -10)
terra::hist(seagrass3, maxcell = 1100000)

## proportion suitability 
under_df <- data.frame(seagrass3)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x=HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

# printed frequencies
seatab <- as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05),
                                  include.lowest = TRUE)))
# values between 0.25 - 0.75
sum(seatab[6:15,2]) # 69207
# total cells: 77789
# proportion 
# 69207/77789 # 0.889
nonzero <- under_df %>%
  filter(HSI_value != 0)

# HSI 1 contributions
all1 <- data.frame(presence[presence == 1]) # 4295
seag1 <- data.frame(seagrass3[seagrass3 == 1]) #4292
# 4292/4295 # 0.9993

#####################################
#####################################

# mask out all HSI below -30 m, above 3 m
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE)
allmax30 <- mask(allmax, 
                   bathmask)
plot(allmax30)
plet(allmax30,
     #tiles = "Stadia.AlidadeSmooth",
     main = "All SAV HSI")

#combined stats
terra::expanse(allmax) # 946.16 km^2 (total bay)
terra::hist(allmax)

terra::expanse(allmax30) # 373.91 km^2 (above -30)
terra::hist(allmax30, maxcell = 1100000)

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


# printed frequencies
alltab <- as.data.frame(table(cut(under_df$HSI_value,breaks=seq(0,1,by=0.05),
                                  include.lowest = TRUE)))
allbay <- data.frame(understorey)
sum(alltab[2:20,2]) # 119611
# 119611/149562 # 0.799
# 119611/378464 # 0.3160


max30DF <- data.frame(allmax30)
max30DF %>%
  filter(HSI_value > 0.5) %>%
  nrow() # 114522
 114522/149562 # 0.765




plot(allmax, col = viridis(nrow(allmax), begin = 0.3))

plet(allmax,
     main = "SAV max HSI")


#####################################
#####################################

## NOTE!!!!
# must change all hist() to include maxcell = large number to get all values!

