###########################################
### 05. Understorey Substrate submodel ####
###########################################

# clear environment
rm(list=setdiff(ls(), c("all_begin", "master_begin")))

# calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, # is replacing the raster package
               viridis)

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### roi directory
roi_dir <- "data/b_intermediate_data/roi"

### output directories
#### submodel directory
dir.create("data/c_submodel_data/understorey_substrate_HSI")
submodel_dir <- "data/c_submodel_data/understorey_substrate_HSI"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect(roi_dir)

#####################################
#####################################

# load data
substrate <- terra::rast("data/b_intermediate_data/understorey_substrate/substrate.tif")

tam_subs <- read_csv("data/x_tam_tables/understorey/understorey_substrate.csv")

#####################################
#####################################

# Create bathymetry HSI model
# mask bathymetry to the roi
subs_mask <- mask(substrate, roi)

# plot to check for correct values
# plot(subs_mask, col = viridis(nrow(subs_mask)))

# extract all values from subs roi
vals1 <- data.frame(values(subs_mask))
vals1 <- vals1 %>%
  rename("value" = "substrate")

# ASSIGN HSI VALUES TO EACH SUBSTRATE CLASS
## Check value <-> level associations
subs_values <- data.frame(levels(subs_mask))
### 0 Boulder
### 1 Cobble/Pebble
### 2 Land
### 3 Mud
### 4 Sand
### 5 Unclassified

## join tam and values tables
tam_subs_new <- tam_subs %>%
  rename("substrate" = "substrate.class")
join_ID <- tam_subs_new %>%
  full_join(subs_values) 

# join HSI and raster values 
index_vals <- vals1 %>%
  left_join(join_ID) %>%
  select(substrate.class.SIV)

# join HSI values with raster
subs_mask[["HSI_value"]] <- index_vals

# check plot
# check plot
plot(subs_mask, col = viridis(nrow(subs_mask), begin = 0.3))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(subs_mask, 
                   filename = file.path(submodel_dir, "substrateHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



