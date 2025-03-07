#########################
### 721. HSI combined ###
#########################

# clear environment
rm(list = setdiff(ls(), "master_begin"))

# calculate start time of code (determine how long it takes to complete all code)
all_begin <- Sys.time()

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

### output directories
#### final directory
dir.create("data/d_suitability_data/all_HSI")
submodel_dir <- "data/d_suitability_data/all_HSI"

#####################################
#####################################

# load data
understorey <- terra::rast("data/d_suitability_data/understorey_HSI/understorey_HSI.tif")

canopy <- terra::rast("data/d_suitability_data/canopy_HSI/canopy_HSI.tif")

seagrass <- terra::rast("data/d_suitability_data/seagrass_HSI/seagrass_HSI.tif")

#####################################
#####################################

# combine and choose max values for each cell
allmax <- terra::mosaic(understorey, 
                        canopy, 
                        seagrass, 
                        fun = "max")

plot(allmax, col = viridis(nrow(allmax), begin = 0.3))

plet(allmax,
     main = "SAV max HSI")

# Export data
## Suitability
terra::writeRaster(allmax, 
                   filename = file.path(submodel_dir, "allmax_HSI.tif"), 
                   overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



