#####################################
### 214. Canopy Velocity submodel ###
#####################################

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
source("code/000_function_interpolate_y.R")

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### roi directory
roi_dir <- "data/b_intermediate_data/roi"

### output directories
#### submodel directory
dir.create("data/c_submodel_data/canopy_velocity_HSI")
submodel_dir <- "data/c_submodel_data/canopy_velocity_HSI"

#####################################
#####################################

# set parameters

# define vector for region of interest
roi <- terra::vect(roi_dir)

#####################################
#####################################

# load data
velocity <- terra::rast("data/b_intermediate_data/velocity/velocity.grd")

tam_velo <- read_csv("data/x_tam_tables/canopy/canopy_velocity.csv")
tam_velo <- tam_velo %>%
  arrange(velocity.ms)

#####################################
#####################################

# Create velocity HSI model
# mask velocity to the roi
velo_mask <- mask(velocity, roi)

# extract all values from bath_roi
vals1 <- data.frame(values(velo_mask))

# FUNCTION TO FIND Y FOR ANY GIVEN X WITH IMPORTED SLOPES

## Calculate slopes
slopes <- diff(tam_velo$velocity.ms.SIV) / diff(tam_velo$velocity.ms)

# calculate index from raster values with function
index_vals <- interpolate_y(vals1$maxSpeed_CIOFS500, tam_velo)

# join HSI values with raster
velo_mask[["HSI_value"]] <- index_vals

# check plot
plot(velo_mask, col = viridis(nrow(velo_mask), begin = 0.3))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(velo_mask, 
                   filename = file.path(submodel_dir, "velocityHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



