#######################################
### 211. Canopy Bathymetry submodel ###
#######################################

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
dir.create("data/c_submodel_data/canopy_bathymetry_HSI")
submodel_dir <- "data/c_submodel_data/canopy_bathymetry_HSI"

#####################################
#####################################

# set parameters

# define vector for region of interest
roi <- terra::vect(roi_dir)

#####################################
#####################################

# load data
bathymetry <- terra::rast("data/b_intermediate_data/bathymetry/bathymetry.grd")

tam_bath <- read_csv("data/x_tam_tables/canopy/canopy_depth.csv")
tam_bath <- tam_bath %>%
  mutate(depth.m = depth.m + 1.546) %>% # correct for NAVD88/MLLW offset
  arrange(depth.m)

#####################################
#####################################

# Create bathymetry HSI model
# mask bathymetry to the roi
bath_mask <- mask(bathymetry, roi)

# extract all values from bath_roi
vals1 <- data.frame(values(bath_mask))

# FUNCTION TO FIND Y FOR ANY GIVEN X WITH IMPORTED SLOPES

## Calculate slopes
slopes <- diff(tam_bath$depth.m.SIV) / diff(tam_bath$depth.m)

# calculate index from raster values with function
index_vals <- interpolate_y(vals1$KBL.bathymetry_GWA.area_50m_EPSG3338, tam_bath)

# join HSI values with raster
bath_mask[["HSI_value"]] <- index_vals

# check plot
plot(bath_mask, col = viridis(nrow(bath_mask), begin = 0.3))  

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(bath_mask, 
                   filename = file.path(submodel_dir, "bathymetryHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



