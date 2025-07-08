##################
### 003. Fetch ###
##################

######### NOT WORKING

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
source("code/000_function_get_fetch.R")
source("code/000_function_get_landwater.R")
source("code/000_function_prep_raster.R")

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### input directories
data_dir <- "data/a_raw_data/bathymetry"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### fetch directory
dir.create(paste0(intermediate_dir, "/",
                  "fetch"))
fetch_dir <- "data/b_intermediate_data/fetch"

#### roi directory
roi_dir <- "data/b_intermediate_data/roi"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
# crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect(roi_dir)
# roi <- project(roi, crs)


#####################################
#####################################

# load data
bathymetry <- terra::rast("data/a_raw_data/KBL-bathymetry_GWA-area_50m_EPSG3338.tiff")

# inspect the data
## coordinate reference system
terra::crs(bathymetry) # EPSG:4269

units(bathymetry) <- "meters"

# mask and crop bathymetry
bath_mask <- mask(bathymetry, roi)
bath_crop <- crop(bathymetry, roi)
#####################################
#####################################

# create land/water raster from bathymetry data
landwater <- bathymetry
landwater[landwater > 0] <- 1
landwater[landwater <= 0] <- NA

## check raster with plot
terra::plot(landwater, col = c("#2e8b57", "#add8e6"))

# run fetchr on entire bathymetry (takes ~4.5 hours)
start <- Sys.time()
fetch_LCI <- get_fetch(
  r = landwater,     # binary land water raster
  max_dist    = 200000,        # maximum distance to calculate fetch in meters (200km)
  in_parallel = FALSE
)
print(Sys.time() - start)

#plot to check entire layer
plot(fetch_LCI)

# constrain bathymetry to roi
fetch_roi <- terra::crop(fetch_LCI, roi)

# plot new raster
plot(fetch_roi, col = viridis(nrow(fetch_roi)))

# save fetch raster
# terra::writeRaster(fetch_roi, filename = file.path(fetch_dir, "fetch.grd"), overwrite = T)


####### TEMPORARY PLACEHOLDER TO CONSTRAIN CURRENT FETCH LAYER TO ROI
# load data
fetch <- terra::rast("data/b_intermediate_data/fetch/fetch.grd")
fetch <- fetch/1000
plot(fetch$lyr.1)

fetch_mask <- mask(fetch, roi)
fetch_crop <- crop(fetch, roi)
plot(fetch_crop)
