#############################
### 14_x1 Canopy Velocity ###
#############################

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
### input directories
data_dir <- "data/a_raw_data/CookInlet_current"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### velocity directory
dir.create(paste0(intermediate_dir, "/",
                  "canopy_velocity"))

velocity_dir <- "data/b_intermediate_data/canopy_velocity"

roi_dir <- "data/b_intermediate_data/roi"


#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# export roi
terra::writeVector(roi, filename = file.path(roi_dir, "roi.shp"), overwrite = T)


#####################################
#####################################

# load data
velocity <- terra::rast("data/a_raw_data/CookInlet_current/maxSpeed_CIOFS500.tif")

# inspect the data
## coordinate reference system
terra::crs(velocity) # EPSG:3338

cat(crs(velocity))

## resolution

terra::res(velocity) # 500 500

### set other aspects of the data
#### units
units(velocity) <- "meters"

### reinspect data
cat(crs(velocity))
velocity
terra::plot(velocity)

#####################################
#####################################

# resample to 50 x 50 m resolution
r2 <- velocity
res(r2) <- 50
velo_roi <- resample(velocity, r2)
terra::res(velo_roi)

# constrain velocity to roi
velo_roi <- terra::crop(velo_roi, roi)

# plot new raster
plot(velo_roi, col = viridis(nrow(velo_roi), begin = 0.3))

#####################################
#####################################

# export raster file
terra::writeRaster(velo_roi, filename = file.path(velocity_dir, "velocity.grd"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate


