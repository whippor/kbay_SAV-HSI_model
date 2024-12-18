###############################
### 23. Seagrass Presence   ###
###############################

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
data_dir <- "data/a_raw_data/a_raw_data/Seagrass_presence/"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### substrate directory
dir.create(paste0(intermediate_dir, "/",
                  "seagrass_presence"))

seagrass_dir <- "data/b_intermediate_data/seagrass_presence"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# import bathymetry as base raster
bathymetry <- terra::rast("data/b_intermediate_data/canopy_bathymetry/bathymetry.grd")

#####################################
#####################################

# load data
seagrass <- terra::vect("data/a_raw_data/Seagrass_presence/KBAY_seagrass.shp")

## reproject into Alaska Albers
seagrass_albers <- project(seagrass, crs)

## make polygons into raster
seagrass_rast <- rasterize(seagrass_albers, bathymetry)

# inspect the data
## coordinate reference system
cat(crs(seagrass_rast)) # EPSG:3338


## resolution
terra::res(seagrass_rast) # 50 50

# check plot
plot(seagrass_rast)

#####################################
#####################################

# expand substrate raster to bathymetry layer
seagrass_expand <- terra::extend(seagrass_rast, bathymetry)
varnames(seagrass_expand) <- "seagrass"

# fill empty space in raster with Unclassified
unclass_rast <- rast(ncol = 1064,
                     nrow = 994,
                     xmin = 119250,
                     xmax = 172450,
                     ymin = 1046527, 
                     ymax = 1096227,
                     crs = crs(seagrass_expand))

seagrass_final <- terra::merge(seagrass_expand, unclass_rast)
names(seagrass_final) <- "seagrass"

# plot new raster
plot(seagrass_final, col = viridis(nrow(seagrass_final), begin = 0.2, end = 0.8))

#####################################
#####################################

# export raster file
terra::writeRaster(seagrass_final, filename = file.path(seagrass_dir, "presence.tif"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate

