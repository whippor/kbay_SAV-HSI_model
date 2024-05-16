#####################
### 2. Substrate  ###
#####################

# clear environment
rm(list = ls())

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
data_dir <- "data/a_raw_data/substrate"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### bathymetry directory
dir.create(paste0(intermediate_dir, "/",
                  "substrate"))

substrate_dir <- "data/b_intermediate_data/substrate"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# import bathymetry as base raster
bathymetry <- terra::rast("~/git/kbay_SAV-HSI_model/data/b_intermediate_data/bathymetry/bathymetry.grd")
bathymetry[] <- NA


#####################################
#####################################

# load data
substrate <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/Kachemak_Subtidal_Benthic_Habitats.SHP/Kachemak_Subtidal_Benthic_Habitats.shp")

## reproject into Alaska Albers
subs_albers <- project(substrate, crs)

## make polygons into raster
subs_rast <- rasterize(subs_albers, bathymetry, "Class")

# inspect the data
## coordinate reference system
terra::crs(subs_rast) # EPSG:4269 - UTM 5N
cat(crs(subs_rast))

## resolution
terra::res(subs_rast) # 50 50

#####################################
#####################################

# constrain bathymetry to roi
bathy_roi <- terra::crop(bathymetry, roi)

# plot new raster
plot(bathy_roi)

#####################################
#####################################

# export raster file
terra::writeRaster(bathy_roi, filename = file.path(bathymetry_dir, "bathymetry.grd"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate













