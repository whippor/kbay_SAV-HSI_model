#######################
### 001. Bathymetry ###
#######################

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
data_dir <- "data/a_raw_data/bathymetry"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### bathymetry directory
dir.create(paste0(intermediate_dir, "/",
                  "roi"))

#### bathymetry directory
dir.create(paste0(intermediate_dir, "/",
                  "bathymetry"))

bathymetry_dir <- "data/b_intermediate_data/bathymetry"

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
bathymetry <- terra::rast("data/a_raw_data/KBL-bathymetry_GWA-area_50m_EPSG3338.tiff")

# inspect the data
## coordinate reference system
terra::crs(bathymetry) # EPSG:4269

cat(crs(bathymetry))

## resolution

terra::res(bathymetry) # 50 50

### set other aspects of the data
#### units
units(bathymetry) <- "meters"

### reinspect data
cat(crs(bathymetry))
bathymetry

#####################################
#####################################

# constrain bathymetry to roi
bathy_roi <- terra::crop(bathymetry, roi)

# plot new raster
plot(bathy_roi, col = viridis(nrow(bathy_roi), begin = 0.3))

#####################################
#####################################

# export raster file
terra::writeRaster(bathy_roi, filename = file.path(bathymetry_dir, "bathymetry.grd"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate













