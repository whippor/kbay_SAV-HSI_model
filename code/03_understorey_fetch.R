#############################
### 03. Understorey Fetch ###
#############################

######### NOT WORKING

# clear environment
rm(list=setdiff(ls(), c("all_begin", "master_begin")))

# calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools,
               tidyverse,
               terra, # is replacing the raster package
               viridis)
install_version("rgdal", version = "1.6.7", repos = "http://cran.us.r-project.org")
install_github("blasee/fetchR")

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

#### fetch directory
dir.create(paste0(intermediate_dir, "/",
                  "understorey_fetch"))
fetch_dir <- "data/b_intermediate_data/understorey_fetch"

#### roi directory
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

units(bathymetry) <- "meters"


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
fetch <- fetchR::get_fetch(
  r           = landwater,     # binary land water raster
  max_dist    = 200000,        # maximum distance to calculate fetch in meters (200km)
  in_parallel = FALSE
)
print(Sys.time() - start)

#plot to check entire layer
plot(fetch)

# constrain bathymetry to roi
fetch_roi <- terra::crop(fetch, roi)

# plot new raster
plot(fetch_roi, col = viridis(nrow(fetch_roi)))

# save fetch raster
terra::writeRaster(fetch_roi, filename = file.path(fetch_dir, "fetch.grd"), overwrite = T)




