##################
### SCRATCHPAD ###
##################


# load data
bathymetry <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")



# load data
substrate <- terra::rast("data/b_intermediate_data/understorey_substrate/substrate.tif")





# constrain bathymetry to 3:-30 m
bathdeep <- bathymetry
bathdeep[bathdeep > -30] <- NA

# constrain bathymetry to 3:-15 m
bathshallow <- bathymetry
bathshallow[bathshallow > -15] <- NA


#mask out each depth range on substrate maps
subsall <- substrate

subsdeep <- mask(subsall, anyNA(bathdeep), maskvalue = TRUE)
plot(subsdeep)
subsshallow <- mask(subsall, anyNA(bathshallow), maskvalue = TRUE)
plot(subsshallow)















########################
### xx. Canopy Fetch ###
########################

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
               viridis,
               fetchr)

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
                  "understorey_bathymetry"))

bathymetry_dir <- "data/b_intermediate_data/understorey_bathymetry"

roi_dir <- "data/b_intermediate_data/roi"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# export roi
terra::writeVector(roi, filename = file.path(roi_dir, "roi.shp"), overwrite = T)


#####################################
#####################################

# load data
bathymetry <- terra::rast("~/git/kbay_SAV-HSI_model/data/a_raw_data/KBL-bathymetry_GWA-area_50m_EPSG3338.tiff")

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
terra::plot(landwater, col = c("#2e8b57", "#add8e6"))

memory.size()
fetch <- fetchr::get_fetch(
  r           = landwater,     # binary land water raster
  max_dist    = 100000,        # maximum distance to calculate fetch in meters (200km)
  in_parallel = FALSE,          # run calculations in parallel
  verbose     = TRUE
)







# constrain bathymetry to roi
bathy_roi <- terra::crop(bathymetry, roi)

# plot new raster
plot(bathy_roi, col = viridis(nrow(bathy_roi), begin = 0.3))
