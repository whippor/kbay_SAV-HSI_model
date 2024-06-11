#############################
### 13. Canopy Presence   ###
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
data_dir <- "data/a_raw_data/a_raw_data/Kelp_presence/"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### substrate directory
dir.create(paste0(intermediate_dir, "/",
                  "canopy_presence"))

kelp_dir <- "data/b_intermediate_data/canopy_presence"

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
bathymetry <- terra::rast("~/git/kbay_SAV-HSI_model/data/b_intermediate_data/canopy_bathymetry/bathymetry.grd")

#####################################
#####################################

# load data
kelp2000 <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/Kelp_presence/kelp2000_final.shp")
kelp2001 <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/Kelp_presence/kelp2001_final.shp")
kelp2002 <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/Kelp_presence/kelp2002_final.shp")

## reproject into Alaska Albers
k00_albers <- project(kelp2000, crs)
k01_albers <- project(kelp2001, crs)
k02_albers <- project(kelp2002, crs)

## make polygons into raster
k00_rast <- rasterize(k00_albers, bathymetry)
k01_rast <- rasterize(k01_albers, bathymetry)
k02_rast <- rasterize(k02_albers, bathymetry)

kelp <- c(k00_rast,
          k01_rast,
          k02_rast)

# inspect the data
## coordinate reference system
cat(crs(kelp)) # EPSG:3338


## resolution
terra::res(kelp) # 50 50

# check plot
plot(kelp)

#####################################
#####################################

# expand substrate raster to bathymetry layer
kelp_expand <- terra::extend(kelp, bathymetry)
varnames(kelp_expand) <- "kelp"

# fill empty space in raster with Unclassified
unclass_rast <- rast(ncol = 1064,
                     nrow = 994,
                     xmin = 119250,
                     xmax = 172450,
                     ymin = 1046527, 
                     ymax = 1096227,
                     crs = crs(kelp_expand))

kelp_final <- terra::merge(kelp_expand, unclass_rast)
names(kelp_final) <- c("kelp", "kelp", "kelp")

# plot new raster
plot(kelp_final, col = viridis(nrow(kelp_final), begin = 0.2, end = 0.8))

#####################################
#####################################

# export raster file
terra::writeRaster(kelp_final, filename = file.path(kelp_dir, "presence.tif"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate

