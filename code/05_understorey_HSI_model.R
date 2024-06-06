####################################
### 5. Understorey HSI submodel ####
####################################

# clear environment
rm(list=setdiff(ls(), "all_begin"))

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
### roi directory
roi_dir <- "data/b_intermediate_data/roi"

### output directories
#### submodel directory
dir.create("data/d_suitability_data/understorey_HSI")
submodel_dir <- "data/d_suitability_data/understorey_HSI"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect(roi_dir)

#####################################
#####################################

# load data
substrate <- terra::rast("data/c_submodel_data/understorey_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/understorey_bathymetry_HSI/bathymetryHSI.grd")

#####################################
#####################################

# create mean-value HSI for understorey (bathymetry, substrate)
under_mean <- mean(substrate[["HSI_value"]], bathymetry[["HSI_value"]])

# create HSI setting any zero to zero 
under_zero <- mean(substrate[["HSI_value"]], bathymetry[["HSI_value"]])
sub_zero <- substrate
subst(sub_zero[["HSI_value"]], 0.01:1, 1)
bath_zero <- bathymetry
subst(bath_zero[["HSI_value"]], 0.01:1, 1)
final_zero <- under_zero[["HSI_value"]] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]]

plot(final_zero, col = viridis(nrow(final_zero), begin = 0.3))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(final_zero, 
                   filename = file.path(submodel_dir, "understorey_HSI.tif"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



