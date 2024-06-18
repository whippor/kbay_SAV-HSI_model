################################
### 19. Canopy HSI submodel ####
################################

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
### roi directory
roi_dir <- "data/b_intermediate_data/roi"

### output directories
#### submodel directory
dir.create("data/d_suitability_data/canopy_HSI")
submodel_dir <- "data/d_suitability_data/canopy_HSI"

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
substrate <- terra::rast("data/c_submodel_data/canopy_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/canopy_bathymetry_HSI/bathymetryHSI.grd")

fetch <- terra::rast("data/c_submodel_data/canopy_fetch_HSI/fetchHSI.grd")

presence <- terra::rast("data/c_submodel_data/canopy_presence_HSI/presenceHSI.grd")

#####################################
#####################################

# create mean-value HSI for canopy (bathymetry, substrate, fetch)
under_mean <- mean(substrate[["HSI_value"]], 
                   bathymetry[["HSI_value"]],
                   fetch[["HSI_value"]])

# weighted mean HSI
# under_mean <- c(substrate[["HSI_value"]], bathymetry[["HSI_value"]])
# under_mean <- terra::weighted.mean(under_mean, c(2, 1))
# names(under_mean) <- "HSI_value"

# create HSI setting any zero to zero 
under_zero <- under_mean
sub_zero <- substrate
sub_zero[sub_zero < 0.01] <- 0
sub_zero[sub_zero > 0.01] <- 1
bath_zero <- bathymetry
bath_zero[bath_zero < 0.01] <- 0
bath_zero[bath_zero > 0.01] <- 1
fetch_zero <- fetch
fetch_zero[fetch_zero < 0.01] <- 0
fetch_zero[fetch_zero > 0.01] <- 1
near_zero <- under_zero[["HSI_value"]] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]] *
  fetch_zero[["HSI_value"]]

# add presence as automatic "1" to raster
presencemask <- presence < 0.99
presence1 <- mask(presence, presencemask, maskvalue = 1)
final_zero <- merge(presence1, near_zero)

plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero,
     main = "Canopy Kelp HSI")

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(final_zero, 
                   filename = file.path(submodel_dir, "canopy_HSI.tif"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate


neardf <- data.frame(near_zero)
notone <- neardf %>%
  filter(HSI_value  == 1)
finaldf <- data.frame(final_zero)
allone <- finaldf %>%
  filter(HSI_value == 1)
