#############################
### 221. Canopy HSI model ###
#############################

## NOTE: Li et al 2024 found that minimum SST and annual SST range
## were strongest predictors of Eualaria presence

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

### output directories
#### submodel directory
dir.create("data/d_suitability_data/canopy_HSI")
submodel_dir <- "data/d_suitability_data/canopy_HSI"

#####################################
#####################################

# load data
substrate <- terra::rast("data/c_submodel_data/canopy_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/canopy_bathymetry_HSI/bathymetryHSI.grd")

fetch <- terra::rast("data/c_submodel_data/canopy_fetch_HSI/fetchHSI.grd")

presence <- terra::rast("data/c_submodel_data/canopy_presence_HSI/presenceHSI.grd")

# velocity <- terra::rast("data/c_submodel_data/canopy_velocity_HSI/velocityHSI.grd")
# match extent of velocity raster to the others
# ext(velocity) <- ext(bathymetry)

#####################################
#####################################

# create mean-value HSI for canopy (bathymetry, substrate, fetch)
under_mean <- mean(substrate[["HSI_value"]], 
                   bathymetry[["HSI_value"]],
                   fetch[["HSI_value"]]
                   # , velocity[["HSI_value"]]
)

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
#velo_zero <- velocity
#velo_zero[velo_zero < 0.01] <- 0
#velo_zero[velo_zero > 0.01] <- 1
near_zero <- under_zero[["HSI_value"]] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]] *
  fetch_zero[["HSI_value"]] 
# * velo_zero[["HSI_value"]]


# add presence as automatic "1" to raster
presencemask <- presence < 0.99
presence1 <- mask(presence, presencemask, maskvalue = 1)
final_zero <- merge(presence1, near_zero)

plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero, 
     col = viridis(nrow(final_zero)),
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


