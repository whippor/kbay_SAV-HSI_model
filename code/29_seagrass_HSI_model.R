#################################
### 29. Seagrass HSI submodel ###
#################################

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
dir.create("data/d_suitability_data/seagrass_HSI")
submodel_dir <- "data/d_suitability_data/seagrass_HSI"

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
substrate <- terra::rast("data/c_submodel_data/seagrass_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/seagrass_bathymetry_HSI/bathymetryHSI.grd")

presence <- terra::rast("data/c_submodel_data/seagrass_presence_HSI/presenceHSI.grd")

#####################################
#####################################

# create mean-value HSI for understorey (bathymetry, substrate)
under_mean <- mean(substrate[["HSI_value"]], bathymetry[["HSI_value"]])

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
near_zero <- under_zero[["HSI_value"]] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]]

# add presence as automatic "1" to raster
presencemask <- presence < 0.99
presence1 <- mask(presence, presencemask, maskvalue = 1)
final_zero <- merge(presence1, near_zero)

plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero,
     col = viridis(nrow(final_zero)),
     main = "Seagrass HSI")

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(final_zero, 
                   filename = file.path(submodel_dir, "seagrass_HSI.tif"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate


 
r <- rast(ncols=5, nrows=5, xmin=0, xmax=5, ymin=0, ymax=5)
r[] <- 1:25
r[1,] <- 5
r[,2] <- 10
r[r>10] <- NA
