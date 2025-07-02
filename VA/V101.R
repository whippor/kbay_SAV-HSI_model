####################################################
### V101. Validation Sampling - Virtual - Canopy ###
####################################################

# reproducibility
renv::restore()

# clear environment
rm(list = ls())

# calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, 
               tidyterra,
               viridis,
               leaflet,
               ggrepel)

################################################
### Canopy HSI model (from 221, no presence) ###
################################################


# set directories

### output directories
#### submodel directory
dir.create("data/e_validation_data/canopy_HSI")
submodel_dir <- "data/e_validation_data/canopy_HSI"

#####################################
#####################################

MLLWcor <- 1.546 # determined w/ https://vdatum.noaa.gov/vdatumweb/

# load data
substrate <- terra::rast("data/c_submodel_data/canopy_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/canopy_bathymetry_HSI/bathymetryHSI.grd")

fetch <- terra::rast("data/c_submodel_data/canopy_fetch_HSI/fetchHSI.grd")

presence <- terra::rast("data/c_submodel_data/canopy_presence_HSI/presenceHSI.grd")

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# create bathmask
# mask out all HSI below -30 m, above 0 m (w/  - 1.546 m correction for NAVD88/MLLW
bathmask <- clamp(bathymetry, lower = -30 - MLLWcor, upper = 0 - MLLWcor, values = FALSE)


#####################################
#####################################

# create geometric mean-value HSI for canopy (bathymetry, substrate, fetch)
under_mean <- exp(mean(log(c((substrate[["HSI_value"]] + 1), 
                             (bathymetry[["HSI_value"]] + 1),
                             (fetch[["HSI_value"]] + 1)))))

# create HSI setting any zero to zero 
# under_zero <- under_mean # for arithmetic
under_zero <- (under_mean - 1) # for geometric
sub_zero <- substrate
sub_zero[sub_zero < 0.01] <- 0
sub_zero[sub_zero > 0.01] <- 1
bath_zero <- bathymetry
bath_zero[bath_zero < 0.01] <- 0
bath_zero[bath_zero > 0.01] <- 1
fetch_zero <- fetch
fetch_zero[fetch_zero < 0.01] <- 0
fetch_zero[fetch_zero > 0.01] <- 1
near_zero <- under_zero[['mean']] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]] *
  fetch_zero[["HSI_value"]] 

# determine overlap of historic presence and model predictions
# crop presence to roi
presence_roi <- terra::mask(presence, roi)

# change 0.5 values to 0
presence_roi[presence_roi < 1] <- 0

# prediction plot
plot(near_zero)

# historic plot
plot(presence_roi)

# extract observed values
obs_values <- data.frame(values(presence_roi))
# extract predicted values
pred_values <- data.frame(values(near_zero))
# join them and calculate model prediction values
joint_values <- obs_values %>%
  bind_cols(pred_values) %>%
  rename(prediction = `mean`) %>%
  mutate(model_performance = case_when(HSI_value < prediction ~ HSI_value + prediction,
                                       HSI_value == prediction ~ 0,
                                       HSI_value > prediction ~ prediction - HSI_value)) %>%
  select(model_performance)

# duplicate raster and replace values
model_diff <- near_zero
model_diff[["model_performance"]] <- joint_values
model_diff <- model_diff[["model_performance"]]

# visualize performance in plot
plot(model_diff, col = viridis(nrow(model_diff), option = "turbo"))

# with bathmask
model_diff30 <- mask(model_diff, bathmask)
model_diff30 <- model_diff30[["lyr1"]]
plot(model_diff30, col = viridis(nrow(model_diff), option = "turbo"))

# histogram of performance values with 30 m cutoff
 canopy_perf_df <- data.frame(model_diff30)
 canopy_perf_df %>%
ggplot(aes(x=lyr1)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw()
 
# binomial regression excluding joint absence
obs_values %>%
   bind_cols(pred_values) %>%
   rename(prediction = `mean`) %>%
   mutate(model_performance = case_when(HSI_value < prediction ~ HSI_value + prediction,
                                        HSI_value == prediction ~ 0,
                                        HSI_value > prediction ~ prediction - HSI_value)) %>%
  filter(!(prediction == 0 & HSI_value == 0)) %>%
  ggplot(aes(x = prediction, y = HSI_value)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  theme_bw()


#####################################
#####################################

# Export data
## Suitability
# terra::writeRaster(final_zero, 
#                   filename = file.path(submodel_dir, "canopy_HSI.tif"), 
#                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



#####################################
#####################################

