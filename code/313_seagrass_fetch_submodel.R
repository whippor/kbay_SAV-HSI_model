####################################
### 313. Seagrass Fetch submodel ###
####################################

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
source("code/000_function_interpolate_y.R")

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### roi directory
roi_dir <- "data/b_intermediate_data/roi"

### output directories
#### submodel directory
dir.create("data/c_submodel_data/seagrass_fetch_HSI")
submodel_dir <- "data/c_submodel_data/seagrass_fetch_HSI"

#####################################
#####################################

# set parameters

# define vector for region of interest
roi <- terra::vect(roi_dir)

#####################################
#####################################

# load data
fetch <- terra::rast("data/b_intermediate_data/fetch/fetch.grd")
fetch <- crop(fetch, roi)
fetch <- fetch/1000

tam_fetch <- read_csv("data/x_tam_tables/seagrass/seagrass_fetch.csv")
tam_fetch <- tam_fetch %>%
  arrange(fetch.km)

#####################################
#####################################

# Create bathymetry HSI model
# mask bathymetry to the roi
fetch_mask <- mask(fetch, roi)

# extract all values from bath_roi
vals1 <- data.frame(values(fetch_mask))

# FUNCTION TO FIND Y FOR ANY GIVEN X WITH IMPORTED SLOPES

## Calculate slopes
slopes <- diff(tam_fetch$fetch.km.SIV) / diff(tam_fetch$fetch.km)

# calculate index from raster values with function
index_vals <- interpolate_y(vals1$lyr.1, tam_fetch)

# join HSI values with raster
fetch_mask[["HSI_value"]] <- index_vals

# check plot
plot(fetch_mask, col = viridis(nrow(fetch_mask)))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(fetch_mask, 
                   filename = file.path(submodel_dir, "fetchHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



