########################################
### 28_x1 Seagrass Velocity submodel ###
########################################

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
dir.create("data/c_submodel_data/seagrass_velocity_HSI")
submodel_dir <- "data/c_submodel_data/seagrass_velocity_HSI"

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
velocity <- terra::rast("data/b_intermediate_data/seagrass_velocity/velocity.grd")

tam_velo <- read_csv("data/x_tam_tables/seagrass/seagrass_velocity.csv")
tam_velo <- tam_velo %>%
  arrange(velocity.ms)

#####################################
#####################################

# Create velocity HSI model
# mask velocity to the roi
velo_mask <- mask(velocity, roi)

# extract all values from bath_roi
vals1 <- data.frame(values(velo_mask))

# FUNCTION TO FIND Y FOR ANY GIVEN X WITH IMPORTED SLOPES

## Calculate slopes
slopes <- diff(tam_velo$velocity.ms.SIV) / diff(tam_velo$velocity.ms)

## Function to interpolate y value for given x values
interpolate_y <- function(x_values) {
  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_velo$x
    if (x_value < min(tam_velo$velocity.ms) || 
        x_value > max(tam_velo$velocity.ms) || 
        is.na(x_value) == TRUE ) {
      interpolated_y_values[i] <- 0  # Set interpolated y value to zero
    } else {
      ##### Find the interval where x_value lies
      idx <- findInterval(x_value, tam_velo$velocity.ms)
      
      ##### Linear interpolation
      y1 <- tam_velo$velocity.ms.SIV[idx]
      y2 <- tam_velo$velocity.ms.SIV[idx + 1]
      slope <- slopes[idx]
      interpolated_y_values[i] <- y1 + slope * (x_value - tam_velo$velocity.ms[idx])
    }
  }
  
  return(interpolated_y_values)
}

# calculate index from raster values with function
index_vals <- interpolate_y(vals1$maxSpeed_CIOFS500)

# join HSI values with raster
velo_mask[["HSI_value"]] <- index_vals

# check plot
plot(velo_mask, col = viridis(nrow(velo_mask), begin = 0.3))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(velo_mask, 
                   filename = file.path(submodel_dir, "velocityHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



