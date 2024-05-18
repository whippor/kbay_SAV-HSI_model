##########################################
### 3. Understorey Bathymetry submodel ###
##########################################

# clear environment
rm(list = ls())

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
submodel_dir <- "data/c_submodel_data"

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
bathymetry <- terra::rast("data/b_intermediate_data/bathymetry/bathymetry.grd")

tam_bath <- read_csv("data/x_tam_tables/understorey/understorey_depth.csv")
tam_bath <- tam_bath %>%
  arrange(depth.m)
xvalue <- names(tam_bath[1])
yvalue <- names(tam_bath[2])

#####################################
#####################################

# Create bathymetry HSI model
# mask bathymetry to the roi
bath_mask <- mask(bathymetry, roi)

# extract all values from bath_roi
vals1 <- data.frame(values(bath_mask))




#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate


# FUNCTION TO FIND Y FOR ANY GIVEN X WITH IMPORTED SLOPES

# Calculate slopes
slopes <- diff(tam_bath$depth.m.SIV) / diff(tam_bath$depth.m)

# Function to interpolate y value for given x values
interpolate_y <- function(x_values) {
  # Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    # Check if x_value is outside the range of tam_bath$x
    if (x_value < min(tam_bath$depth.m) || x_value > max(tam_bath$depth.m)) {
      interpolated_y_values[i] <- 0  # Set interpolated y value to zero
    } else {
      # Find the interval where x_value lies
      idx <- findInterval(x_value, tam_bath$depth.m)
      
      # Linear interpolation
      y1 <- tam_bath$depth.m.SIV[idx]
      y2 <- tam_bath$depth.m.SIV[idx + 1]
      slope <- slopes[idx]
      interpolated_y_values[i] <- y1 + slope * (x_value - tam_bath$depth.m[idx])
    }
  }
  
  return(interpolated_y_values)
}

# Example usage:
user_x_values <- c(2, 5, 0, -1, -6, -6.9, -8)  # User-supplied x values
user_y_values <- interpolate_y_values(user_x_values)
print(data.frame(x = user_x_values, y = user_y_values))
