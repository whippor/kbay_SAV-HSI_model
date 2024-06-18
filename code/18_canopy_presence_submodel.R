#####################################
### 18. Canopy Presence submodel ####
#####################################

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
dir.create("data/c_submodel_data/canopy_presence_HSI")
submodel_dir <- "data/c_submodel_data/canopy_presence_HSI"

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
presence <- terra::rast("data/b_intermediate_data/canopy_presence/presence.tif")

#####################################
#####################################

# Create presence HSI model
# mask presence to the roi
subs_mask <- mask(presence, roi)

# plot to check for correct values
# plot(subs_mask, col = viridis(nrow(subs_mask)))

# extract all values from subs roi and give 50/50 chance of kelp in unknown areas
vals1 <- data.frame(values(subs_mask))
index_vals <- vals1 %>%
  mutate(sumkelp = rowSums(select(.,c("kelp", "kelp.1", "kelp.2")), na.rm = TRUE)) %>%
  mutate(presence = case_when(sumkelp == 1 ~ 1,
                              sumkelp == 2 ~ 1,
                              sumkelp == 3 ~ 1,
                              sumkelp == 0 ~ 0.5)) %>%
  select(presence)

# join HSI values with raster
subs_mask[["HSI_value"]] <- index_vals

# subset to just HSI values
hsi_mask <- subs_mask[["HSI_value"]]

# check plot
# check plot
plot(hsi_mask, col = viridis(nrow(hsi_mask), begin = 0.3))

#####################################
#####################################

# Export data
## Suitability
terra::writeRaster(hsi_mask, 
                   filename = file.path(submodel_dir, "presenceHSI.grd"), 
                   overwrite = T)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate



