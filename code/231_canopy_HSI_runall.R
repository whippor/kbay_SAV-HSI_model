################################
### 19x1. Canopy HSI run all ###
################################

# RUNS ALL BUT FETCH

# reproducibility
renv::restore()

# clear environment
rm(list = setdiff(ls(), "master_begin"))

# calculate start time of code (determine how long it takes to complete all code)
all_begin <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, # is replacing the raster package
               viridis)

#####################################
#####################################

# run scripts
source("code/11_canopy_bathymetry.R")
source("code/12_canopy_substrate.R")
# source("code/13_canopy_fetch.R")
source("code/14_canopy_presence.R")
source("code/14_x1_canopy_velocity.R")
source("code/15_canopy_bathymetry_submodel.R")
source("code/16_canopy_substrate_submodel.R")
source("code/17_canopy_fetch_submodel.R")
source("code/18_canopy_presence_submodel.R")
source("code/18_x1_canopy_velocity_submodel.R")
source("code/19_canopy_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero, 
     col = viridis(nrow(final_zero)),
     main = "Canopy Kelp HSI")


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



