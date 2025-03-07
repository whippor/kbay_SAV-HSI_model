#################################
### 331. Seagrass HSI run all ###
#################################

# reproducibility
# renv::snapshot()
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
source("code/001_bathymetry.R")
source("code/002_substrate.R")
# source("code/003_fetch.R")
source("code/004_velocity.R")
source("code/006_seagrass_presence.R")
source("code/311_seagrass_bathymetry_submodel.R")
source("code/312_seagrass_substrate_submodel.R")
source("code/313_seagrass_fetch_submodel.R")
# source("code/314_seagrass_velocity_submodel.R")
source("code/316_seagrass_presence_submodel.R")
source("code/321_seagrass_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero,
     col = viridis(nrow(final_zero)),
     main = "Seagrass HSI")


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



