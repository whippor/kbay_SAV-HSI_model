#################################
### 230. Seagrass HSI run all ###
#################################

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
source("code/21_seagrass_bathymetry.R")
source("code/22_seagrass_substrate.R")
# source("code/23_seagrass_fetch.R")
source("code/24_seagrass_presence.R")
source("code/25_seagrass_bathymetry_submodel.R")
source("code/26_seagrass_substrate_submodel.R")
source("code/27_seagrass_fetch_submodel.R")
source("code/28_seagrass_presence_submodel.R")
source("code/29_seagrass_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero,
     col = viridis(nrow(final_zero)),
     main = "Seagrass HSI")


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



