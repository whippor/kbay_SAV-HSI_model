###############################
### 231. Canopy HSI run all ###
###############################

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
# source("code/004_velocity.R")
source("code/005_canopy_presence.R")
source("code/211_canopy_bathymetry_submodel.R")
source("code/212_canopy_substrate_submodel.R")
source("code/213_canopy_fetch_submodel.R")
# source("code/214_canopy_velocity_submodel.R")
source("code/215_canopy_presence_submodel.R")
source("code/221_canopy_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero)))

plet(final_zero, 
     col = viridis(nrow(final_zero)),
     main = "Canopy Kelp HSI")


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



