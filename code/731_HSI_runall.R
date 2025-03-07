########################
### 731. HSI run all ###
########################

# reproducibility
# renv::snapshot()
renv::restore()

# clear environment
rm(list = ls())

# calculate start time of code (determine how long it takes to complete all code)
master_begin <- Sys.time()

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
source("code/006_seagrass_presence.R")
source("code/131_understorey_HSI_runall.R")
source("code/231_canopy_HSI_runall.R")
source("code/331_seagrass_HSI_runall.R")
source("code/721_combined_HSI_model.R")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - master_begin) # print how long it takes to calculate



