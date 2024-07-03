#######################
### 32. HSI run all ###
#######################

# reproducibility
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
source("code/00_directory_creation.R")
source("code/08_understorey_HSI_runall.R")
source("code/19_x1_canopy_HSI_runall.R")
source("code/29_x1_seagrass_HSI_runall.R")
source("code/31_combined_HSI_model.R")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - master_begin) # print how long it takes to calculate



