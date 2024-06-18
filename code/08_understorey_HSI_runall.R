###################################
### 08. Understorey HSI run all ###
###################################

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
source("code/01_understorey_bathymetry.R")
source("code/02_understorey_substrate.R")
# source("code/03_understorey_fetch.R")
source("code/04_understorey_bathymetry_submodel.R")
source("code/05_understorey_substrate_submodel.R")
source("code/06_understorey_fetch_submodel.R")
source("code/07_understorey_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero), begin = 0.3))


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



