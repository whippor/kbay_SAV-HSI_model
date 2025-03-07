####################################
### 131. Understorey HSI run all ###
####################################

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
source("code/111_understorey_bathymetry_submodel.R")
source("code/112_understorey_substrate_submodel.R")
source("code/113_understorey_fetch_submodel.R")
# source("code/114_understorey_velocity_submodel.R")
source("code/121_understorey_HSI_model.R")


plot(final_zero, col = viridis(nrow(final_zero), begin = 0.3))

plet(final_zero, 
     col = viridis(nrow(final_zero)),
     main = "Understorey Kelp HSI")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - all_begin) # print how long it takes to calculate



