##################
### SCRATCHPAD ###
##################


# load data
bathymetry <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")



# load data
substrate <- terra::rast("data/b_intermediate_data/understorey_substrate/substrate.tif")





# constrain bathymetry to 3:-30 m
bathdeep <- bathymetry
bathdeep[bathdeep > -30] <- NA

# constrain bathymetry to 3:-15 m
bathshallow <- bathymetry
bathshallow[bathshallow > -15] <- NA


#mask out each depth range on substrate maps
subsall <- substrate

subsdeep <- mask(subsall, anyNA(bathdeep), maskvalue = TRUE)
plot(subsdeep)
subsshallow <- mask(subsall, anyNA(bathshallow), maskvalue = TRUE)
plot(subsshallow)
