#######################
### 999. SCRATCHPAD ###
#######################


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







# spring Seldovia 5 and -1
# neap Seldovia 3.5 and 1.5

library(tidyverse)
library(terra)

crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("C:/Users/Ross.Whippo/Desktop/jakolof_polygon.kml")
roi <- project(roi, crs)

bathymetry <- terra::rast("C:/Users/Ross.Whippo/Desktop/Kachemak_bathymetry_8m.tif/Kachemak_bathymetry_8m.tif")
bathymetry <- project(bathymetry, crs)
bathymetry <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")

plot(bathymetry)
bath_mask <- terra::mask(bathymetry, roi)
plot(bath_mask)
bath_clip <- terra::crop(bath_mask, roi)
plot(bath_clip)

# 5 to -1
springhigh <- bath_mask
springlow <- bath_mask
springhigh[springhigh > 5] <- NA
springhighdf <- data.frame(springhigh)


springlow[springlow > -1] <- NA
springlowdf <- data.frame(springlow)
plot(springlow)

8.003082*8.003082 # 64.04932
#springhigh
64.04932*245509 # 15724685
#springlow
64.04932*237557 #15215364

# 3.141259r^2 = 15724685
15724685/3.14159 # 5005327
sqrt(5005327) # 2237.259 m radius circle (high)

15215364/3.14159 # 4843205
sqrt(4843205) # 2200.728

# truncated cone volume
# V = (1/3) * pi * h (r^2 + r * R + R^2)
# V = (1/3) * 3.14159 * 6 (2200.728^2 + 2200.728 * 2237.259 + 2237.728^2)
3.14159*6 # 18.84954
2200.728^2 # 4843204
2237.728^2 #5007427
2200.728*2237.259 # 4923599
4843204+5007427+4923599 # 14774230
14774230*18.84954 # 278487439
278487439*(1/3) # 92829146
# ~ 93 million cubic meters of water on a spring tide


sum(abs(springhighdf)) * 64.04932 # 776373227
sum(abs(springlowdf)) * 64.04932 # 776215528
776373227 - 776215528




# 3.5 to 1.5
neaphigh <- bath_mask
neaplow <- bath_mask
neaphigh[neaphigh > 3.5] <- NA
neaphighdf <- data.frame(neaphigh)

neaplow[neaplow > 1.5] <- NA
neaplowdf <- data.frame(neaplow)

#neaphigh
64.04932*245509 # 15724685
#neaplow
60.04932*245509 # 14742649

# 3.141259r^2 = 15724685










# UA STUFF


interpolate_y(vals_bath, under_tam_bath)





## Function to interpolate y value for given x values
interpolate_y <- function(x_values, tam_table) {
  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_table$x
    ifelse(x_value < min(as.numeric(unlist(tam_table[,1]))), 0, x_value)  
    ifelse(x_value > max(as.numeric(unlist(tam_table[,1]))), 0, x_value)  
    ifelse(is.na(x_value) == TRUE, 0, x_value) 
    # Set interpolated y value to zero
    
    ##### Find the interval where x_value lies
    idx <- findInterval(x_value, as.numeric(unlist(tam_table[,1])))
    
    ##### Linear interpolation
    y1 <- as.numeric(unlist(tam_table[,2]))[idx]
    y2 <- as.numeric(unlist(tam_table[,2]))[idx + 1]
    slope <- slopes[idx]
    interpolated_y_values[i] <- y1 + slope * (x_value - as.numeric(unlist(tam_table[,1]))[idx])
  }
  return(interpolated_y_values)
}



for (j in seq_along(soboldf[1])) {
  x_value <- soboldf[j, 1]
  
  
}





colnames(x)[z]
colnames(x[[1]])[1]

