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


for (i in 1:length(all_list))


all_list <- list(vals_bath, vals_fetch)

for (i in 1:length(all_list)) {
  
  assign(paste(i, names(all_list[[i]]), sep = "_"), all_list[[i]])

  }


















source(sobol_interpolate_y())
debug(sobol_interpolate_y)


sobol_interpolate_y <- function(x_values, 
                                tam_table, 
                                samples = 10,
                                sobol_order = "first") {
  # name parameters
  all_list <- x_values
  
  
  for (i in 1:length(all_list)) {
    assign(paste("x", i, sep = "_"), all_list[[i]])
  }
  
  
  x <<- tam_table
  
  # check that there is a tam table for each input
 
    # add adjusted SIV values to TAM tables
    x <<- map(x, \(x) rename(x, "SIV" = contains("SIV")))
    for (i in seq_along(x)) {
      x[[i]]$unif_rand <- runif(nrow(x[[i]]), min = -0.2, max = 0.2)
      x[[i]]$TAM_adj <- x[[i]]$SIV + x[[i]]$unif_rand
      x[[i]]$TAM_adj[1] <- x[[i]]$SIV[1]
      x[[i]]$TAM_adj[length(x[[i]]$TAM_adj)] <- x[[i]]$SIV[length(x[[i]]$SIV)]
      x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj > 1, 1, x[[i]]$TAM_adj)
      x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj < 0, 0, x[[i]]$TAM_adj)
    }
    
    # create sobol sequence for x values
    z <<- length(all_list)
    
    
    for (i in 1:z) {
      N <- samples
      df <- paste0("x", sep = "_", i)
      params <- paste(colnames(get(df)), "sobol", sep = "_")
      mat <- sobol_matrices(N = N,
                            params = params,
                            order = sobol_order)
      mat[, 1] <- qunif(mat[, 1], min(get(df)), max(get(df)))
      
      soboldf <<- data.frame(mat, matrix(NA, nrow = length(mat), ncol = 1))
      soboldf <<- rename(soboldf, !!paste(unlist(colnames(get(df))), "newSIV", sep = "_") := 2)
      
      for (j in 1:nrow(soboldf)) {
        x_value <- soboldf[j, 1]
        
        #### Check if x_value is outside the range of tam_table$x
        ifelse(x_value < min(as.numeric(unlist(x[[1]][1]))), 0, x_value)
        ifelse(x_value > max(as.numeric(unlist(x[[1]][1]))), 0, x_value)
        ifelse(is.na(x_value) == TRUE, 0, x_value)
        # Set interpolated y value to zero
        
        ##### Find the interval where x_value lies
        idx <- findInterval(x_value, as.numeric(unlist(x[[1]][1])))
        ##### Find slopes
        slopes <- diff(x[[1]]$TAM_adj) / diff(unlist(x[[1]][1]))
        
        ##### Linear interpolation
        y1 <- as.numeric(unlist(x[[1]]$TAM_adj))[idx]
        y2 <- as.numeric(unlist(x[[1]]$TAM_adj))[idx + 1]
        slope <- slopes[idx]
        soboldf[j, 2] <<- y1 + slope * (x_value - as.numeric(unlist(x[[1]][1]))[idx])
        
      }
      
      assign(paste(unlist(colnames(x[[i]][1])), i, "sobol_df", sep = "_"), value = soboldf)
      
    }
  
  
}
















lotka_volterra_fun <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- r * X * (1 - X / K) - alpha * X * Y
    dY <- -m * Y + theta * X * Y
    list(c(dX, dY))
  })
}

# Define the settings of the sensitivity analysis
N <- 2 ^ 5 # Sample size of sample matrix
params <- c("r", "alpha", "m", "theta", "K", "X", "Y") # Parameters

# Define the times
times <- seq(5, 20, 1)

# Define the times at which the output is wanted
timeOutput <- c(10, 15)

# Construct the sample matrix
mat <- sobol_matrices(N = N, params = params)

# Transform to appropriate distributions
mat[, "r"] <- qunif(mat[, "r"], 0.8, 1.8)
mat[, "alpha"] <- qunif(mat[, "alpha"], 0.2, 1)
mat[, "m"] <- qunif(mat[, "m"], 0.6, 1)
mat[, "theta"] <- qunif(mat[, "theta"], 0.05, 0.15)
mat[, "K"] <- qunif(mat[, "K"], 47, 53)
mat[, "X"] <- floor(mat[, "X"] * (15 - 8 + 1) + 8)
mat[, "Y"] <- floor(mat[, "Y"] * (2 - 6 + 1) + 6)

# Run the model
y <- list()
for (i in 1:nrow(mat)) {
  y[[i]] <- sobol_ode(d = mat[i, ],
                      times = times,
                      timeOutput = timeOutput,
                      state = c(X = mat[[i, "X"]], Y = mat[[i, "Y"]]),
                      func = lotka_volterra_fun)
}




for (i in 1:z) {
  df <- paste("x", i, sep = "_")
  N <- length(get(df))
  params <- paste0(names(x[[i]][1]), sep = "_", "sobol")
  mat <- sobol_matrices(N = N/3,
                        params = params,
                        order = sobol_order,
                        type = "R")
  mat[, 1] <- qunif(mat[, 1], min(get(df)), max(get(df)))
  
  soboldf <- data.frame(mat, matrix(NA, nrow = length(mat), ncol = 1))
  soboldf <- rename(soboldf, !!paste(names(x[[i]][1]), "newSIV", sep = "_") := 2)
  sobol_vec_x <- as.vector(soboldf[,1])
  sobol_vec_SIV <- c()
  min_x <- min(as.numeric(unlist(x[[i]][1])))
  max_x <- max(as.numeric(unlist(x[[i]][1])))
  ##### Find slopes
  slopes <- diff(x[[i]]$TAM_adj) / diff(unlist(x[[i]][1]))
  
  for (j in 1:length(sobol_vec_x)) {
    x_value <- sobol_vec_x[j]
    
    #### Check if x_value is outside the range of tam_table$x
    ifelse(x_value < min_x, 0, x_value)
    ifelse(x_value > max_x, 0, x_value)
    ifelse(is.na(x_value) == TRUE, 0, x_value)
    # Set interpolated y value to zero
    
    ##### Find the interval where x_value lies
    idx <- findInterval(x_value, as.numeric(unlist(x[[i]][1])))
    
    ##### Linear interpolation
    y1 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx]
    y2 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx + 1]
    slope <- slopes[idx]
    sobol_vec_SIV <- c(sobol_vec_SIV, y1 + slope * (x_value - as.numeric(unlist(x[[i]][1]))[idx]))
    
  }
  
  soboldf[,2] <- sobol_vec_SIV
  output_list[[paste(unlist(colnames(x[[i]][1])))]] <- list(soboldf)

  }



library(sensobol)
library(dplyr)
library(purrr)

x_values <- c(vals_bath, vals_fetch) # for testing function
tam_table <- list(under_tam_bath, under_tam_fetch) # for testing function
sobol_order = "first"
  
###
N = 10000
params = c("depth", "fetch")
perm = 2

MC_SIV <- c()

for (k in 1:perm) {

SIV_set <- c()

for (i in 1:length(params)) {
  
  x <- tam_table
  
  # add adjusted SIV values to TAM tables
  x <- map(x, \(x) rename(x, "SIV" = contains("SIV")))
  for (m in seq_along(x)) {
    x[[m]]$unif_rand <- runif(nrow(x[[m]]), min = -0.2, max = 0.2)
    x[[m]]$TAM_adj <- x[[m]]$SIV + x[[m]]$unif_rand
    x[[m]]$TAM_adj[1] <- x[[m]]$SIV[1]
    x[[m]]$TAM_adj[length(x[[m]]$TAM_adj)] <- x[[m]]$SIV[length(x[[m]]$SIV)]
    x[[m]]$TAM_adj <- ifelse(x[[m]]$TAM_adj > 1, 1, x[[m]]$TAM_adj)
    x[[m]]$TAM_adj <- ifelse(x[[m]]$TAM_adj < 0, 0, x[[m]]$TAM_adj)
  }

mat <- sobol_matrices(N = N, params = params, order = sobol_order)
mat[,1] <- qunif(mat[,1], -30, 3)
mat[,2] <- qunif(mat[,2], 0.01, 50)

##### Find slopes
slopes <- diff(x[[i]]$TAM_adj) / diff(unlist(x[[i]][1]))

min_x <- min(as.numeric(unlist(x[[i]][1])))
max_x <- max(as.numeric(unlist(x[[i]][1])))

sobol_vec_SIV <- c()

for (j in 1:nrow(mat)) {
  x_value <- mat[j,i]
  
  #### Check if x_value is outside the range of tam_table$x
  ifelse(x_value < min_x, 0, x_value)
  ifelse(x_value > max_x, 0, x_value)
  ifelse(is.na(x_value) == TRUE, 0, x_value)
  # Set interpolated y value to zero
  
  ##### Find the interval where x_value lies
  idx <- findInterval(x_value, as.numeric(unlist(x[[i]][1])))
  
  ##### Linear interpolation
  y1 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx]
  y2 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx + 1]
  slope <- slopes[idx]
  sobol_vec_SIV <- c(sobol_vec_SIV, y1 + slope * (x_value - as.numeric(unlist(x[[i]][1]))[idx]))
  
}

SIV_set[[paste(unlist(colnames(x[[i]][1])), "SIV", sep = "_")]] <- sobol_vec_SIV

}

MC_SIV[[paste("run", k, sep = "_")]] <- SIV_set

rm(SIV_set)

}
###







for (i in 1:ncol(mat)) {
  
  N <- length(get(df))
  params <- paste0(names(x[[i]][1]), sep = "_", "sobol")
#  mat <- sobol_matrices(N = N/3,
#                        params = params,
#                        order = sobol_order,
#                        type = "R")
#  mat[, 1] <- qunif(mat[, 1], min(get(df)), max(get(df)))
  
  soboldf <- data.frame(mat, matrix(NA, nrow = length(mat), ncol = 1))
  soboldf <- rename(soboldf, !!paste(names(x[[i]][1]), "newSIV", sep = "_") := 2)
  sobol_vec_x <- as.vector(soboldf[,1])
  sobol_vec_SIV <- c()
  min_x <- min(as.numeric(unlist(x[[i]][1])))
  max_x <- max(as.numeric(unlist(x[[i]][1])))
  ##### Find slopes
  slopes <- diff(x[[i]]$TAM_adj) / diff(unlist(x[[i]][1]))
  
  for (j in 1:length(mat)) {
    x_value <- mat[i,j]
    
    #### Check if x_value is outside the range of tam_table$x
    ifelse(x_value < min_x, 0, x_value)
    ifelse(x_value > max_x, 0, x_value)
    ifelse(is.na(x_value) == TRUE, 0, x_value)
    # Set interpolated y value to zero
    
    ##### Find the interval where x_value lies
    idx <- findInterval(x_value, as.numeric(unlist(x[[i]][1])))
    
    ##### Linear interpolation
    y1 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx]
    y2 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx + 1]
    slope <- slopes[idx]
    sobol_vec_SIV <- c(sobol_vec_SIV, y1 + slope * (x_value - as.numeric(unlist(x[[i]][1]))[idx]))
    
  }
  
  soboldf[,2] <- sobol_vec_SIV
  output_list[[paste(unlist(colnames(x[[i]][1])))]] <- list(soboldf)
  
}


















