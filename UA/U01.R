#################################
### U01. Uncertainty Analysis ###
#################################

# reproducibility
renv::restore()

# clear environment
rm(list = ls())

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, # is replacing the raster package
               viridis,
               EnvStats,
               sensobol,
               foreach,
               parallel,
               doParallel,
               sn)
source("code/000_function_interpolate_y.R")
source("code/000_function_unif_vals.R")

#####################################
#####################################

# set parameters
roi_dir <- "data/b_intermediate_data/roi"

#####################################
#####################################

# import region of interest
roi <- terra::vect(roi_dir)

# load data
bathymetry <- terra::rast("data/b_intermediate_data/bathymetry/bathymetry.grd")
bath_mask <- mask(bathymetry, roi)
bath_mask <- clamp(bath_mask, lower = -30, upper = 3, values = FALSE)

substrate <- terra::rast("data/b_intermediate_data/substrate/substrate.tif")
subs_mask <- mask(substrate, roi)
subs_mask <- mask(subs_mask, bath_mask)

fetch <- terra::rast("data/b_intermediate_data/fetch/fetch.grd")
fetch <- fetch/1000
fetch_mask <- mask(fetch, roi)
fetch_mask <- mask(fetch_mask, bath_mask)

## TAM TABLES
under_tam_bath <- read_csv("data/x_tam_tables/understorey/understorey_depth.csv")
under_tam_bath <- under_tam_bath %>%
  arrange(depth.m)

under_tam_subs <- read_csv("data/x_tam_tables/understorey/understorey_substrate.csv")

under_tam_fetch <- read_csv("data/x_tam_tables/understorey/understorey_fetch.csv")
under_tam_fetch <- under_tam_fetch %>%
  arrange(fetch.km)

can_tam_bath <- read_csv("data/x_tam_tables/canopy/canopy_depth.csv")
can_tam_bath <- can_tam_bath %>%
  arrange(depth.m)

can_tam_subs <- read_csv("data/x_tam_tables/canopy/canopy_substrate.csv")

can_tam_fetch <- read_csv("data/x_tam_tables/canopy/canopy_fetch.csv")
can_tam_fetch <- can_tam_fetch %>%
  arrange(fetch.km)

sea_tam_bath <- read_csv("data/x_tam_tables/seagrass/seagrass_depth.csv")
sea_tam_bath <- sea_tam_bath %>%
  arrange(depth.m)

sea_tam_subs <- read_csv("data/x_tam_tables/seagrass/seagrass_substrate.csv")

sea_tam_fetch <- read_csv("data/x_tam_tables/seagrass/seagrass_fetch.csv")
sea_tam_fetch <- sea_tam_fetch %>%
  arrange(fetch.km)

#####################################
#####################################

# extract values from base layers
vals_bath <- data.frame(values(bath_mask, na.rm = TRUE))
vals_bath <- vals_bath %>%
  rename("depth" = "KBL.bathymetry_GWA.area_50m_EPSG3338")

vals_sub <- data.frame(values(subs_mask, na.rm = TRUE))
vals_sub <- vals_sub %>%
  rename("value" = "subclass")

vals_fetch <- data.frame(values(fetch_mask, na.rm = TRUE))
vals_fetch <- vals_fetch %>%
  rename("fetch" = "lyr.1")



x_values <- c(vals_bath, vals_fetch) # for testing function
tam_table <- list(under_tam_bath, under_tam_fetch) # for testing function
samples = 10000
sobol_order = "first"
#####################################
#####################################
              
# create multivariate input variables based on sobol sequences and find SIV 
# from adjusted TAM tables

## Function 
sobol_interpolate_y <- function(x_values, 
                                tam_table, 
                                samples = 10,
                                sobol_order = "first") {
  # output list 
  output_list <<- list()
  
  # name parameters
  all_list <- x_values
  
  
  for (i in 1:length(all_list)) {
    assign(paste("x", i, sep = "_"), all_list[[i]])
  }
  

  x <- tam_table
  
  # check that there is a tam table for each input
  if (length(all_list) == length(x)) {
    # add adjusted SIV values to TAM tables
    x <- map(x, \(x) rename(x, "SIV" = contains("SIV")))
    for (i in seq_along(x)) {
      x[[i]]$unif_rand <- runif(nrow(x[[i]]), min = -0.2, max = 0.2)
      x[[i]]$TAM_adj <- x[[i]]$SIV + x[[i]]$unif_rand
      x[[i]]$TAM_adj[1] <- x[[i]]$SIV[1]
      x[[i]]$TAM_adj[length(x[[i]]$TAM_adj)] <- x[[i]]$SIV[length(x[[i]]$SIV)]
      x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj > 1, 1, x[[i]]$TAM_adj)
      x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj < 0, 0, x[[i]]$TAM_adj)
    }
    
    # create sobol sequence for x values
    z <- length(all_list)
    
    
    for (i in 1:z) {
      N <- samples
      df <- paste("x", i, sep = "_")
      params <- paste0(names(x[[i]][1]), sep = "_", "sobol")
      mat <- sobol_matrices(N = N,
                            params = params,
                            order = sobol_order,
                            type = "R")
      mat[, 1] <- qunif(mat[, 1], min(get(df)), max(get(df)))
      
      soboldf <- data.frame(mat, matrix(NA, nrow = length(mat), ncol = 1))
      soboldf <- rename(soboldf, !!paste(names(x[[i]][1]), "newSIV", sep = "_") := 2)
      
      for (j in 1:nrow(soboldf)) {
        x_value <- soboldf[j, 1]
        
        #### Check if x_value is outside the range of tam_table$x
        ifelse(x_value < min(as.numeric(unlist(x[[i]][1]))), 0, x_value)
        ifelse(x_value > max(as.numeric(unlist(x[[i]][1]))), 0, x_value)
        ifelse(is.na(x_value) == TRUE, 0, x_value)
        # Set interpolated y value to zero
        
        ##### Find the interval where x_value lies
        idx <- findInterval(x_value, as.numeric(unlist(x[[i]][1])))
        ##### Find slopes
        slopes <- diff(x[[i]]$TAM_adj) / diff(unlist(x[[i]][1]))
        
        ##### Linear interpolation
        y1 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx]
        y2 <- as.numeric(unlist(x[[i]]$TAM_adj))[idx + 1]
        slope <- slopes[idx]
        soboldf[j, 2] <- y1 + slope * (x_value - as.numeric(unlist(x[[i]][1]))[idx])
        
      }
      
      output_list[[paste(unlist(colnames(x[[i]][1])))]] <<- list(soboldf)
     # assign(paste(unlist(colnames(x[[i]][1])), i, "sobol_df", sep = "_"), value = soboldf)
      
    }
    
    
    
  } else {
    
    print("Number of parameters do not match number of TAM tables")
  
    }
  
}


sobol_interpolate_y(c(vals_bath, vals_fetch), list(under_tam_bath, under_tam_fetch), 10)
output2 <- output_list
output3 <- output_list

n <- max(length(vals_bath$depth), length(vals_fetch$fetch))
length(vals_bath$depth) <- n                      
length(vals_fetch$fetch) <- n
cbind(vals_bath, vals_fetch)
test_result <- sobol_interpolate_y(vals_bath, under_tam_bath, samples = 1000)


# Add rand uniform to all imported TAM tables

x <- list(sea_tam_bath, sea_tam_fetch) # x for testing function

# keep outside the loop that separates dataframes
## COMPLETE BLOCK - ADDS ADJUSTED SIV VALUE TO TAM TABLES
x <- map(x, \(x) rename(x, "SIV" = contains("SIV")))
for (i in seq_along(x)) {
  x[[i]]$unif_rand <- runif(nrow(x[[i]]), min = -0.2, max = 0.2)
  x[[i]]$TAM_adj <- x[[i]]$SIV + x[[i]]$unif_rand
  x[[i]]$TAM_adj[1] <- x[[i]]$SIV[1]
  x[[i]]$TAM_adj[length(x[[i]]$TAM_adj)] <- x[[i]]$SIV[length(x[[i]]$SIV)]
  x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj > 1, 1, x[[i]]$TAM_adj)
  x[[i]]$TAM_adj <- ifelse(x[[i]]$TAM_adj < 0, 0, x[[i]]$TAM_adj)
}
##


x <- mapply(cbind, x, "unif_rand" = runif(sum(lengths(x), min = -0.2, max = 0.2)))

for (i in 1:length(x)){
  require(tidyverse)
  x <- x %>%
    
  x <- lapply(fun = cbind(x, "unif_rand" = runif(nrow(x[[i]]), min = -0.2, max = 0.2)))
  assign(paste("TAM", i, sep = "_"), as.data.frame(x[[i]]))
  
}



  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  colnames(tam_table)[2] <- "SIV"
  unif_rand <- data.frame(TAM_adj = runif(nrow(tam_table), min = -0.2, max = 0.2))
  tam_20 <- tam_table$SIV + unif_rand 
  tam_20[1,] <- tam_table[1, 2]
  tam_20[nrow(tam_20),] <- tam_table[nrow(tam_table), 2]
  tam_20$TAM_adj <- ifelse(tam_20$TAM_adj > 1, 1, tam_20$TAM_adj)
  tam_20$TAM_adj <- ifelse(tam_20$TAM_adj < 0, 0, tam_20$TAM_adj)
  tam_new <- data.frame(tam_table[,1], tam_20)
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_table$x
    ifelse(x_value < min(as.numeric(unlist(tam_new[,1]))), 0, x_value)  
    ifelse(x_value > max(as.numeric(unlist(tam_new[,1]))), 0, x_value)  
    ifelse(is.na(x_value) == TRUE, 0, x_value) 
    # Set interpolated y value to zero
    
    ##### Find the interval where x_value lies
    idx <- findInterval(x_value, as.numeric(unlist(tam_new[,1])))
    
    ##### Linear interpolation
    y1 <- as.numeric(unlist(tam_new[,2]))[idx]
    y2 <- as.numeric(unlist(tam_new[,2]))[idx + 1]
    slope <- slopes[idx]
    interpolated_y_values[i] <- y1 + slope * (x_value - as.numeric(unlist(tam_new[,1]))[idx])
    
  }
  
  return(interpolated_y_values)
}




test_function <- function(...) {
  x <- list(...)
  for (i in 1:length(x)){
    assign(paste("TAM", i, sep = "_"), as.data.frame(x[[i]]))
    print(x[[i]])
  }
}



test_function(sea_tam_bath, sea_tam_subs)





























# create 20% variation all TAM values

## Function to interpolate y value for given x values
UA_interpolate_y <- function(x_values, tam_table) {
  # ensure dataframe
  x_values <- data.frame(x_values)
  ### Initialize an empty vector for interpolated y values
  interpolated_y_values <- numeric(length(x_values))
  
  colnames(tam_table)[2] <- "SIV"
  unif_rand <- data.frame(TAM_adj = runif(nrow(tam_table), min = -0.2, max = 0.2))
  tam_20 <- tam_table$SIV + unif_rand 
  tam_20[1,] <- tam_table[1, 2]
  tam_20[nrow(tam_20),] <- tam_table[nrow(tam_table), 2]
  tam_20$TAM_adj <- ifelse(tam_20$TAM_adj > 1, 1, tam_20$TAM_adj)
  tam_20$TAM_adj <- ifelse(tam_20$TAM_adj < 0, 0, tam_20$TAM_adj)
  tam_new <- data.frame(tam_table[,1], tam_20)
  
  for (i in seq_along(x_values)) {
    x_value <- x_values[i]
    
    #### Check if x_value is outside the range of tam_table$x
    ifelse(x_value < min(as.numeric(unlist(tam_new[,1]))), 0, x_value)  
    ifelse(x_value > max(as.numeric(unlist(tam_new[,1]))), 0, x_value)  
    ifelse(is.na(x_value) == TRUE, 0, x_value) 
      # Set interpolated y value to zero
    
      ##### Find the interval where x_value lies
      idx <- findInterval(x_value, as.numeric(unlist(tam_new[,1])))
      
      ##### Linear interpolation
      y1 <- as.numeric(unlist(tam_new[,2]))[idx]
      y2 <- as.numeric(unlist(tam_new[,2]))[idx + 1]
      slope <- slopes[idx]
      interpolated_y_values[i] <- y1 + slope * (x_value - as.numeric(unlist(tam_new[,1]))[idx])
    
  }
  
  return(interpolated_y_values)
}



UA_interpolate_y(mat, under_tam_bath)



x_values <- mat

tam_table <- under_tam_bath
#####################################
#####################################


# extract values from base layers
vals_bath <- data.frame(values(bath_mask, na.rm = TRUE))
vals_bath <- vals_bath %>%
  rename("value" = "KBL.bathymetry_GWA.area_50m_EPSG3338")

vals_sub <- data.frame(values(subs_mask, na.rm = TRUE))
vals_sub <- vals_sub %>%
  rename("value" = "subclass")

vals_fetch <- data.frame(values(fetch_mask, na.rm = TRUE))
vals_fetch <- vals_fetch %>%
  rename("value" = "lyr.1")


#####################################
#####################################

# sobol for bathymetry and fetch



# test for possible distributions of data
bath_samp <- sample(vals_bath$value, 3000) # subsample so function can handle
hist(bath_samp) # visually inspect for candidate dists
ggplot(data.frame(bath_samp), aes(x = bath_samp)) +
  geom_boxplot()
distChoose(bath_samp, choices = c("gevd", "logis", "norm"))
# nonparametric - use uniform dist

fetch_samp <- sample(vals_fetch$value, 3000) # subsample so function can handle
hist(fetch_samp) # visually inspect for candidate dists
ggplot(data.frame(fetch_samp), aes(x = fetch_samp)) +
  geom_boxplot()
distChoose(fetch_samp, choices = c("gevd", "logis", "norm"))
# nonparametric - use uniform dist

# set up parameters for sobol matrix
N <- nrow(vals_bath) * 0.3
params <- c("bath")
order <- "first"

mat <- sobol_matrices(N = N, params = params, order = order)

mat[, "bath"] <- qunif(mat[, "bath"], min(vals_bath), max(vals_bath))
mat[, "fetch"] <- qunif(mat[, "fetch"], min(vals_fetch), max(vals_fetch))


library(foreach)
library(parallel)
library(doParallel)


n.cores <- makeCluster(floor(detectCores() * 0.75)) 
registerDoParallel(n.cores)  

stopCluster(n.cores)

hist(mat)



install.packages("EnvStats")
library(EnvStats)



mat <- sobol_matrices(N = N, params = params, order = order)
mat[, "bath"] <- qgevd(mat[, "bath"], location = -10, shape = 0.9)

hist(mat)

median(vals_bath$value, na.rm = TRUE)









# sobol for fetch

N <- 2^8
params <- c("fetch")
order <- "first"

fetch_mat <- sobol_matrices(N = N, params = params, order = order)
fetch_mat[, "fetch"] <- qunif(fetch_mat[, "fetch"], min(vals_fetch), max(vals_fetch))




hist(fetch_mat)




