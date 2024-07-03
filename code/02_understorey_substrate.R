#################################
### 2. Understorey Substrate  ###
#################################

# clear environment
rm(list=setdiff(ls(), c("all_begin", "master_begin")))

# calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, # is replacing the raster package
               viridis)

#####################################
#####################################

# set directories
## define data directory (as this is an R Project, pathnames are simplified)
### input directories
data_dir <- "data/a_raw_data/substrate"

### output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### substrate directory
dir.create(paste0(intermediate_dir, "/",
                  "understorey_substrate"))

substrate_dir <- "data/b_intermediate_data/understorey_substrate"

#####################################
#####################################

# set parameters

## coordinate reference system
### EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)

# import bathymetry as base raster
bathymetry <- terra::rast("~/git/kbay_SAV-HSI_model/data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")

#####################################
#####################################

# load data
substrate <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/Kachemak_Subtidal_Benthic_Habitats.SHP/Kachemak_Subtidal_Benthic_Habitats.shp")
intertidal_segs <- terra::vect("~/git/kbay_SAV-HSI_model/data/a_raw_data/tidalbands_shore_information_translated/tidalbands_shore_information_translatedPolygon.shp")

## reproject into Alaska Albers
subs_albers <- project(substrate, crs)
segs_albers <- project(intertidal_segs, crs)

## reduce to fewer categories
subs_df <- data.frame(subs_albers)
subs_df <- subs_df %>% 
  mutate(A = coalesce(Sub_subgr, Sub_grp)) %>%
  mutate(substrate = coalesce(A, Class))
avoid <- c("Void", 
           "Shell 50-90%, Cobble/gravel 0-10%",
           "Shell 90-100%")
subs_df$substrate <- replace(subs_df$substrate, 
                             subs_df$substrate %in% avoid, 
                             "Unclassified")
subs_albers[["substrate"]] <- subs_df$substrate

segs_df <- data.frame(segs_albers)
segs_df <- segs_df %>%
  mutate(subclass = case_when(subclass == "Rubble" ~ "Boulder",
                              subclass == "Cobble/Gravel" ~ "Cobble/Pebble",
                              subclass == "Mud/Organic" ~ "Mud",
                              is.na(subclass) ~ "Unclassified",
                              .default = subclass))
segs_albers[["subclass"]] <- segs_df$subclass


## make polygons into raster
subs_rast <- rasterize(subs_albers, bathymetry, "substrate")

segs_rast <- rasterize(segs_albers, bathymetry, "subclass")

# inspect the data
## coordinate reference system
terra::crs(subs_rast) # EPSG:3338
cat(crs(subs_rast))

terra::crs(segs_rast)
cat(crs(segs_rast))

## resolution
terra::res(subs_rast) # 50 50
terra::res(segs_rast)

#####################################
#####################################

# expand substrate rasters to bathymetry layer
subs_expand <- terra::extend(subs_rast, bathymetry)
varnames(subs_expand) <- "substrate"


# fill empty space in raster with Unclassified
unclass_rast <- rast(ncol = 1064,
                     nrow = 994,
                     xmin = 119250,
                     xmax = 172450,
                     ymin = 1046527, 
                     ymax = 1096227,
                     crs = crs(subs_expand))
levels(unclass_rast) <- "Unclassified" 
values(unclass_rast) <- 5
names(unclass_rast) <- "substrate"
subs_final <- terra::merge(subs_expand, unclass_rast)

# merge substrate with intertidal segment raster and maintain land
plot(subs_final)


# fill empty space in raster with Unclassified
names(unclass_rast) <- "subclass"
segs_expanded <- terra::merge(segs_rast, unclass_rast)


segs_df <- data.frame(segs_expanded)
subs_df <- data.frame(subs_final)
joined_df <- bind_cols(segs_df, subs_df)
joined_df <- joined_df %>%
  mutate(substrate = case_when(subclass == "Mud" ~ "Mud",
                               subclass == "Sand" ~ "Sand",
                               subclass == "Boulder" ~ "Boulder",
                               subclass == "Cobble/Pebble" ~ "Cobble/Pebble",
                               subclass == "Bedrock" ~ "Bedrock",
                               .default = substrate))
values(segs_rast) <- joined_df$substrate
plot(segs_rast, label = TRUE)
plot(subs_final, add = TRUE)

plet(segs_rast)

click()
#####################################
#####################################

# export raster file
terra::writeRaster(all_final, filename = file.path(substrate_dir, "substrate.tif"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate


