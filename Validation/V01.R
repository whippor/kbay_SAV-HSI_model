################################
### V01. Validation Sampling ###
################################

# reproducibility
renv::restore()

# clear environment
rm(list = ls())

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               terra, 
               tidyterra,
               viridis,
               leaflet,
               ggrepel)

#####################################
#####################################

# load data
understorey <- terra::rast("data/d_suitability_data/understorey_HSI/understorey_HSI.tif")

canopy <- terra::rast("data/d_suitability_data/canopy_HSI/canopy_HSI.tif")

seagrass <- terra::rast("data/d_suitability_data/seagrass_HSI/seagrass_HSI.tif")

bathymetry_under <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")

bathymetry_canop <- terra::rast("data/b_intermediate_data/canopy_bathymetry/bathymetry.grd")

bathymetry_seagr <- terra::rast("data/b_intermediate_data/seagrass_bathymetry/bathymetry.grd")

MLLWcor <- 1.546 # determined w/ https://vdatum.noaa.gov/vdatumweb/

#####################################
#####################################


# generate set of stratified random points for each SAV and save to csv
# if points already saved as csv, skip to next section

# seagrass

bathmask <- clamp(bathymetry_seagr, lower = -10 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
clamped_final <- mask(seagrass, 
                      bathmask)
pts <- as.points(clamped_final, na.rm = TRUE)
pts <- project(pts, "+proj=longlat")
lonlatpts <- crds(pts)
seagrass_latlon <- data.frame(lonlatpts)
seagrass_HSI <- as.data.frame(clamped_final)
HSI_coords <- seagrass_latlon %>%
  mutate(HSI = seagrass_HSI$HSI_value)

rounded_HSI <- HSI_coords %>% 
  mutate(roundHSI = round(HSI, 1)) %>% # round HSI to nearest 0.1
  mutate(bayside = case_when(y < 59.599 ~ "South", # roughly stratify the bay N-S
                             y >= 59.599 ~ "North"))

sampled_HSI <- rounded_HSI %>%
  group_by(bayside, roundHSI) %>%
  sample_n(3, replace = TRUE) %>% # there are only 2 0.4 rounded HSI cell in South
  bind_rows(data.frame("x" = c(rep(NA, 11)),
                       "y" = c(rep(NA, 11)),
                       "HSI" = c(rep(NA, 11)),
                       "roundHSI" = c(seq(from = 0, to = 1, by = 0.1))))

sampled_HSI$site <- 1:nrow(sampled_HSI)
locations <- vect(sampled_HSI, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(sampled_HSI, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Seagrass HSI - Validation Sites")

# write_csv(filter(sampled_HSI, !is.na(HSI)), "data/e_validation_data/seagrass_validation/SeagrassHSIValidation.csv")  


# understorey

bathmask <- clamp(bathymetry_under, lower = -30 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
clamped_final <- mask(understorey, 
                      bathmask)
pts <- as.points(clamped_final, na.rm = TRUE)
pts <- project(pts, "+proj=longlat")
lonlatpts <- crds(pts)
under_latlon <- data.frame(lonlatpts)
under_HSI <- as.data.frame(clamped_final)
HSI_coords <- under_latlon %>%
  mutate(HSI = under_HSI$mean)

rounded_HSI <- HSI_coords %>%
  mutate(roundHSI = round(HSI, 1)) %>% # round HSI to nearest 0.1
  mutate(bayside = case_when(y < 59.599 ~ "South", # roughly stratify the bay N-S
                             y >= 59.599 ~ "North"))

sampled_HSI <- rounded_HSI %>%
  group_by(bayside, roundHSI) %>%
  sample_n(3)

sampled_HSI$site <- 1:nrow(sampled_HSI)
locations <- vect(sampled_HSI, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(sampled_HSI, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Understorey HSI - Validation Sites")

# write_csv(filter(sampled_HSI, !is.na(HSI)), "data/e_validation_data/understorey_validation/UnderstoreyHSIValidation.csv")


# canopy

bathmask <- clamp(bathymetry_canop, lower = -30 - MLLWcor, upper = 0 - MLLWcor, values = FALSE)
clamped_final <- mask(canopy, 
                      bathmask)
pts <- as.points(clamped_final, na.rm = TRUE)
pts <- project(pts, "+proj=longlat")
lonlatpts <- crds(pts)
canop_latlon <- data.frame(lonlatpts)
canop_HSI <- as.data.frame(clamped_final)
HSI_coords <- canop_latlon %>%
  mutate(HSI = canop_HSI$HSI_value)

rounded_HSI <- HSI_coords %>%
  mutate(roundHSI = round(HSI, 1)) %>% # round HSI to nearest 0.1
  mutate(bayside = case_when(y < 59.599 ~ "South", # roughly stratify the bay N-S
                             y >= 59.599 ~ "North"))

sampled_HSI <- rounded_HSI %>%
  group_by(bayside, roundHSI) %>%
  sample_n(3, replace = TRUE) %>% # there is only one 0.3 rounded HSI cell in North
  bind_rows(data.frame("x" = c(rep(NA, 11)),
                       "y" = c(rep(NA, 11)),
                       "HSI" = c(rep(NA, 11)),
                       "roundHSI" = c(seq(from = 0, to = 1, by = 0.1)))) %>%
  distinct()

sampled_HSI$site <- 1:nrow(sampled_HSI)
locations <- vect(sampled_HSI, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(sampled_HSI, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Canopy HSI - Validation Sites")

# write_csv(filter(sampled_HSI, !is.na(HSI)), "data/e_validation_data/canopy_validation/CanopyHSIValidation.csv")


#####################################
#####################################


# all locations

SGval <- read_csv("data/e_validation_data/seagrass_validation/SeagrassHSIValidation.csv")
SGval$type <- "S"
Canval <- read_csv("data/e_validation_data/canopy_validation/CanopyHSIValidation.csv")
Canval$type <- "C"
Undval <- read_csv("data/e_validation_data/understorey_validation/UnderstoreyHSIValidation.csv")
Undval$type <- "U"


# canopy

bathmask <- clamp(bathymetry_canop, lower = -30 - MLLWcor, upper = 0 - MLLWcor, values = FALSE)
clamped_final <- mask(canopy, 
                      bathmask)

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")
locations <- vect(Canval, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(Canval, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Canopy HSI - Validation Sites")


# understorey

bathmask <- clamp(bathymetry_under, lower = -30 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
clamped_final <- mask(understorey, 
                      bathmask)

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")
locations <- vect(Undval, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(Undval, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Understorey HSI - Validation Sites")


# seagrass

bathmask <- clamp(bathymetry_seagr, lower = -10 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
clamped_final <- mask(seagrass, 
                      bathmask)

clamped_proj <- project(clamped_final, "+proj=longlat +datum=WGS84")
locations <- vect(SGval, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = clamped_proj) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(SGval, !is.na(HSI)),
                  aes(x = x, y = y, label = site), 
                  col = "red",
                  size = 5,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("HSI", na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "Seagrass HSI - Validation Sites")


# combined

# EPSG:3338 is NAD83 / Alaska Albers (https://epsg.io/3338)
crs <- "EPSG:3338"

# define vector for region of interest
roi <- terra::vect("data/a_raw_data/LDA_2016.kml")
roi <- project(roi, crs)
bathymetry <- terra::rast("data/b_intermediate_data/bathymetry/bathymetry.grd")

# constrain bathymetry to roi
bathy_roi <- terra::crop(bathymetry, roi, mask = TRUE) 
bathmask <- clamp(bathy_roi, lower = -30 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)


# join master sampling map

Validation_master <- SGval %>%
  bind_rows(Canval, Undval) %>%
  unite("location", c(type, site), sep = "_")

locations <- vect(Validation_master, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

ggplot() +
  geom_spatraster(data = project(bathmask, "+proj=longlat +datum=WGS84")) +
  geom_spatvector(data = locations, col = "red", size = 3) +
  geom_spatvector(data = filter(locations, !is.na(HSI)), aes(color = roundHSI)) +
  geom_text_repel(data = filter(Validation_master, !is.na(HSI)),
                  aes(x = x, y = y, label = location), 
                  col = "red",
                  size = 3,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  scale_fill_viridis("Depth (m)", 
                     option = "mako",
                     na.value = NA) +
  scale_color_viridis(guide = "none") +
  theme_minimal() +
  labs(title = "ALL SAV HSI - Validation Sites")

# write_csv(Validation_master, "data/e_validation_data/AllHSIValidation.csv")

