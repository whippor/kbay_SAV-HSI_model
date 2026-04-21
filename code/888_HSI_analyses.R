#########################
### 888. HSI analyses ###
#########################

# reproducibility
renv::restore()

# clear environment
rm(list = ls())

#####################################
#####################################

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,
  terra, # is replacing the raster package
  viridis,
  leaflet,
  ggrepel,
  tidyterra
)

#####################################
#####################################

# load data
roi_dir <- "data/b_intermediate_data/roi"
roi <- terra::vect(roi_dir)

understorey <- terra::rast("data/d_suitability_data/understorey_HSI/understorey_HSI.tif")

canopy <- terra::rast("data/d_suitability_data/canopy_HSI/canopy_HSI.tif")

seagrass <- terra::rast("data/d_suitability_data/seagrass_HSI/seagrass_HSI.grd")

allmax <- terra::rast("data/d_suitability_data/all_HSI/allmax_HSI.tif")

bathymetry_under <- terra::rast("data/b_intermediate_data/understorey_bathymetry/bathymetry.grd")

bathymetry_canop <- terra::rast("data/b_intermediate_data/canopy_bathymetry/bathymetry.grd")

bathymetry_seagr <- terra::rast("data/b_intermediate_data/seagrass_bathymetry/bathymetry.grd")

presence <- terra::rast("data/c_submodel_data/seagrass_presence_HSI/presenceHSI.grd")

MLLWcor <- 1.573 # determined w/ https://dggs.alaska.gov/hazards/coastal/ak-tidal-datum-faq.html

#####################################
#####################################

# canopy stats
# mask out all HSI below -30 m, above 0 m (w/  - 1.546 m correction for NAVD88/MLLW
# bathmask <- clamp(bathymetry_canop, lower = -30 - MLLWcor, upper = 0 - MLLWcor, values = FALSE)
bathmask <- clamp(bathymetry_canop, lower = -30, upper = 0, values = FALSE) # Correction already made in submodel?
bathmask <- crop(
  bathmask,
  roi
)
canopy30 <- mask(
  canopy,
  bathmask
)
plot(canopy30)


plet(canopy30,
  # tiles = "Stadia.AlidadeSmooth",
  # tiles = "USGS.USTopo",
  tiles = "Stadia.StamenTerrain",
  main = "Canopy Kelp HSI"
)

pal <- colorNumeric(palette = "viridis", domain = values(canopy30), na.color = "transparent")
leaflet(options = leafletOptions(zoomSnap = 0.001, zoomDelta = 0.001)) |>
  addTiles() |>
  addRasterImage(x = canopy30, colors = pal) |>
  addLegend(pal = pal, values = values(canopy30), title = "Canopy Kelp HSI", position = "bottomright", opacity = 1)

# canopy stats
terra::expanse(canopy) # 926.84 km^2 (total bay)
terra::hist(canopy)

terra::expanse(canopy30) # 350.25 km^2 (above -30)


# canopy area with HSI > 0

# Calculate the area of each cell in km^2
cell_areas <- cellSize(canopy30, unit = "km")

# Define the specific value you are interested in (e.g., values > 400)
# Create a mask for cells that meet your criteria
specific_value_mask <- ifel(canopy30 > 0, 1, NA)

# Mask the cell areas raster to keep only the areas of the cells you want to sum
masked_areas <- mask(cell_areas, specific_value_mask)

# Sum the areas of the remaining cells
total_area <- global(masked_areas, "sum", na.rm = TRUE)

# Print the total area
print(total_area)


terra::hist(canopy30, maxcell = 1100000)

## proportion suitability
under_df <- data.frame(canopy30)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x = HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

# printed frequencies
cantab <- as.data.frame(table(cut(under_df$HSI_value,
  breaks = seq(0, 1, by = 0.05),
  include.lowest = TRUE
)))
# values between 0.10 - 0.65
sum(cantab[3:13, 2]) # 2262
sum(cantab[, 2])
# total cells: 140098
# proportion
# 2262/140098 # 0.02
allones <- under_df %>%
  filter(HSI_value == 1)

# 15466 canopy HSI 1 from presence

#####################################
#####################################

# understorey stats
# mask out all HSI below -30 m, above 3 m
# bathmask <- clamp(bathymetry_under, lower = -30 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE) # alrerady corrected in submodel?
bathmask <- crop(
  bathmask,
  roi
) # recrop to same extent
understorey30 <- mask(
  understorey,
  bathmask
)
plot(understorey30)
plet(understorey30,
  # tiles = "Stadia.AlidadeSmooth",
  tiles = "Stadia.StamenTerrain",
  main = "Understorey Kelp HSI",
  range = c(-0.001, 1.0001),
  breaks = seq(-0.001, 1.001, by = 0.01)
)
pal <- colorNumeric(palette = "viridis", domain = values(understorey30), na.color = "transparent")
leaflet(options = leafletOptions(zoomSnap = 0.001, zoomDelta = 0.001)) |>
  addTiles() |>
  addRasterImage(x = understorey30, colors = pal)
# addLegend(pal = pal, values = values(understorey30), title = "Understorey Kelp HSI", position = "bottomright", opacity = 1)

# understorey stats
terra::expanse(understorey) # 926.29 km^2 (total bay)
terra::hist(understorey, maxcell = 1100000)

terra::expanse(understorey30) # 354.85 km^2 (above -30)

# understorey area with HSI > 0

# Calculate the area of each cell in km^2
cell_areas <- cellSize(understorey30, unit = "km")

# Define the specific value you are interested in (e.g., values > 400)
# Create a mask for cells that meet your criteria
specific_value_mask <- ifel(understorey30 > 0, 1, NA)

# Mask the cell areas raster to keep only the areas of the cells you want to sum
masked_areas <- mask(cell_areas, specific_value_mask)

# Sum the areas of the remaining cells
total_area <- global(masked_areas, "sum", na.rm = TRUE)

# Print the total area
print(total_area)


terra::hist(understorey30, maxcell = 1100000)

## proportion suitability
under_df <- data.frame(understorey30)
under_df %>%
  mutate(HSIround = round(`HSI_value`, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x = `HSI_value`)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

# printed frequencies
undtab <- as.data.frame(table(cut(under_df$`HSI_value`,
  breaks = seq(0, 1, by = 0.05),
  include.lowest = TRUE
)))

# overlap of HSI = 1 understorey with canopy HSI = 1
kelp_overlap <- c(canopy, understorey)
kelpDF <- data.frame(kelp_overlap)
kelpDF |>
  pivot_longer(HSI_value:HSI_value.1, names_to = "category", values_to = "HSI") |>
  mutate(category = case_when(
    category == "HSI_value" ~ "canopy",
    category == "HSI_value.1" ~ "understorey"
  )) |>
  filter(HSI == 1) %>%
  tally(HSI == 1) # 16746 cells
canopyDF <- data.frame(canopy)
test <- canopyDF %>%
  filter(HSI_value == 1) %>%
  nrow() # 16746 cells
# 5186/19430  ~27%

#####################################
#####################################

# seagrass stats
# mask out all HSI below -7 m, above 3 m
# bathmask <- clamp(bathymetry_seagr, lower = -10 - MLLWcor, upper = 3 - MLLWcor, values = FALSE)
bathmask <- clamp(bathymetry_seagr, lower = -10, upper = 3, values = FALSE) # alrerady corrected in submodel?
bathcrop <- crop(bathmask, seagrass)
seagrass3 <- mask(
  seagrass,
  bathcrop
)
plot(seagrass3)
plet(seagrass3,
  # tiles = "Stadia.AlidadeSmooth",
  tiles = "Stadia.StamenTerrain",
  main = "Seagrass HSI"
)
pal <- colorNumeric(palette = "viridis", domain = values(seagrass3), na.color = "transparent")
leaflet(options = leafletOptions(zoomSnap = 0.001, zoomDelta = 0.001)) |>
  addTiles() |>
  addRasterImage(x = seagrass3, colors = pal)
# addLegend(pal = pal, values = values(seagrass3), title = "Seagrass HSI", position = "bottomright", opacity = 1)


# seagrass stats
terra::expanse(seagrass) # 926.28 km^2 (total bay)
terra::hist(seagrass)


terra::expanse(seagrass3) # 174.90 km^2 (3 to -10)
terra::expanse(seagrass3)

# seagrass area with HSI > 0

# Calculate the area of each cell in km^2
cell_areas <- cellSize(seagrass3, unit = "km")

# Define the specific value you are interested in (e.g., values > 400)
# Create a mask for cells that meet your criteria
specific_value_mask <- ifel(seagrass3 > 0, 1, NA)

# Mask the cell areas raster to keep only the areas of the cells you want to sum
masked_areas <- mask(cell_areas, specific_value_mask)

# Sum the areas of the remaining cells
total_area <- global(masked_areas, "sum", na.rm = TRUE)

# Print the total area
print(total_area)


terra::hist(seagrass3, maxcell = 1100000)

## proportion suitability
under_df <- data.frame(seagrass3)
under_df %>%
  mutate(HSIround = round(HSI_value, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x = HSI_value)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()

# printed frequencies
seatab <- as.data.frame(table(cut(under_df$HSI_value,
  breaks = seq(0, 1, by = 0.05),
  include.lowest = TRUE
)))
# values between 0.25 - 0.75
sum(seatab[6:15, 2]) # 17769
sum(seatab[, 2])
# total cells: 69959
# proportion
# 17769/69959 # 0.25
nonzero <- under_df %>%
  filter(HSI_value != 0)
# 42183


allones <- under_df %>%
  filter(HSI_value == 1)

# 6057 seagrass HSI 1 from presence


# HSI 1 contributions
all1 <- data.frame(presence[presence == 1]) # 4295
seag1 <- data.frame(seagrass3[seagrass3 == 1]) # 6057
# 4295/6057 # 0.71

#####################################
#####################################

# mask out all HSI below -30 m, above 3 m
bathmask <- clamp(bathymetry_under, lower = -30, upper = 3, values = FALSE)
allmax30 <- mask(
  allmax,
  bathmask
)
plot(allmax30)
plet(allmax30,
  # tiles = "Stadia.AlidadeSmooth",
  tiles = "Stadia.StamenTerrain",
  main = "All SAV HSI"
)
pal <- colorNumeric(palette = "viridis", domain = values(allmax30), na.color = "transparent")
leaflet(options = leafletOptions(zoomSnap = 0.001, zoomDelta = 0.001)) |>
  addTiles() |>
  addRasterImage(x = allmax30, colors = pal)
# addLegend(pal = pal, values = values(allmax30), title = "All SAV HSI", position = "bottomright", opacity = 1)


# combined stats
terra::expanse(allmax) # 946.16 km^2 (total bay)
terra::hist(allmax)

terra::expanse(allmax30) # 383.09 km^2 (above -30)


# ALL area with HSI > 0

# Calculate the area of each cell in km^2
cell_areas <- cellSize(allmax30, unit = "km")

# Define the specific value you are interested in (e.g., values > 400)
# Create a mask for cells that meet your criteria
specific_value_mask <- ifel(allmax30 > 0, 1, NA)

# Mask the cell areas raster to keep only the areas of the cells you want to sum
masked_areas <- mask(cell_areas, specific_value_mask)

# Sum the areas of the remaining cells
total_area <- global(masked_areas, "sum", na.rm = TRUE)

# Print the total area
print(total_area)


terra::hist(allmax30, maxcell = 1100000)

## proportion suitability
under_df <- data.frame(allmax30)
under_df %>%
  mutate(HSIround = round(`mean`, digits = 2)) %>%
  group_by(HSIround) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggplot(under_df, aes(x = `mean`)) +
  geom_histogram(binwidth = 0.05) +
  theme_bw()


# printed frequencies
alltab <- as.data.frame(table(cut(under_df$HSI_value,
  breaks = seq(0, 1, by = 0.05),
  include.lowest = TRUE
)))
allbay <- data.frame(understorey)
sum(alltab[2:20, 2]) # 119611
# 119611/149562 # 0.799
# 119611/378464 # 0.3160


max30DF <- data.frame(allmax30)
max30DF %>%
  filter(HSI_value > 0.5) %>%
  nrow() # 114522
114522 / 149562 # 0.765




plot(allmax, col = viridis(nrow(allmax), begin = 0.3))

plet(allmax,
  main = "SAV max HSI"
)


#####################################
#####################################

## NOTE!!!!
# must change all hist() to include maxcell = large number to get all values!


############### HIST of diff between arithmetic and geometric means for understorey

# load data
substrate <- terra::rast("data/c_submodel_data/understorey_substrate_HSI/substrateHSI.grd")

bathymetry <- terra::rast("data/c_submodel_data/understorey_bathymetry_HSI/bathymetryHSI.grd")

fetch <- terra::rast("data/c_submodel_data/understorey_fetch_HSI/fetchHSI.grd")



# create arithmetic mean-value HSI for understorey (bathymetry, substrate, fetch)
under_mean <- mean(
  substrate[["HSI_value"]],
  bathymetry[["HSI_value"]],
  fetch[["HSI_value"]]
  # , velocity[["HSI_value"]]
)


# create geometric mean-value HSI for understorey (bathymetry, substrate, fetch)
under_mean <- exp(mean(log(c(
  (substrate[["HSI_value"]] + 1),
  (bathymetry[["HSI_value"]] + 1),
  (fetch[["HSI_value"]] + 1)
))))


mean_plot <- data.frame(under_mean)
ggplot(mean_plot, aes(x = `HSI_value`)) +
  geom_histogram(binwidth = 0.05) +
  ggtitle("Arithmetic") +
  theme_bw()


zero_plot <- data.frame(under_zero)
ggplot(zero_plot, aes(x = `mean`)) +
  geom_histogram(binwidth = 0.05) +
  ggtitle("Geometric") +
  theme_bw()

final_plot <- data.frame(final_zero)
ggplot(final_plot, aes(x = `mean`)) +
  geom_histogram(binwidth = 0.05) +
  ggtitle("Geometric, final") +
  theme_bw()


# create HSI setting any zero to zero
under_zero <- (under_mean - 1)
sub_zero <- substrate
sub_zero[sub_zero < 0.01] <- 0
sub_zero[sub_zero > 0.01] <- 1
bath_zero <- bathymetry
bath_zero[bath_zero < 0.01] <- 0
bath_zero[bath_zero > 0.01] <- 1
fetch_zero <- fetch
fetch_zero[fetch_zero < 0.01] <- 0
fetch_zero[fetch_zero > 0.01] <- 1
final_zero <- under_zero[["mean"]] *
  sub_zero[["HSI_value"]] *
  bath_zero[["HSI_value"]] *
  fetch_zero[["HSI_value"]]


#####################
#####################

library(ggplot2)
library(hexbin)

dat <- data.frame(x = rnorm(10000), y = rnorm(10000))

ggplot(dat, aes(x = x, y = y)) +
  geom_hex() +
  coord_fixed() +
  scale_fill_gradientn(colours = viridis(256, option = "B"))

# using code from RColorBrewer to demo the palette
n <- 200
image(
  1:n, 1, as.matrix(1:n),
  col = viridis(n, option = "D"),
  xlab = "viridis n", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)


tamdata <- data.frame(
  continuous = c(10, 5, 3, 1),
  value = c(0, 1, 1, 0)
)

ggplot(tamdata, aes(x = continuous, y = value)) +
  geom_line(linewidth = 8) +
  geom_point(color = "blue", size = 20) +
  ylab("habitat suitability") +
  theme(text = element_text(size = 64))


tamdata2 <- data.frame(
  discrete = c("A", "B", "C", "D"),
  value = c(0.5, 1, 0.3, 0.15)
)

ggplot(tamdata2, aes(x = discrete, y = value, fill = discrete)) +
  geom_col() +
  scale_fill_viridis(
    discrete = TRUE,
    option = "mako",
    begin = 0.2,
    end = 0.8
  ) +
  theme_linedraw() +
  ylab("habitat suitability") +
  theme(text = element_text(size = 64)) +
  guides(fill = "none")


ggplot() +
  geom_spatraster(data = bathymetry, aes(fill = HSI_value)) + # Map raster values to fill
  scale_fill_viridis_c(
    option = "rocket",
    begin = 0.65,
    end = 0,
    direction = -1,
    na.value = "white"
  ) + # Use a viridis color scale
  theme_minimal() +
  guides(fill = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

ggplot() +
  geom_spatraster(data = substrate, aes(fill = HSI_value)) + # Map raster values to fill
  scale_fill_viridis_c(
    option = "rocket",
    begin = 0.65,
    end = 0,
    direction = -1,
    na.value = "white"
  ) + # Use a viridis color scale
  theme_minimal() +
  guides(fill = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

ggplot() +
  geom_spatraster(data = fetch, aes(fill = HSI_value)) + # Map raster values to fill
  scale_fill_viridis_c(
    option = "rocket",
    begin = 0.65,
    end = 0,
    direction = -1,
    na.value = "white"
  ) + # Use a viridis color scale
  theme_minimal() +
  guides(fill = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

ggplot() +
  geom_spatraster(data = presence, aes(fill = HSI_value)) + # Map raster values to fill
  scale_fill_viridis_c(
    option = "rocket",
    begin = 0.65,
    end = 0,
    direction = -1,
    na.value = "white"
  ) + # Use a viridis color scale
  theme_minimal() +
  guides(fill = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )


pal <- colorNumeric(
  palette = viridis::viridis(100), # Generate 100 colors from the viridis palette
  domain = values(canopy30),
  na.color = "transparent" # Handle NA values
)


my_map <- leaflet() %>%
  addProviderTiles(providers$Stadia.StamenTerrain) %>%
  setView(lng = -151.41, lat = 59.60, zoom = 10) %>%
  addRasterImage(
    x = canopy30,
    colors = pal,
    opacity = 0.8,
    project = TRUE # Ensure it's projected correctly for Leaflet
  )

mapshot(
  x = my_map,
  file = "canopy_large.png",
  vwidth = 8000,
  vheight = 6000,
  zoom = 2
)


cellSize(seagrass)
