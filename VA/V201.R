########################################
### V201. Field Validation Analysis ####
########################################

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


#####################################
#####################################

# import all HSI maps
canopy <- terra::rast("data/d_suitability_data/canopy_HSI/canopy_HSI.tif")
understorey <- terra::rast("data/d_suitability_data/understorey_HSI/understorey_HSI.tif")
seagrass <- terra::rast("data/d_suitability_data/seagrass_HSI/seagrass_HSI.tif")

# import SAV observations
kbay_seagrass_pa <- read_csv("~/git/kbay_seagrass_monitoring/data/seagrass_density.csv")

SAV_obs <- read_csv("data/e_validation_data/field_obs_SAV.csv") |>
  full_join(kbay_seagrass_pa)

# convert observations to spat vector


SAV_vect <- vect(SAV_obs, geom = c("lon", "lat"), crs = "epsg:4326")
values(SAV_vect) <- SAV_obs[, "qual_dens"]
crs <- crs(seagrass)
SAV_vect <- project(SAV_vect, crs)

plot(seagrass)
plot(SAV_vect, add = TRUE)

#####################################
#####################################

# SEAGRASS validation

# extract HSI values for points
HSI_seagrass_pts <- extract(seagrass, SAV_vect)
HSI_latlon <- geom(SAV_vect)
HSI_joined <- HSI_latlon |>
  bind_cols(HSI = HSI_seagrass_pts$HSI_value) |>
  mutate(lon = x, lat = y) |>
  select(lon, lat, HSI) 


# join with field obs
HSI_obs <- SAV_obs |>
  full_join(HSI_joined) |>
  mutate(quant_dens = case_when(qual_dens == "none" ~ 0,
                                qual_dens == "trace" ~ 1,
                                qual_dens == "sparse" ~ 2,
                                qual_dens == "thin" ~ 3,
                                qual_dens == "moderate" ~ 4,
                                qual_dens == "thick" ~ 5)) |>
  mutate(presabs = case_when(qual_dens == "none" ~ 0,
                             .default = 1))

ggplot(HSI_obs, aes(x = quant_dens, y = HSI)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(HSI_obs, aes(x = HSI, y = presabs)) +
  geom_jitter() +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE) 

ggplot(HSI_obs, aes(x = as.character(presabs), y = HSI)) +
  geom_count() 

HSI_obs |>
  mutate(overestimate = case_when(HSI > 0 & presabs == 0 ~ 1,
                                  .default = 0)) |>
  mutate(underestimate = case_when(HSI == 0 & presabs > 0 ~ 1,
                                   .default = 0)) |>
  mutate(matched_abs = case_when(HSI == 0 & presabs == 0 ~ 1,
                                 .default = 0)) |>
  mutate(matched_pres = case_when(HSI > 0 & presabs > 0 ~ 1,
                                  .default = 0)) |>
  pivot_longer(overestimate:matched_pres, names_to = "estimate", values_to = "count") |>
  group_by(estimate) |>
  summarise(count = sum(count)) |>
  ggplot() +
  geom_col(aes(x = estimate, y = count)) 



# PULL MAX QUAL VALUE FOR EACH RASTER CELL

# Extract the raster cell ID for each point in the spatvector
cell_ids <- extract(seagrass, SAV_vect, cell=TRUE)[, "cell"]

# Get the associated values from the spatvector
point_values <- values(SAV_vect)

# Combine cell IDs and point values into a dataframe
data <- data.frame(cell_id = cell_ids, value = point_values)

data <- data |>
  mutate(qual_dens = case_when(qual_dens == "none" ~ 0,
                                qual_dens == "trace" ~ 1,
                                qual_dens == "sparse" ~ 2,
                                qual_dens == "thin" ~ 3,
                                qual_dens == "moderate" ~ 4,
                                qual_dens == "thick" ~ 5))

# Calculate the maximum value for each cell ID
max_values_per_cell <- aggregate(qual_dens ~ cell_id, data = data, FUN = max)

# Create a new spatraster to store the maximum values
# First, create a raster with the same properties as the original spat_raster
max_raster <- rast(seagrass)

# Assign the maximum values to the corresponding cells in the new raster
# Need to match the cell IDs from max_values_per_cell to the raster cell IDs
# This requires careful indexing. A common way is to create a vector
# where the index is the cell_id and the value is the max_value
max_value_vector <- rep(NA, ncell(seagrass))
max_value_vector[max_values_per_cell$cell_id] <- max_values_per_cell$qual_dens

values(max_raster) <- max_value_vector

# Now 'max_raster' is a spatraster where each cell contains the maximum
# value of the spatvector points that fall within it.
plet(max_raster)

# Pull out cells that both have values (not NA)
        
# Convert both rasters to dataframes with coordinates
df1 <- as.data.frame(seagrass, xy=TRUE)
df2 <- as.data.frame(max_raster, xy=TRUE)

# Rename the value columns to make them distinct before merging
colnames(df1)[ncol(df1)] <- "value_HSI"
colnames(df2)[ncol(df2)] <- "value_obs"

# Merge the two dataframes based on their x and y coordinates
# Using 'all.x = TRUE' ensures all rows from df1 are kept,
# and matching rows from df2 are added. Non-matching rows from df2 are ignored.
merged_df <- merge(df1, df2, by.x = c("x", "y"), by.y = c("x", "y"), all.x = TRUE)

# Identify rows where raster2's value is not NA
all_values <- merged_df[!is.na(merged_df$value_obs), ]

# plot
presabs <- all_values |>
  mutate(overestimate = case_when(value_HSI > 0 & value_obs == 0 ~ 1,
                                  .default = 0)) |>
  mutate(underestimate = case_when(value_HSI == 0 & value_obs > 0 ~ 1,
                                   .default = 0)) |>
  mutate(matched_abs = case_when(value_HSI == 0 & value_obs == 0 ~ 1,
                                 .default = 0)) |>
  mutate(matched_pres = case_when(value_HSI > 0 & value_obs > 0 ~ 1,
                                  .default = 0)) |>
  pivot_longer(overestimate:matched_pres, names_to = "estimate", values_to = "count") |>
  group_by(estimate) |>
  summarise(count = sum(count)) 

sum(presabs$count)
176/214 

 ggplot(presabs) +
  geom_col(aes(x = estimate, y = count)) 

max_raster
plet(seagrass)
