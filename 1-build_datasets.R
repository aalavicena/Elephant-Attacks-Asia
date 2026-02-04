# =============================================================================
# 1. Build datasets
# Creates occurrence and background point datasets for HEC modelling.
# =============================================================================

# Configure base path: set to your data root (e.g. "E:/Spatial Analysis" or "/Volumes/MEC/Spatial Analysis")
BASE_PATH <- "E:/Spatial Analysis"
PATH_VARS <- file.path(BASE_PATH, "_Working Variables")
PATH_DATA <- file.path(BASE_PATH, "data")
PATH_OUT  <- BASE_PATH

library(terra)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Load and clean occurrence data
# -----------------------------------------------------------------------------
nstack <- rast(file.path(PATH_VARS, "stack/stack-250717.tif"))
df <- read.csv(file.path(PATH_DATA, "250529_data.csv"))

df <- df %>%
  select(longitude, latitude, human_death, human_injury) %>%
  drop_na(longitude, latitude)

df$human_death <- as.factor(df$human_death)
levels(df$human_death)[levels(df$human_death) == ""] <- NA
df$human_injury <- as.factor(df$human_injury)

deaths <- df %>%
  select(longitude, latitude, human_death) %>%
  filter(!is.na(human_death))

injuries <- df %>%
  select(longitude, latitude, human_injury) %>%
  filter(!is.na(human_injury))

# -----------------------------------------------------------------------------
# Extract raster values at occurrence points
# -----------------------------------------------------------------------------
occ_deaths   <- terra::extract(nstack, as.matrix(deaths[, c("longitude", "latitude")]), na.rm = TRUE, xy = TRUE)
occ_injuries <- terra::extract(nstack, as.matrix(injuries[, c("longitude", "latitude")]), na.rm = TRUE, xy = TRUE)
occ_all      <- terra::extract(nstack, as.matrix(df[, c("longitude", "latitude")]), na.rm = TRUE, xy = TRUE)

# Drop ID column from extract(); bind outcomes for occ_all
occ_deaths   <- occ_deaths[, -1L, drop = FALSE]
occ_injuries <- occ_injuries[, -1L, drop = FALSE]
occ_all      <- occ_all[, -1L, drop = FALSE]
occ_all      <- bind_cols(occ_all, df[, c("human_death", "human_injury")])

# -----------------------------------------------------------------------------
# Background points
# -----------------------------------------------------------------------------
boundary_path <- file.path(PATH_VARS, "Boundary/IUCN-redlist-100km-dissolved.shp")
vec <- vect(boundary_path)
vec <- project(vec, nstack)

set.seed(123)
bg_sample <- spatSample(vec, 7000L, method = "random")
bg <- terra::extract(nstack, bg_sample, na.rm = TRUE, xy = TRUE)
bg <- bg[, -1L, drop = FALSE]
bg <- na.omit(bg)

n_deaths   <- nrow(occ_deaths)
n_injuries <- nrow(occ_injuries)
n_presence <- nrow(occ_all)
n_bg       <- 3L * max(n_deaths, n_injuries, n_presence)

set.seed(123)
bg_d   <- bg[sample(nrow(bg), 3L * n_deaths), ]
bg_i   <- bg[sample(nrow(bg), 3L * n_injuries), ]
bg_all <- bg[sample(nrow(bg), n_bg), ]
bg_all$human_death <- NA
bg_all$human_injury <- NA

# -----------------------------------------------------------------------------
# Combine presence and background; save
# -----------------------------------------------------------------------------
d_d   <- rbind(cbind(att = 1L, occ_deaths),   cbind(att = 0L, bg_d))
d_i   <- rbind(cbind(att = 1L, occ_injuries), cbind(att = 0L, bg_i))
d_all <- rbind(cbind(att = 1L, occ_all),      cbind(att = 0L, bg_all))

write.csv(d_d,   file.path(PATH_OUT, "250617_deaths_datasets.csv"),   row.names = FALSE)
write.csv(d_i,   file.path(PATH_OUT, "250617_injuries_datasets.csv"), row.names = FALSE)
write.csv(d_all, file.path(PATH_OUT, "250617_all_datasets.csv"),       row.names = FALSE)
