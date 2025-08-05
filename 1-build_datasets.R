#############################
## Creating datasets

library(terra)
library(dplyr)
library(tidyverse)

############################

nstack <- rast("E:/Spatial Analysis/_Working Variables/stack/stack-250717.tif")
df <- read.csv("E:/Spatial Analysis/data/250529_data.csv")

# clean data
df <- df %>%
  select(longitude, latitude, human_death, human_injury) %>%
  drop_na(longitude, latitude)
summary(df)

df$human_death <- as.factor(df$human_death)
levels(df$human_death)[levels(df$human_death) == ""] <- NA
summary(df)

df$human_injury <- as.factor(df$human_injury)

deaths <- df %>%
  select(longitude, latitude, human_death) %>%
  filter(is.factor(human_death) & !is.na(human_death))
summary(deaths) # 1110

injuries <- df %>%
  select(longitude, latitude, human_injury) %>%
  filter(is.factor(human_injury) & !is.na(human_injury))
summary(injuries) # 347

# Sample raster
set.seed(123)
occ_deaths <- terra::extract(nstack, deaths[,1:2], na.rm = TRUE, xy = TRUE)
occ_deaths <- occ_deaths[,-c(1)]
#plot(occ_deaths[, c('x', 'y')])

occ_injuries <- terra::extract(nstack, injuries[,1:2], na.rm = TRUE, xy = TRUE)
occ_injuries <- occ_injuries[,-c(1)]
#plot(occ_injuries[, c('x', 'y')])

occ_all <- terra::extract(nstack, df[,1:2], na.rm = TRUE, xy = TRUE)
occ_all <- occ_all[,-c(1)]
occ_all <- bind_cols(occ_all, df[, c("human_death", "human_injury")])
#plot(occ_all[, c('x', 'y')])

# Generate background points
set.seed(123)
vec <- vect("E:/Spatial Analysis/_Working Variables/Boundary/IUCN-redlist-100km-dissolved.shp")
vec <- project(vec, nstack)

sample <- spatSample(vec, 7000, method = "random")
bg <- terra::extract(nstack, sample, na.rm = TRUE, xy = TRUE)
bg <- bg[,-c(1)]
bg <- na.omit(bg)
is.na(bg)

# Sample the background according to 3x number of attacks
bg_d <- bg[sample(nrow(bg), 3330),]
bg_i <- bg[sample(nrow(bg), 1041),]
bg_all <- bg[sample(nrow(bg), 4371),]
bg_all$human_death <- NA
bg_all$human_injury <- NA

# Combine attacks and background points
d_d <- rbind(cbind(att = 1, occ_deaths), cbind(att = 0, bg_d))
d_i <- rbind(cbind(att = 1, occ_injuries), cbind(att = 0, bg_i))
d_all <- rbind(cbind(att = 1, occ_all), cbind(att = 0, bg_all))
summary(d_d)
summary(d_i)
summary(d_all)

# Save all
write.csv(d_d, "F:/Spatial Analysis/250617_deaths_datasets.csv")
write.csv(d_i, "F:/Spatial Analysis/250617_injuries_datasets.csv")
write.csv(d_all, "F:/Spatial Analysis/250617_all_datasets.csv")
