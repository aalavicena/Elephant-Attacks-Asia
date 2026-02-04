# =============================================================================
# 4. Results exploration
# Distribution of points, risk proportions, baseline/projected comparisons.
# Set BASE_PATH and paths below to your data locations.
# =============================================================================

BASE_PATH <- "/Volumes/MEC/Spatial Analysis"
PATH_RESULTS <- file.path(BASE_PATH, "_Supplementary Codebase/_Results")
PATH_VARS    <- file.path(BASE_PATH, "_Working Variables")
PATH_BOUNDARY <- file.path(PATH_VARS, "Boundary/IUCN-redlist-100km-dissolved.shp")

library(sf)
library(terra)
library(exactextractr)
library(dplyr)

# Data
df <- read.csv(file.path(BASE_PATH, "250717_all_datasets.csv"))
ras <- rast(file.path(PATH_RESULTS, "baseline-250720.tif"))
elerange <- st_as_sf(PATH_BOUNDARY)
elerange <- st_transform(elerange, st_crs(ras))

# Load polygons
setwd("C:/Users/lenovo/OneDrive/Research/_PhD/2 - Wild Elephant Attacks on People/Data/GIS/National Boundaries")

bgd <- read_sf("Bangladesh/gadm41_BGD_shp/gadm41_BGD_0.shp")
ind <- read_sf("India/gadm41_IND_shp/gadm41_IND_0.shp")
mmr <- read_sf("Myanmar/gadm41_MMR_shp/gadm41_MMR_0.shp")
slk <- read_sf("Sri Lanka/gadm41_LKA_shp/gadm41_LKA_0.shp")
chn <- read_sf("China/NAME_1_Yunnan.shp")
tha <- read_sf("Thailand/gadm41_THA_shp/gadm41_THA_0.shp")
lao <- read_sf("Laos/gadm41_LAO_shp/gadm41_LAO_0.shp")
mas_p <- read_sf("Malaysia/MYS (Peninsular)_GADM_0.shp")
mas_s <- read_sf("Malaysia/MYS (Borneo)_GADM.shp")
ina <- read_sf("Indonesia/IDN_Sumatra_GADM (Level 0).shp")
btn <- read_sf("Bhutan/gadm41_BTN_shp/gadm41_BTN_0.shp")
cam <- read_sf("Cambodia/gadm41_KHM_shp/gadm41_KHM_0.shp")
vnm <- read_sf("Vietnam/gadm41_VNM_shp/gadm41_VNM_0.shp")
npl <- read_sf("Nepal/gadm41_NPL_shp/gadm41_NPL_0.shp") # change to 0 for country level (distribution of data points), 1 for province level (benchmark comparison)

countries <- list(bgd,ind,mmr,slk,chn,tha,lao,mas_p,mas_s,ina,btn,cam,vnm,npl)
countries <- bind_rows(countries)
countries <- st_transform(countries, crs = st_crs(ras))

### Exploration
# 1. Distribution of data points
df.sf <- df %>% 
  select(att, x, y, human_death, human_injury) %>%
  filter(att == 1)

df.sf <- df.sf[-1437,]
df.sf <- st_as_sf(df.sf, coords = c("x", "y"), crs = 4326)

sum.df <- countries %>% 
  mutate(counts = lengths(st_intersects(., df.sf)))

summary <- sum.df %>%
  group_by(COUNTRY) %>%
  summarise(count = sum(counts))

# 2. In relation to IUCN Range Map
iucn <- read_sf("C:/Users/lenovo/OneDrive/Research/_PhD/2 - Wild Elephant Attacks on People/Data/GIS/Asian Elephant Distribution/IUCN_redlist.shp")
plot(iucn)

sum.iucn <- iucn %>%
  mutate(counts = lengths(st_intersects(., df.sf)))

# 3. Predicted risk and benchmark data comparison
ext <- exact_extract(ras, countries, fun = "weighted_mean", weights = "area")
print(ext)

ext <- as.data.frame(ext)
ext <- cbind(ext, countries$NAME_1)
colnames(ext) <- c("weighted_mean_risk", "province")
write.csv(ext, "E:/Spatial Analysis/_Supplementary Codebase/_Results/weighted-mean-risk.csv")

bm <- read.csv("E:/Spatial Analysis/_Supplementary Codebase/_Results/risk-corr.csv")
summary(bm)

cor.test(bm$weighted_mean_risk, bm$crude_standard, method = 'kendall')

# 4. Baseline risk proportion and distribution
bgd <- vect("National Boundaries/gadm41_BGD_0.shp")
ind <- vect("National Boundaries/gadm41_IND_0.shp")
mmr <- vect("National Boundaries/gadm41_MMR_0.shp")
slk <- vect("National Boundaries/gadm41_LKA_0.shp")
chn <- vect("National Boundaries/Yunnan.shp")
tha <- vect("National Boundaries/gadm41_THA_0.shp")
lao <- vect("National Boundaries/gadm41_LAO_0.shp")
mas_p <- vect("National Boundaries/peninsula-malaysia.shp")
mas_s <- vect("National Boundaries/sabah.shp")
ina <- vect("National Boundaries/idn_sumatra.shp")
btn <- vect("National Boundaries/gadm41_BTN_0.shp")
cam <- vect("National Boundaries/gadm41_KHM_0.shp")
vnm <- vect("National Boundaries/gadm41_VNM_0.shp")
npl <- vect("National Boundaries/gadm41_NPL_0.shp")
countries <- list(bgd,ind,mmr,slk,chn,tha,lao,mas_p,mas_s,ina,btn,cam,vnm,npl)


# Categorise raster and polygonize to reduce commission error
m <- c(-Inf, 0.4, 1,
       0.41, 0.59, 2,
       0.6, 0.79, 3,
       0.8, Inf, 4)

m <- matrix(m, ncol = 3, byrow = TRUE)

ras <- classify(ras, m, others = NA)
levels(ras) <- data.frame(ID = c(1, 2, 3, 4), 
                          Category = c("Low", "Moderate", "High", 
                                       "Severe"))
plot(ras)

at.risk <- as.polygons(ras)
countries <- lapply(countries, project, at.risk)

rbgd <- expanse(crop(at.risk, countries[[1]]), unit = 'km')
rind <- expanse(crop(at.risk, countries[[2]]), unit = 'km')
rmmr <- expanse(crop(at.risk, countries[[3]]), unit = 'km')
rslk <- expanse(crop(at.risk, countries[[4]]), unit = 'km')
rchn <- expanse(crop(at.risk, countries[[5]]), unit = 'km')
rtha <- expanse(crop(at.risk, countries[[6]]), unit = 'km')
rlao <- expanse(crop(at.risk, countries[[7]]), unit = 'km')
rmas <- expanse(crop(at.risk, countries[[8]]), unit = 'km')
rmas_s <- expanse(crop(at.risk, countries[[9]]), unit = 'km')
rina <- expanse(crop(at.risk, countries[[10]]), unit = 'km')
rbtn <- expanse(crop(at.risk, countries[[11]]), unit = 'km')
rcam <- expanse(crop(at.risk, countries[[12]]), unit = 'km')
rvnm <- expanse(crop(at.risk, countries[[13]]), unit = 'km')
rnpl <- expanse(crop(at.risk, countries[[14]]), unit = 'km')

risk.prop <- list(bgd=rbgd,ind=rind,mmr=rmmr,slk=rslk,chn=rchn,tha=rtha,lao=rlao,mas=rmas,mas_s=rmas_s,ina=rina,btn=rbtn,cam=rcam,vnm=rvnm,npl=rnpl)
max_len <- max(sapply(risk.prop, length))
pad <- lapply(risk.prop, function(v) {
  c(v, rep(NA, max_len - length(v)))
})
risk.prop <- as.data.frame(pad)
risk.prop <- as.data.frame(t(risk.prop))

#risk.prop <- risk.prop %>%
 # group_by(value) %>%
  #mutate(row_id = row_number()) %>%
  #ungroup() %>%
  #tidyr::pivot_wider(names_from = value,
   #           values_from = count,
    #          id_cols = row_id)

risk.prop$country <- c("Bangladesh", "India", "Myanmar", "Sri Lanka", "China", "Thailand", "Laos", "Peninsula Malaysia", "Sabah (Malaysia)", "Indonesia", 
             "Bhutan", "Cambodia", "Vietnam", "Nepal")
colnames(risk.prop) <- c('Low', 'Moderate', 'High', 'Severe', 'Country')
risk.prop <- relocate(risk.prop, Country, .before = Low)
rownames(risk.prop) <- NULL

risk.prop$total_area_at_risk <- rowSums(risk.prop[,2:5], na.rm = TRUE)

# Percentages
categories <- c("Low", "Moderate", "High", "Severe")

risk.prop$total_area <- risk.prop$Extended_Buffer_Area
for (cat in categories) {
  risk.prop[, paste0(cat, "_pct")] <- (risk.prop[, cat] / risk.prop$total_area) * 100
}

risk.prop$total_area_pct <- (risk.prop$total_area_at_risk / sum(risk.prop$total_area_at_risk)) * 100

# Area of extended buffer
bbgd <- expanse(crop(elerange, countries[[1]]), unit = 'km')
bind <- expanse(crop(elerange, countries[[2]]), unit = 'km')
bmmr <- expanse(crop(elerange, countries[[3]]), unit = 'km')
bslk <- expanse(crop(elerange, countries[[4]]), unit = 'km')
bchn <- expanse(crop(elerange, countries[[5]]), unit = 'km')
btha <- expanse(crop(elerange, countries[[6]]), unit = 'km')
blao <- expanse(crop(elerange, countries[[7]]), unit = 'km')
bmas <- expanse(crop(elerange, countries[[8]]), unit = 'km')
bmas_s <- expanse(crop(elerange, countries[[9]]), unit = 'km')
bina <- expanse(crop(elerange, countries[[10]]), unit = 'km')
bbtn <- expanse(crop(elerange, countries[[11]]), unit = 'km')
bcam <- expanse(crop(elerange, countries[[12]]), unit = 'km')
bvnm <- expanse(crop(elerange, countries[[13]]), unit = 'km')
bnpl <- expanse(crop(elerange, countries[[14]]), unit = 'km')

b.countries <- list(bgd=bbgd,ind=bind,mmr=bmmr,slk=bslk,chn=bchn,tha=btha,lao=blao,mas=bmas,mas_s=bmas_s,ina=bina,btn=bbtn,cam=bcam,vnm=bvnm,npl=bnpl)
b.countries <- as.data.frame(b.countries)
b.countries <- as.data.frame(t(b.countries))
colnames(b.countries) <- "Extended_Buffer_Area"
rownames(b.countries) <- NULL

risk.prop <- cbind(risk.prop, b.countries)
risk.prop <- relocate(risk.prop, Extended_Buffer_Area, .before = Low)
risk.prop$total_area <- risk.prop$Extended_Buffer_Area
risk.prop$pct_at_risk <- (risk.prop$total_area_at_risk / risk.prop$total_area) * 100

write.csv(risk.prop, file.path(PATH_RESULTS, "baseline-risk-proportion.csv"))

# 6. Baseline human population at risk
path <- "/Volumes/MEC/Spatial Analysis/Datasets/WorldPop/counts_2015/"
files <- list.files(path = path, pattern = "\\.tif$", full.names = TRUE)
r <- lapply(files, rast)
r <- sprc(r)
r <- merge(r)

bgd <- st_read("National Boundaries/gadm41_BGD_0.shp", crs = 4326)
ind <- st_read("National Boundaries/gadm41_IND_0.shp", crs = 4326)
mmr <- st_read("National Boundaries/gadm41_MMR_0.shp", crs = 4326)
slk <- st_read("National Boundaries/gadm41_LKA_0.shp", crs = 4326)
chn <- st_read("National Boundaries/Yunnan.shp", crs = 4326)
tha <- st_read("National Boundaries/gadm41_THA_0.shp", crs = 4326)
lao <- st_read("National Boundaries/gadm41_LAO_0.shp", crs = 4326)
mas_p <- st_read("National Boundaries/peninsula-malaysia.shp", crs = 4326)
mas_s <- st_read("National Boundaries/sabah.shp", crs = 4326)
ina <- st_read("National Boundaries/idn_sumatra.shp", crs = 4326)
btn <- st_read("National Boundaries/gadm41_BTN_0.shp", crs = 4326)
cam <- st_read("National Boundaries/gadm41_KHM_0.shp", crs = 4326)
vnm <- st_read("National Boundaries/gadm41_VNM_0.shp", crs = 4326)
npl <- st_read("National Boundaries/gadm41_NPL_0.shp", crs = 4326)
countries <- list(bgd,ind,mmr,slk,chn,tha,lao,mas_p,mas_s,ina,btn,cam,vnm,npl)

rbgd <- crop(at.risk, bgd)
rind <- crop(at.risk, ind)
rmmr <- crop(at.risk, mmr)
rslk <- crop(at.risk, slk)
rchn <- crop(at.risk, chn)
rtha <- crop(at.risk, tha)
rlao <- crop(at.risk, lao)
rmas <- crop(at.risk, mas_p)
rmas_s <- crop(at.risk, mas_s)
rina <- crop(at.risk, ina)
rbtn <- crop(at.risk, btn)
rnpl <- crop(at.risk, npl)
rvnm <- crop(at.risk, vnm)
rkhm <- crop(at.risk, cam)

rbgd <- st_as_sf(rbgd)
rind <- st_as_sf(rind)
rmmr <- st_as_sf(rmmr)
rslk <- st_as_sf(rslk)
rchn <- st_as_sf(rchn)
rtha <- st_as_sf(rtha)
rlao <- st_as_sf(rlao)
rmas <- st_as_sf(rmas)
rmas_s <- st_as_sf(rmas_s)
rina <- st_as_sf(rina)
rbtn <- st_as_sf(rbtn)
rnpl <- st_as_sf(rnpl)
rvnm <- st_as_sf(rvnm)
rkhm <- st_as_sf(rkhm)

zbgd <- exact_extract(r, rbgd, fun = "sum", weights = "area")
zind <- exact_extract(r, rind, fun = "sum", weights = "area")
zmmr <- exact_extract(r, rmmr, fun = "sum", weights = "area")
zslk <- exact_extract(r, rslk, fun = "sum", weights = "area")
zchn <- exact_extract(r, rchn, fun = "sum", weights = "area")
ztha <- exact_extract(r, rtha, fun = "sum", weights = "area")
zlao <- exact_extract(r, rlao, fun = "sum", weights = "area")
zmas <- exact_extract(r, rmas, fun = "sum", weights = "area")
zmas_s <- exact_extract(r, rmas_s, fun = "sum", weights = "area")
zina <- exact_extract(r, rina, fun = "sum", weights = "area")
zbtn <- exact_extract(r, rbtn, fun = "sum", weights = "area")
znpl <- exact_extract(r, rnpl, fun = "sum", weights = "area")
zvnm <- exact_extract(r, rvnm, fun = "sum", weights = "area")
zkhm <- exact_extract(r, rkhm, fun = "sum", weights = "area")

risk.counts <- list(bgd=zbgd,ind=zind,mmr=zmmr,slk=zslk,chn=zchn,tha=ztha,lao=zlao,mas=zmas,mas_s=zmas_s,
                  ina=zina,btn=zbtn,cam=zkhm,vnm=zvnm,npl=znpl)
max_len <- max(sapply(risk.counts, length))
pad <- lapply(risk.counts, function(v) {
  c(v, rep(NA, max_len - length(v)))
})
risk.counts <- as.data.frame(pad)
risk.counts <- as.data.frame(t(risk.counts))

write.csv(risk.counts, "_Supplementary Codebase/_Results/baseline-popcounts-at-risk.csv")

# 7. Projected Risks
bc <- rast("_Supplementary Codebase/_Results/baseline-250720.tif")
ssp1 <- rast("_Supplementary Codebase/_Results/ssp1-250720.tif")
ssp3 <- rast("_Supplementary Codebase/_Results/ssp3-250720.tif")
ssp5 <- rast("_Supplementary Codebase/_Results/ssp5-250720.tif")

# Calculate differences
sc1 <- ssp1 - bc
sc3 <- ssp3 - bc
sc5 <- ssp5 - bc

# Average Change
avg <- (sc1 + sc3 + sc5) / 3
plot(avg)

# use SD to detect changes
m <- c(-0.5, -0.06, 1,
       -0.07, -0.02, 2,  
       -0.03, 0.03, 3,
       0.04, 0.08, 4,
       0.09, 0.5, 5)

m <- matrix(m, ncol = 3, byrow = TRUE)

avg <- classify(avg, m, others = NA)
plot(avg)

levels(avg) <- data.frame(ID = c(1,2,3,4,5),
                          Category = c("Strong Decrease", "Moderate Decrease", "No Change", "Moderate Increase", "Strong Increase"))
plot(avg)
writeRaster(avg, "_Supplementary Codebase/_Results/ssp-avg-250729.tif", overwrite = TRUE)

avg.b <- as.polygons(avg)
ch <- as.data.frame(expanse(avg.b, unit = 'km'))
colnames(ch) <- 'area'
rownames(ch) <- c('Strong Decrease', 'Moderate Decrease', 'No Change', 'Moderate Increase', 'Strong Increase')
ch$pct <- (ch[,1]/sum(ch[,1]))*100
ch

# Changes per scenario
# SSP1
ch1 <- classify(sc1, m, others = NA)
plot(ch1)

levels(ch1) <- data.frame(ID = c(1,2,3,4,5),
                          Category = c("Strong Decrease", "Moderate Decrease", "No Change", "Moderate Increase", "Strong Increase"))
plot(ch1)
writeRaster(avg, "_Supplementary Codebase/_Results/ssp1-changes-250729.tif", overwrite = TRUE)

ssp1.ch <- freq(ch1)
ssp1.ch$pct <- (ssp1.ch[,3]/sum(ssp1.ch[,3]))*100
ssp1.ch

# SSP3
ch3 <- classify(sc3, m, others = NA)
plot(ch3)

levels(ch3) <- data.frame(ID = c(1,2,3,4,5),
                          Category = c("Strong Decrease", "Moderate Decrease", "No Change", "Moderate Increase", "Strong Increase"))
plot(ch3)
writeRaster(avg, "_Supplementary Codebase/_Results/ssp3-changes-250729.tif")

ssp3.ch <- freq(ch3)
ssp3.ch$pct <- (ssp3.ch[,3]/sum(ssp3.ch[,3]))*100
ssp3.ch

# SSP5
ch5 <- classify(sc5, m, others = NA)
plot(ch5)

levels(ch5) <- data.frame(ID = c(1,2,3,4,5),
                          Category = c("Strong Decrease", "Moderate Decrease", "No Change", "Moderate Increase", "Strong Increase"))
plot(ch5)
writeRaster(avg, "_Supplementary Codebase/_Results/ssp5-changes-250729.tif")

ssp5.ch <- freq(ch5)
ssp5.ch$pct <- (ssp5.ch[,3]/sum(ssp5.ch[,3]))*100
ssp5.ch

# 8. Projected Human Population at Risk
p1 <- rast("/Volumes/MEC/Spatial Analysis/Datasets/projected-popcounts/SSP1_2050.tif")
p3 <- rast("/Volumes/MEC/Spatial Analysis/Datasets/projected-popcounts/SSP3_2050.tif")
p5 <- rast("/Volumes/MEC/Spatial Analysis/Datasets/projected-popcounts/SSP5_2050.tif")
elerange <- vect("/Volumes/MEC/Spatial Analysis/_Working Variables/Boundary/tropical-asia.gpkg")
elerange <- project(elerange, p1)
p1 <- crop(p1, elerange, mask = TRUE)
p3 <- crop(p3, elerange, mask = TRUE)
p5 <- crop(p5, elerange, mask = TRUE)

pavg <- (p1 + p3 + p5) / 3

r.avg <- (ssp1 + ssp3 + ssp5) / 3
m <- c(-Inf, 0.4, 1,
       0.41, 0.59, 2,
       0.6, 0.79, 3,
       0.8, Inf, 4)
m <- matrix(m, ncol = 3, byrow = TRUE)
r.avg <- classify(r.avg, m, others = NA)
levels(r.avg) <- data.frame(ID = c(1, 2, 3, 4), 
                          Category = c("Low", "Moderate", "High", 
                                       "Severe"))
plot(r.avg)

r.avg <- as.polygons(r.avg)

bgd <- st_read("National Boundaries/gadm41_BGD_0.shp", crs = 4326)
ind <- st_read("National Boundaries/gadm41_IND_0.shp", crs = 4326)
mmr <- st_read("National Boundaries/gadm41_MMR_0.shp", crs = 4326)
slk <- st_read("National Boundaries/gadm41_LKA_0.shp", crs = 4326)
chn <- st_read("National Boundaries/Yunnan.shp", crs = 4326)
tha <- st_read("National Boundaries/gadm41_THA_0.shp", crs = 4326)
lao <- st_read("National Boundaries/gadm41_LAO_0.shp", crs = 4326)
mas_p <- st_read("National Boundaries/peninsula-malaysia.shp", crs = 4326)
mas_s <- st_read("National Boundaries/sabah.shp", crs = 4326)
ina <- st_read("National Boundaries/idn_sumatra.shp", crs = 4326)
btn <- st_read("National Boundaries/gadm41_BTN_0.shp", crs = 4326)
cam <- st_read("National Boundaries/gadm41_KHM_0.shp", crs = 4326)
vnm <- st_read("National Boundaries/gadm41_VNM_0.shp", crs = 4326)
npl <- st_read("National Boundaries/gadm41_NPL_0.shp", crs = 4326)

rbgd <- crop(r.avg, bgd)
rind <- crop(r.avg, ind)
rmmr <- crop(r.avg, mmr)
rslk <- crop(r.avg, slk)
rchn <- crop(r.avg, chn)
rtha <- crop(r.avg, tha)
rlao <- crop(r.avg, lao)
rmas <- crop(r.avg, mas_p)
rmas_s <- crop(r.avg, mas_s)
rina <- crop(r.avg, ina)
rbtn <- crop(r.avg, btn)
rnpl <- crop(r.avg, npl)
rvnm <- crop(r.avg, vnm)
rkhm <- crop(r.avg, cam)

rbgd <- st_as_sf(rbgd)
rind <- st_as_sf(rind)
rmmr <- st_as_sf(rmmr)
rslk <- st_as_sf(rslk)
rchn <- st_as_sf(rchn)
rtha <- st_as_sf(rtha)
rlao <- st_as_sf(rlao)
rmas <- st_as_sf(rmas)
rmas_s <- st_as_sf(rmas_s)
rina <- st_as_sf(rina)
rbtn <- st_as_sf(rbtn)
rnpl <- st_as_sf(rnpl)
rvnm <- st_as_sf(rvnm)
rkhm <- st_as_sf(rkhm)

zbgd <- exact_extract(pavg, rbgd, fun = "sum", weights = "area")
zind <- exact_extract(pavg, rind, fun = "sum", weights = "area")
zmmr <- exact_extract(pavg, rmmr, fun = "sum", weights = "area")
zslk <- exact_extract(pavg, rslk, fun = "sum", weights = "area")
zchn <- exact_extract(pavg, rchn, fun = "sum", weights = "area")
ztha <- exact_extract(pavg, rtha, fun = "sum", weights = "area")
zlao <- exact_extract(pavg, rlao, fun = "sum", weights = "area")
zmas <- exact_extract(pavg, rmas, fun = "sum", weights = "area")
zmas_s <- exact_extract(pavg, rmas_s, fun = "sum", weights = "area")
zina <- exact_extract(pavg, rina, fun = "sum", weights = "area")
zbtn <- exact_extract(pavg, rbtn, fun = "sum", weights = "area")
znpl <- exact_extract(pavg, rnpl, fun = "sum", weights = "area")
zvnm <- exact_extract(pavg, rvnm, fun = "sum", weights = "area")
zkhm <- exact_extract(pavg, rkhm, fun = "sum", weights = "area")

risk.counts <- list(bgd=zbgd,ind=zind,mmr=zmmr,slk=zslk,chn=zchn,tha=ztha,lao=zlao,mas=zmas,mas_s=zmas_s,
                    ina=zina,btn=zbtn,cam=zkhm,vnm=zvnm,npl=znpl)

max_len <- max(sapply(risk.counts, length))
pad <- lapply(risk.counts, function(v) {
  c(v, rep(NA, max_len - length(v)))
})
risk.counts <- as.data.frame(pad)
risk.counts <- as.data.frame(t(risk.counts))
colnames(risk.counts) <- c("Low", "Moderate", "High", "Severe")
risk.counts

write.csv(risk.counts, "_Supplementary Codebase/_Results/ssp-avg-popcounts-at-risk.csv")
