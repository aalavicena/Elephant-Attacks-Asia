####################
### Analysis
setwd("E:/Spatial Analysis")

# library(rphylopic)
library(terra)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(gridExtra)
library(conflicted)
library(DHARMa)
library(spaMM)

# Load and scale data
df <- read.csv("250717_all_datasets.csv")
summary(df)

df$att <- as.factor(df$att)
df$dist_forest <- df$dist_forest * 1000 # convert to meter

final.data <- df %>%
  dplyr::select(att,ed,mesh,hfp,dist_pa,gdp,popdens,tri,dist_forest,dist_cropland,dist_fp,c.heights,c.cover,wetness) %>%
  mutate(ed = c(scale(ed)),
         mesh = c(scale(mesh)),
         hfp = c(scale(hfp)),
         dist_pa = c(scale(dist_pa)),
         gdp = c(scale(gdp)),
         popdens = c(scale(popdens)),
         tri = c(scale(tri)),
         dist_forest = c(scale(dist_forest)),
         dist_cropland = c(scale(dist_cropland)),
         dist_fp = c(scale(dist_fp)),
         c.heights = c(scale(c.heights)),
         c.cover = c(scale(c.cover)),
         wetness = c(scale(wetness)))

# Correlation between variables
plot(final.data)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; n8 <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.12345678), digits = digits)[1]
  txt <- paste0(prefix, txt)
  
  if(missing(cex.cor)) 
    cex.cor <- 0.8/strwidth
  text(0.5, 0.5, txt, cex.cor = cex.cor)
}
pairs(final.data, upper.panel = panel.cor, diag.panel = panel.hist)

# Spearman is used when relationship between variables is not linear
final.data %>%
  as.matrix() %>%
  Hmisc::rcorr(., type = "spearman")
rm(panel.hist,panel.cor) # c.cover removed in favour of wetness

final.data <- df %>%
  dplyr::select(att,ed,mesh,hfp,dist_pa,gdp,popdens,tri,dist_forest,dist_cropland,dist_fp,wetness,x,y) %>%
  mutate(ed = c(scale(ed)),
         mesh = c(scale(mesh)),
         hfp = c(scale(hfp)),
         dist_pa = c(scale(dist_pa)),
         gdp = c(scale(gdp)),
         popdens = c(scale(popdens)),
         tri = c(scale(tri)),
         dist_forest = c(scale(dist_forest)),
         dist_cropland = c(scale(dist_cropland)),
         dist_fp = c(scale(dist_fp)),
         wetness = c(scale(wetness)))

weights <- rep(NA, nrow(final.data))
weights[final.data[, "att"] == 1] <- 1
weights[final.data[, "att"] == 0] <- 1456 / 4371 # presence:absence
weights
sum(weights[final.data[ , "att"] == 1])
sum(weights[final.data[ , "att"] == 0])

final.data$weights <- weights # to make sure spaMM works

final.data <- distinct(final.data, x, y, .keep_all = TRUE)
final.data <- na.omit(final.data)

### Modelling
# Parallel processing
nthr <- parallel::detectCores(logical = FALSE) - 1L

## Test each sets of variables
# Fragmentation model
frag.m.w <- fitme(att ~ ed + mesh + Matern(1|y+x),
                data = final.data,
                family = binomial(link = "logit"),
                weights.form = ~weights,
                control.HLfit = list(NbThreads = max(nthr, 1L)))

# Habitat and natural resources model
hbr.m.w <- fitme(att ~ wetness + Matern(1|y+x),
                data = final.data,
                family = binomial(link = "logit"),
                weights.form = ~weights,
                control.HLfit = list(NbThreads = max(nthr, 1L)))

# Socioeconomic model
soc.m.w <- fitme(att ~ popdens + gdp + Matern(1|y+x),
               data = final.data,
               family = binomial(link = "logit"),
               weights.form = ~weights,
               control.HLfit = list(NbThreads = max(nthr, 1L)))

# Potential HEC model
hec.m.w <- fitme(att ~ dist_pa + dist_forest + dist_cropland + dist_fp + Matern(1|y+x),
               data = final.data,
               family = binomial(link = "logit"),
               weights.form = ~weights,
               control.HLfit = list(NbThreads = max(nthr, 1L)))

# Human Disturbance model
hd.m.w <- fitme(att ~ hfp + Matern(1|y+x),
              data = final.data,
              family = binomial(link = "logit"),
              weights.form = ~weights,
              control.HLfit = list(NbThreads = max(nthr, 1L)))

# Potential HEC with quadratic terms
hec_quad.m.w <- fitme(att ~ dist_pa + dist_forest + dist_cropland + dist_fp + I(dist_pa^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + Matern(1|y+x),
               data = final.data,
               family = binomial(link = "logit"),
               weights.form = ~weights,
               control.HLfit = list(NbThreads = max(nthr, 1L)))

# Null model
null.w <- fitme(att ~ 1 + Matern(1|y+x),
              data = final.data,
              family = binomial(link = "logit"),
              weights.form = ~weights,
              control.HLfit = list(NbThreads = max(nthr, 1L)))

# Topography model
topo.m <- fitme(att ~ tri + Matern(1|y+x),
                data = final.data,
                family = binomial(link = "logit"),
                weights.form = ~weights,
                control.HLfit = list(NbThreads = max(nthr, 1L)))

# Full model without quadratic terms
f1.m.w <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp + dist_forest + dist_cropland + dist_fp + dist_pa + Matern(1|y+x), 
            data = final.data,
            family = binomial(link = "logit"),
            weights.form = ~weights,
            control.HLfit = list(NbThreads = max(nthr, 1L)))

# Full model with quadratic terms
f2.m_wetness <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp + dist_forest + dist_cropland + dist_fp + dist_pa + I(wetness^2) + I(dist_forest^2)
                      + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2) + Matern(1|y+x), 
                      data = final.data,
                      family = binomial(link = "logit"),   
                      weights.form = ~weights,
                      control.HLfit = list(NbThreads = max(nthr, 1L)))

# Model without SAC
f3.no_sac <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp + dist_forest + dist_cropland + dist_fp + dist_pa + I(wetness^2) + I(dist_forest^2)
              + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2), 
              data = final.data,
              family = binomial(link = "logit"),   
              weights.form = ~weights)
