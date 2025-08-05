## Model Selection and Inference
setwd("E:/Spatial Analysis/_Supplementary Codebase")

library(DHARMa)
library(spaMM)
library(dplyr)
library(ggplot2)

f2.m <- readRDS("_Models/sets 3/f2.m_wetness.rds")
f2.m.no_quad <- readRDS("_Models/sets 2/f2.wetness_noquad.rds")
extractAIC(f2.m)
extractAIC(f2.m.no_quad)
simulateResiduals(f2.m.no_quad, n = 1000, plot = TRUE)

summary(f2.m)

path <- "E:/Spatial Analysis/_Supplementary Codebase/_Models"
files <- list.files(path = path, full.names = TRUE)

models <- lapply(files, readRDS)
file_names <- tools::file_path_sans_ext(basename(files))
names(models) <- file_names
print(names(models))

aic <- lapply(models, extractAIC)
aic <- do.call(rbind, aic)
aic <- as.data.frame(aic)
aic$model <- rownames(aic) 
rownames(aic) <- NULL

ranked_aic <- aic[order(aic$AIC), ]

# calculate delta AIC
min <- ranked_aic$AIC[1]
ranked_aic$delta_aic <- ranked_aic$AIC - min 

print(ranked_aic) # f2.m performs best

# Check assumptions
f2.m <- readRDS("E:/Spatial Analysis/_Supplementary Codebase/_Models/f2.m.rds")
res <- simulateResiduals(f2.m, n = 1000)
plot(res)

pseudoR2(f2.m) # 0.35 P.S. Not so important, only provide if reviewers asked! -spaMM authors

# SAC
nbc <- 75
cor_f2 <- pgirmess::correlog(coords = final.data[,c("x", "y")],
                          z = residuals(f2.m),
                          method = "Moran", nbclass = nbc)
cor_f2 <- as.data.frame(cor_f2)

# compare with structured model
cor_f3 <- pgirmess::correlog(coords = final.data[,c("x", "y")],
                             z = residuals(f3.no_sac),
                             method = "Moran", nbclass = nbc)
cor_f3 <- as.data.frame(cor_f3)

cor_f2$variable <- "Spatial Model"
cor_f3$variable <- "Non-Spatial Model"
cor <- rbind(cor_f2, cor_f3)

correlogram <- ggplot(cor, aes(x = dist.class, y = coef, color = variable, shape = variable)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_line(linewidth = 1) + 
  
  geom_point(size = 3, fill = "white") +
  
  scale_color_manual(values = c("Non-Spatial Model" = "#d73027", "Spatial Model" = "#4575b4")) +
  scale_shape_manual(values = c("Non-Spatial Model" = 16, "Spatial Model" = 17)) +
  scale_fill_manual(values = c("Non-Spatial Model" = "#d73027", "Spatial Model" = "#4575b4")) +
  
  labs(
    title = "Residual Spatial Autocorrelation",
    x = "Distance Class",
    y = "Moran's I Coefficient",
    color = "Model Type",
    shape = "Model Type"
  ) +
  
  theme_bw() +
  theme(legend.position = "bottom", legend.title.align = 0.5)

ggsave("_Plots/correlogram-250721.jpeg",
       correlogram,
       jpeg,
       width = 15,
       height = 7,
       units = "in",
       dpi = 300)

### Inference
# Odds Ratio and 95% Confidence Interval
coefs <- as.data.frame(summary(f2.m)$beta_table)
coefs

var.list <- c("ed", "mesh", "tri", "wetness", "hfp", "popdens", "gdp",
              "dist_forest", "dist_cropland", "dist_fp", "dist_pa")

row <- row.names(coefs) %in% var.list
lower <- coefs[row, 'Estimate'] - 1.96*coefs[row, 'Cond. SE']
upper <- coefs[row, 'Estimate'] + 1.96*coefs[row, 'Cond. SE']
c(lower, upper)

odds <- exp(coefs$Estimate)

odds.ci <- data.frame(variable = var.list,
                      odds = odds[2:12],
                      lower.odds = exp(lower),
                      upper.odds = exp(upper),
                      lower.ci = lower,
                      upper.ci = upper)
gt::gt(odds.ci)


# ------------------------------------------------------------------------------
## Plot Odds Ratio
var.names <- c("Edge Density", "Effective Mesh Size", "Terrain Ruggedness Index","Tassled Cap Wetness Index",
               "Human Footprint", "Human Population Density", "Gross Domestic Product",
               "Distance to Forest", "Distance to Cropland", "Distance to Forest Plantations",
               "Distance to Protected Areas")

# Dataframe for graphical purposes
odds.df <- data.frame(yAxis = length(var.names):1,
                      odds = odds[2:12],
                      lower = exp(lower),
                      upper = exp(upper),
                      row.names = var.names)
odds.df

# Define variable sets
frag <- c("Edge Density", "Effective Mesh Size")
soc <- c("Human Population Density", "Gross Domestic Product")
hab <- "Tassled Cap Wetness Index"
hec <- c("Distance to Forest", "Distance to Cropland", "Distance to Forest Plantations",
         "Distance to Protected Areas")
hd <- "Human Footprint"
tcm <- "Terrain Ruggedness Index"

# Create new column to categorize variables
odds.df <- odds.df %>%
  mutate(
    variable_set = case_when(
      rownames(odds.df) %in% frag ~ "Fragmentation",
      rownames(odds.df) %in% soc ~ "Socioeconomic",
      rownames(odds.df) %in% hab ~ "Habitat and Natural Resources",
      rownames(odds.df) %in% hec ~ "Potential HEC",
      rownames(odds.df) %in% hd ~ "Human Disturbance",
      rownames(odds.df) %in% tcm ~ "Topography",
      TRUE ~ "Other"
    )
  )

odds.df$variable_set <- factor(odds.df$variable_set, levels = c("Fragmentation", "Socioeconomic", "Habitat and Natural Resources", 
                                                                "Potential HEC", "Human Disturbance", "Topography"))
odds.df

odds.plot <- 
  ggplot(data = odds.df,
         aes(x = odds,
             y = reorder(rownames(odds.df), yAxis))) +
  geom_vline(
    aes(xintercept = 1),
    linewidth = 1,
    linetype = 'dashed') +
  geom_errorbarh(
    aes(xmax = upper,
        xmin = lower),
    linewidth = 1,
    height = .4,
    color = 'gray50') +
  geom_point(size = 3.5) +
  coord_cartesian(xlim = c(0, 3)) +
  scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  facet_wrap(~ variable_set, scales = "free_y", ncol = 1) + 
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 28, face="bold", hjust = 0), 
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 22, hjust = 1), 
    axis.title.x = element_text(size = 28),
    panel.spacing = unit(1, "lines") 
  ) +
  ylab('') +
  xlab('Odds Ratios')

odds.plot
ggsave('E:/Spatial Analysis/_Supplementary Codebase/_Plots/250720-odds_plot.jpeg',
       odds.plot,
       jpeg,
       width = 15,
       height = 16,
       units = 'in',
       dpi = 300)


# ------------------------------------------------------------------------------
## Plot Effect Sizes

# 1. Edge Density
ed_scaled <- seq(min(final.data$ed),
                 max(final.data$ed),
                 by = 0.01)

ed_new <- data.frame(ed = rep(ed_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_ed <- predict(f2.m, newdata = ed_new, re.form = NA, type = "response")
att_ed <- as.data.frame(att_ed)

range(df$ed, na.rm = TRUE)
ed_raw <- seq(from = 0, to = 10, length.out = 482) # raw values

ed_new <- cbind(ed_raw, ed_new)

ed_graph <- 
  ggplot() +
  geom_line(data = ed_new, 
            aes(x = ed_raw, 
                y = att_ed$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,10)) +
  labs(x = "Edge Density",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(ed_graph)

# 2. Effective Mesh Size
mesh_scaled <- seq(min(final.data$mesh),
                 max(final.data$mesh),
                 by = 0.01)

mesh_new <- data.frame(mesh = rep(mesh_scaled),
                     hfp = mean(final.data$hfp),
                     ed = mean(final.data$ed),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_mesh <- predict(f2.m, newdata = mesh_new, re.form = NA, type = "response")
att_mesh <- as.data.frame(att_mesh)

range(df$mesh, na.rm = TRUE)
mesh_raw <- seq(from = 100, to = 12100, length.out = 319) # raw values

mesh_new <- cbind(mesh_raw, mesh_new)

mesh_graph <- 
  ggplot() +
  geom_line(data = mesh_new, 
            aes(x = mesh_raw, 
                y = att_mesh$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,12100)) +
  labs(x = "Effective Mesh Size",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(mesh_graph)

# 3. Terrain Ruggedness Index
tri_scaled <- seq(min(final.data$tri),
                 max(final.data$tri),
                 by = 0.01)

tri_new <- data.frame(tri = rep(tri_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     ed = mean(final.data$ed),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_tri <- predict(f2.m, newdata = tri_new, re.form = NA, type = "response")
att_tri <- as.data.frame(att_tri)

range(df$tri, na.rm = TRUE)
tri_raw <- seq(from = 0, to = 1745, length.out = 859) # raw values

tri_new <- cbind(tri_raw, tri_new)

tri_graph <- 
  ggplot() +
  geom_line(data = tri_new, 
            aes(x = tri_raw, 
                y = att_tri$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,1745)) +
  labs(x = "Terrain Ruggedness Index",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(tri_graph)

# 5. Wetness
w_scaled <- seq(min(final.data$wetness),
                  max(final.data$wetness),
                  by = 0.01)

w_new <- data.frame(wetness = rep(w_scaled),
                      ed = mean(final.data$ed),
                      mesh = mean(final.data$mesh),
                      tri = mean(final.data$tri),
                      gdp = mean(final.data$gdp),
                      dist_forest = mean(final.data$dist_forest),
                      dist_pa = mean(final.data$dist_pa),
                      dist_cropland = mean(final.data$dist_cropland),
                      dist_fp = mean(final.data$dist_fp),
                    hfp = mean(final.data$hfp),
                    popdens = mean(final.data$popdens))

att_w <- predict(f2.m, newdata = w_new, re.form = NA, type = "response")
att_w <- as.data.frame(att_w)

range(df$wetness, na.rm = TRUE)
w_raw <- seq(from = -0.26, to = 0.49, length.out = 1463) # raw values

w_new <- cbind(w_raw, w_new)

w_graph <- 
  ggplot() +
  geom_line(data = w_new, 
            aes(x = w_raw, 
                y = att_w$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(-0.26, 0.5)) +
  labs(x = "Tassled Cap Wetness Index",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(w_graph)


# 8. Human Footprint
hfp_scaled <- seq(min(final.data$hfp),
                 max(final.data$hfp),
                 by = 0.01)

hfp_new <- data.frame(hfp = rep(hfp_scaled),
                     ed = mean(final.data$ed),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_hfp <- predict(f2.m, newdata = hfp_new, re.form = NA, type = "response")
att_hfp <- as.data.frame(att_hfp)

range(df$hfp, na.rm = TRUE)
hfp_raw <- seq(from = 0, to = 128, length.out = 1530) # raw values

hfp_new <- cbind(hfp_raw, hfp_new)

hfp_graph <- 
  ggplot() +
  geom_line(data = hfp_new, 
            aes(x = hfp_raw, 
                y = att_hfp$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,128)) +
  labs(x = "Human Footprint",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(hfp_graph)

# 9. Population Density
popdens_scaled <- seq(min(final.data$popdens),
                 max(final.data$popdens),
                 by = 0.01)

popdens_new <- data.frame(popdens = rep(popdens_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     ed = mean(final.data$ed),
                     wetness = mean(final.data$wetness))

att_popdens <- predict(f2.m, newdata = popdens_new, re.form = NA, type = "response")
att_popdens <- as.data.frame(att_popdens)

range(df$popdens, na.rm = TRUE)
popdens_raw <- seq(from = 0, to = 32265, length.out = 4315) # raw values

popdens_new <- cbind(popdens_raw, popdens_new)

popdens_graph <- 
  ggplot() +
  geom_line(data = popdens_new, 
            aes(x = popdens_raw, 
                y = att_popdens$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,32265)) +
  labs(x = "Human Population Density",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(popdens_graph)

# 10. Gross Domestic Product
gdp_scaled <- seq(min(final.data$gdp),
                 max(final.data$gdp),
                 by = 0.01)

gdp_new <- data.frame(gdp = rep(gdp_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     ed = mean(final.data$ed),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_gdp <- predict(f2.m, newdata = gdp_new, re.form = NA, type = "response")
att_gdp <- as.data.frame(att_gdp)

range(df$gdp, na.rm = TRUE)
gdp_raw <- seq(from = 0, to = 87, length.out = 1668) # raw values

gdp_new <- cbind(gdp_raw, gdp_new)

gdp_graph <- 
  ggplot() +
  geom_line(data = gdp_new, 
            aes(x = gdp_raw, 
                y = att_gdp$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,87)) +
  labs(x = "Gross Domestic Product",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(gdp_graph)

# 11. Distance to Forest
forest_scaled <- seq(min(final.data$dist_forest),
                 max(final.data$dist_forest),
                 by = 0.01)

forest_new <- data.frame(dist_forest = rep(forest_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     ed = mean(final.data$ed),
                     dist_pa = mean(final.data$dist_pa),
                     dist_cropland = mean(final.data$dist_cropland),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_forest <- predict(f2.m, newdata = forest_new, re.form = NA, type = "response")
att_forest <- as.data.frame(att_forest)

range(df$dist_forest, na.rm = TRUE)
forest_raw <- seq(from = 0, to = 645, length.out = 953) # raw values

forest_new <- cbind(forest_raw, forest_new)

forest_graph <- 
  ggplot() +
  geom_line(data = forest_new, 
            aes(x = forest_raw, 
                y = att_forest$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,645)) +
  labs(x = "Distance to Forest (m)",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(forest_graph)

# 12. Distance to Cropland
crp_scaled <- seq(min(final.data$dist_cropland),
                 max(final.data$dist_cropland),
                 by = 0.01)

crp_new <- data.frame(dist_cropland = rep(crp_scaled),
                     hfp = mean(final.data$hfp),
                     mesh = mean(final.data$mesh),
                     tri = mean(final.data$tri),
                     gdp = mean(final.data$gdp),
                     dist_forest = mean(final.data$dist_forest),
                     dist_pa = mean(final.data$dist_pa),
                     ed = mean(final.data$ed),
                     dist_fp = mean(final.data$dist_fp),
                     popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_crp <- predict(f2.m, newdata = crp_new, re.form = NA, type = "response")
att_crp <- as.data.frame(att_crp)

range(df$dist_cropland, na.rm = TRUE)
crp_raw <- seq(from = 0, to = 204185, length.out = 1276) # raw values

crp_new <- cbind(crp_raw, crp_new)

crp_graph <- 
  ggplot() +
  geom_line(data = crp_new, 
            aes(x = crp_raw, 
                y = att_crp$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,204185)) +
  labs(x = "Distance to Cropland (m)",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(crp_graph)

# 12. Distance to Forest Plantations
fp_scaled <- seq(min(final.data$dist_fp),
                  max(final.data$dist_fp),
                  by = 0.01)

fp_new <- data.frame(dist_fp = rep(fp_scaled),
                      hfp = mean(final.data$hfp),
                      mesh = mean(final.data$mesh),
                      tri = mean(final.data$tri),
                      gdp = mean(final.data$gdp),
                      dist_forest = mean(final.data$dist_forest),
                      dist_pa = mean(final.data$dist_pa),
                      ed = mean(final.data$ed),
                      dist_cropland = mean(final.data$dist_cropland),
                      popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_fp <- predict(f2.m, newdata = fp_new, re.form = NA, type = "response")
att_fp <- as.data.frame(att_fp)

range(df$dist_fp, na.rm = TRUE)
fp_raw <- seq(from = 0, to = 21824, length.out = 1050) # raw values

fp_new <- cbind(fp_raw, fp_new)

fp_graph <- 
  ggplot() +
  geom_line(data = fp_new, 
            aes(x = fp_raw, 
                y = att_fp$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,21824)) +
  labs(x = "Distance to Forest Plantations (m)",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(fp_graph)

# 12. Distance to Protected Areas
pa_scaled <- seq(min(final.data$dist_pa),
                  max(final.data$dist_pa),
                  by = 0.01)

pa_new <- data.frame(dist_pa = rep(pa_scaled),
                      hfp = mean(final.data$hfp),
                      mesh = mean(final.data$mesh),
                      tri = mean(final.data$tri),
                      gdp = mean(final.data$gdp),
                      dist_forest = mean(final.data$dist_forest),
                      dist_cropland = mean(final.data$dist_cropland),
                      ed = mean(final.data$ed),
                      dist_fp = mean(final.data$dist_fp),
                      popdens = mean(final.data$popdens),
                     wetness = mean(final.data$wetness))

att_pa <- predict(f2.m, newdata = pa_new, re.form = NA, type = "response")
att_pa <- as.data.frame(att_pa)

range(df$dist_pa, na.rm = TRUE)
pa_raw <- seq(from = 0, to = 149313, length.out = 578) # raw values

pa_new <- cbind(pa_raw, pa_new)

pa_graph <- 
  ggplot() +
  geom_line(data = pa_new, 
            aes(x = pa_raw, 
                y = att_pa$V1),
            linetype = 'solid',
            size = 1) +
  coord_cartesian(ylim = c(0,1),
                  xlim = c(0,149313)) +
  labs(x = "Distance to Protected Areas (m)",
       y = "Probability of Human Casualties") +
  theme(legend.position = "NONE",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

plot(pa_graph)

# Save prediction plots
pred_plot <- 
  gridExtra::grid.arrange(ed_graph, 
               mesh_graph, 
               tri_graph,
               w_graph,
               hfp_graph,
               popdens_graph,
               gdp_graph,
               forest_graph,
               crp_graph,
               fp_graph,
               pa_graph,
               nrow = 3)

ggsave('E:/Spatial Analysis/_Supplementary Codebase/_Plots/250720-pred_plot.jpeg',
       pred_plot,
       jpeg,
       width = 30,
       height = 15,
       units = 'in',
       dpi = 300)
