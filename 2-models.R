# =============================================================================
# 2. Fit models
# Prepares data, fits candidate GLMMs with spaMM, and saves model objects.
# Run after 1-build_datasets.R (expects 250617_all_datasets.csv or 250717).
# =============================================================================

BASE_PATH <- "E:/Spatial Analysis"
setwd(BASE_PATH)

library(terra)
library(dplyr)
library(ggplot2)
library(DHARMa)
library(spaMM)

# Optional: Spearman correlation (uncomment if needed)
# library(Hmisc)

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------
df <- read.csv(file.path(BASE_PATH, "250717_all_datasets.csv"))
df$att <- as.factor(df$att)
df$dist_forest <- df$dist_forest * 1000  # convert to metres

# Scale predictors and build analysis dataset (no c.cover / c.heights)
predictors <- c(
  "ed", "mesh", "hfp", "dist_pa", "gdp", "popdens", "tri",
  "dist_forest", "dist_cropland", "dist_fp", "wetness", "x", "y"
)
final.data <- df %>%
  select(att, all_of(predictors)) %>%
  mutate(across(
    c(ed, mesh, hfp, dist_pa, gdp, popdens, tri, dist_forest, dist_cropland, dist_fp, wetness),
    ~ c(scale(.))
  ))

# Case weights: presence:absence ratio
n_presence <- sum(final.data$att == 1L)
n_absence  <- sum(final.data$att == 0L)
final.data$weights <- ifelse(final.data$att == 1L, 1, n_presence / n_absence)

final.data <- distinct(final.data, x, y, .keep_all = TRUE)
final.data <- na.omit(final.data)

# -----------------------------------------------------------------------------
# Correlation (optional): uncomment to run
# -----------------------------------------------------------------------------
# panel.hist <- function(x, ...) {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(usr[1:2], 0, 1.5))
#   h <- hist(x, plot = FALSE)
#   breaks <- h$breaks; nB <- length(breaks)
#   y <- h$counts / max(h$counts)
#   rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
# }
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y, use = "pairwise.complete.obs"))
#   txt <- format(c(r, 0.12345678), digits = digits)[1]
#   if (missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
#   text(0.5, 0.5, paste0(prefix, txt), cex = cex.cor)
# }
# pairs(final.data, upper.panel = panel.cor, diag.panel = panel.hist)
# final.data %>% as.matrix() %>% Hmisc::rcorr(., type = "spearman")

# -----------------------------------------------------------------------------
# Model fitting
# -----------------------------------------------------------------------------
nthr <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
hl_opt <- list(NbThreads = nthr)
wform <- ~ weights

# Candidate models
frag.m.w     <- fitme(att ~ ed + mesh + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
hbr.m.w      <- fitme(att ~ wetness + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
soc.m.w      <- fitme(att ~ popdens + gdp + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
hec.m.w      <- fitme(att ~ dist_pa + dist_forest + dist_cropland + dist_fp + Matern(1|y+x),
                     data = final.data, family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
hd.m.w       <- fitme(att ~ hfp + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
hec_quad.m.w <- fitme(att ~ dist_pa + dist_forest + dist_cropland + dist_fp +
                       I(dist_pa^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + Matern(1|y+x),
                     data = final.data, family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
null.w       <- fitme(att ~ 1 + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)
topo.m       <- fitme(att ~ tri + Matern(1|y+x), data = final.data,
                     family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)

# Full models
f1.m.w <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp +
                  dist_forest + dist_cropland + dist_fp + dist_pa + Matern(1|y+x),
                data = final.data, family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)

f2.m_wetness <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp +
                        dist_forest + dist_cropland + dist_fp + dist_pa +
                        I(wetness^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2) +
                        Matern(1|y+x),
                      data = final.data, family = binomial(link = "logit"), weights.form = wform, control.HLfit = hl_opt)

f3.no_sac <- fitme(att ~ ed + mesh + tri + wetness + hfp + popdens + gdp +
                     dist_forest + dist_cropland + dist_fp + dist_pa +
                     I(wetness^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2),
                   data = final.data, family = binomial(link = "logit"), weights.form = wform)

# Save key objects for downstream scripts (create dirs if needed)
dir_models <- file.path(BASE_PATH, "_Models", "sets 3")
dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
saveRDS(f2.m_wetness, file.path(dir_models, "f2.m_wetness.rds"))
saveRDS(f3.no_sac,   file.path(dir_models, "f3.no_sac.rds"))
# Save final.data for 3-model_selection-inference.R and 5-cross validation.R
saveRDS(final.data, file.path(BASE_PATH, "_Models", "final_data.rds"))
