# =============================================================================
# 5. Cross-validation
# Spatial cross-validation for the selected model. Run after 2-models.R
# (sources final.data and formula) or load final.data from _Models/final_data.rds.
# =============================================================================

BASE_PATH <- "E:/Spatial Analysis"
# If not run in same session as 2-models.R, load: final.data <- readRDS(file.path(BASE_PATH, "_Models", "final_data.rds"))

library(spaMM)
library(sperrorest)

predict_spaMM <- function(object, newdata, ...) {
  as.vector(spaMM::predict(object, newdata = newdata, re.form = NA, type = "response"))
}

model_formula <- att ~ ed + mesh + tri + wetness + hfp + popdens + gdp +
  dist_forest + dist_cropland + dist_fp + dist_pa +
  I(wetness^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2) +
  Matern(1|y+x)

# Cross-validation
sp_cv <- sperrorest(
  data = final.data,
  formula = model_formula,
  coords = c("x", "y"),
  model_fun = spaMM::fitme,
  model_args = list(family = binomial(link='logit'),
                    weights.form = ~weights),
  pred_fun = predict_spaMM,
  smp_fun = partition_kmeans,
  smp_args = list(
    repetition = 1:2,
    nfold = 5),
  progress = TRUE
)

cv <- summary(sp_cv$error_rep)
