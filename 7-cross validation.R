library(sperrorest)

predict_spaMM <- function(object, newdata, ...) {
  if (!requireNamespace("spaMM", quietly = TRUE)) {
    stop("Package 'spaMM' is required for this function.")
  }
  preds <- spaMM::predict(
    object,
    newdata = newdata,
    re.form = NA,
    type = "response"
  )
  return(as.vector(preds))
}

model_formula <- att ~ ed + mesh + tri + wetness + hfp + popdens + gdp + dist_forest + dist_cropland + dist_fp + dist_pa + I(wetness^2) + I(dist_forest^2) + I(dist_cropland^2) + I(dist_fp^2) + I(dist_pa^2) + Matern(1|y+x)

predict_spaMM <- function(object, newdata, ...) {
  preds <- predict(
    object,
    newdata = newdata,
    re.form = NA, 
    type = "response"
  )
  return(as.vector(preds))
}

nthr <- parallel::detectCores(logical = FALSE) - 1L

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