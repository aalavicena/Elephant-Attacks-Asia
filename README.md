# Rising risk of elephant-caused human casualties in Tropical Asia

Supplementary code for **Avicena et al.**, *Rising risk of elephant-caused human casualties in Tropical Asia*. Currently under review in *Communications Sustainability*. Preprint version of the paper is available at ResearchSquare: https://www.researchsquare.com/article/rs-8113717/v1

This repository contains the R code used to build datasets, fit spatial models, perform model selection and inference, explore results, and run cross-validation for mapping elephant-induced human casualty risk in tropical Asia (2015–2050).

---

## Requirements

- **R** (≥ 4.0 recommended)
- **R packages:**  
  `terra`, `dplyr`, `tidyr`, `ggplot2`, `spaMM`, `DHARMa`, `pgirmess`, `gridExtra`, `sf`, `exactextractr`, `sperrorest`  
  Optional: `Hmisc` (Spearman correlation in script 2), `gt` (odds-ratio table in script 3).

Install dependencies from CRAN:

```r
install.packages(c(
  "terra", "dplyr", "tidyr", "ggplot2", "spaMM", "DHARMa",
  "pgirmess", "gridExtra", "sf", "exactextractr", "sperrorest"
))
```

---

## Data (not in this repo)

Scripts expect your own data in a **base directory** (e.g. `E:/~` or `/Volumes/~`). At the top of each script, set:

```r
BASE_PATH <- "E:/~"   # or your path
```

All spatial datasets and GLMM models will be made available in Dryad. You need:

- **Script 1:**  
  - Raster stack: `stacks.tif`  
  - Occurrence CSV: `data.csv` (columns including `longitude`, `latitude`, `human_death`, `human_injury`)  
  - Boundary polygon: `boundary.shp`
- **Script 2:**  
  - Dataset produced by script 1: `data.csv`.
- **Scripts 3–4:**  
  - Model RDS and (for script 4) result rasters and boundary/country shapefiles as referenced in the scripts (paths can be adjusted via `BASE_PATH` and the same folder layout). We provide the 11 models in Dryad since it takes very long time to run each model (~4 days/model on a high-end workstation)

Exact paths and folder names are set inside each script; edit the `BASE_PATH` and any path variables at the top of the file to match your setup.

---

## How to run

Run in order:

1. **`1-build_datasets.R`**  
   Builds occurrence and background point datasets and writes CSV files (e.g. `deahts.csv`, `injuries.csv`, `data.csv`) to `BASE_PATH`.

2. **`2-models.R`**  
   Loads the full dataset (e.g. `data.csv`), prepares and scales predictors, fits candidate and full models with **spaMM**, and saves:
   - `_Models/f2.m_wetness.rds`
   - `_Models/f3.no_sac.rds`
   - `_Models/final_model.rds`

3. **`3-model_selection-inference.R`**  
   Loads the saved model and `final_model.rds`, optionally compares AIC across all `.rds` in `_Models`, runs DHARMa and SAC checks, and produces odds-ratio and effect plots. Outputs go to `_Plots/` (e.g. correlogram, odds plot, prediction plot).

4. **`4-results_exploration.R`**  
   Explores spatial results (distribution of points, risk proportions, baseline vs projected). Set `BASE_PATH` and any script-specific paths (e.g. National Boundaries, WorldPop) to your machine.

5. **`5-cross validation.R`**  
   Spatial cross-validation for the selected model. Expects `final.data` in the R session (run `2-models.R` in the same session) or load it from `_Models/final_model.rds` before running.

---

## Project structure (after running)

```
Elephant-Attacks-Asia/
├── 1-build_datasets.R
├── 2-models.R
├── 3-model_selection-inference.R
├── 4-results_exploration.R
├── 5-cross validation.R
├── README.md
└── REMOVAL_PLAN.md
```

Outputs (created on your machine, not in the repo):

- `BASE_PATH/*.csv` (datasets from script 1)
- `BASE_PATH/_Models/` (model RDS and `final_model.rds`)
- `BASE_PATH/_Plots/` (figures from script 3)
- `BASE_PATH/_Supplementary Codebase/_Results/` (if used by script 4)

---

## Citation

If you use this code, please cite the paper:

*Avicena et al., Rising risk of elephant-caused human casualties in Tropical Asia* (repository: [github.com/aalavicena/Elephant-Attacks-Asia](https://github.com/aalavicena/Elephant-Attacks-Asia)).

---

## License

See repository or paper (pending) for license and terms of use.
