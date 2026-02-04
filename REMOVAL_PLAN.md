# Plan: Unused Files and Functions

## Summary

- **Files:** All 5 R scripts and `README.md` are part of the pipeline; no files are recommended for removal.
- **Functions:** No user-defined functions are unused; the only custom functions (`predict_spaMM` in script 5 and optional `panel.*` in script 2) are used or optional.
- **Package loads:** A few libraries are loaded but not used in a given script; these can be removed or kept for consistency (see below).
- **Model objects:** Several fitted models in `2-models.R` are never saved or used downstream; you can drop them if you only need the selected model.

---

## 1. Package loads (optional cleanup)

| Script | Package | Status | Suggestion |
|--------|---------|--------|------------|
| `2-models.R` | `ggplot2` | Not used in script 2 | Remove if you never plan to add plots here; otherwise keep. |
| `2-models.R` | `terra` | Loaded but not used in script 2 | Safe to remove. |
| `3-model_selection-inference.R` | `pgirmess` | Used for correlogram | Keep. |

All other `library()` calls are used.

---

## 2. Fitted models in `2-models.R`

Only **f2.m_wetness** and **f3.no_sac** are saved and used in scripts 3, 4, and 5.

The following are fitted but **not** saved or used elsewhere:

- `frag.m.w`, `hbr.m.w`, `soc.m.w`, `hec.m.w`, `hd.m.w`, `hec_quad.m.w`, `null.w`, `topo.m`, `f1.m.w`

**Options:**

- **Option A (minimal):** Remove the fits you don’t need (e.g. keep only `f2.m_wetness` and `f3.no_sac`) to shorten run time and simplify the script.
- **Option B (keep for comparison):** Keep all fits and add `saveRDS()` for each (or for a subset) into `_Models/`, so script 3’s AIC comparison (over all `.rds` in `_Models`) includes them. No removal.

---

## 3. Script 3: Optional / redundant code

- **f2.m.no_quad:** Only referenced in commented-out optional comparison. No change needed; leave commented if you don’t use it.
- **Duplicate read of f2.m:** The previous version read `f2.m.rds` and then `f2.m_wetness.rds`; the refactor uses only `f2.m_wetness.rds`. No further removal.

---

## 4. Script 4: Paths and variables

- **Hardcoded paths:** Replaced with `BASE_PATH` and derived paths. Remaining absolute paths (e.g. `National Boundaries/`, WorldPop paths) are noted in the README; consider moving them into a single path config.
- **risk.prop$total_area:** Previously undefined; now set to `Extended_Buffer_Area` before use. No removal.

---

## 5. Repository-level

- **Data/code not in repo:** Scripts expect external data (raster stack, CSV, shapefiles, WorldPop, etc.) in folders like `E:/Spatial Analysis`, `/Volumes/MEC/`, etc. The repo only contains R code and README; no data files to remove.
- **Workspace file:** `avicena.github.io.code-workspace` lives in the parent site repo, not in Elephant-Attacks-Asia; no action in this repo.

---

## Recommended actions (short)

1. **Optional:** In `2-models.R`, remove `library(terra)` and `library(ggplot2)` if you confirm they are unused.
2. **Optional:** In `2-models.R`, remove unused candidate model fits (see list above) if you only need the selected model and want a shorter run.
3. **No file removals:** Keep all 5 R scripts and the README.

After any change, re-run the pipeline (1 → 2 → 3, and 4/5 as needed) to confirm nothing breaks.
