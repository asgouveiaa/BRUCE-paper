# BRUCE: supporting code

Code to reproduce the modelling, validation and exposure analysis of the paper (Bayesian spatiotemporal
SPDE-INLA interpolation of air temperature in Rio de Janeiro, compared with LCZ4r). The `bruce` package
(data retrieval and preprocessing) is not distributed. **The data are not included**: one of the sensor
networks belongs to a private project. They can be requested from the authors (contact in the paper).

## Scripts

| Script | Purpose | Reproduces |
|---|---|---|
| `01_fit_model.R` | Model functions; fits the monthly windows (median and maximum) | Fixed effects and hyperparameters (Tables 2-3) |
| `02_validation.R` | Leave-region-out validation, SPDE-INLA vs LCZ4r (RMSE, bootstrap CI), CRPS | Validation tables, S1 |
| `03_exposure.R` | Exposed population and APS territory classification on 2023-11-18 | Figure 4 |

```
Rscript 01_fit_model.R [window ...]       # e.g. median_2025-01  (~30 min per window)
Rscript 02_validation.R tables            # tables from existing predictions (no INLA)
Rscript 02_validation.R refit [window ...]  # refit the 15 folds (region x month) and CRPS
Rscript 03_exposure.R
```
Env variables: `BRUCE_DATA` (default `data/bruce_data.RData`), `BRUCE_OUT` (default `results`).
Packages: INLA, dplyr, Matrix, sf, terra, ggplot2, patchwork, scoringRules, LCZ4r.

## Data file (`.RData`)

All objects are plain data frames (no geometry column) unless stated. Each element of `WINDOWS` is already preprocessed.

- `AIRPORT_STATIONS`: `station_id` of airport stations.
- `WINDOWS`: named list `"<median|max>_<YYYY-MM>"`; each element has
  - `stations`: station-days with `station_id`, `source` (`official`/`aura`/`wunder`), `date`, `x`, `y` (UTM 23S, m), `longitude`, `latitude`, `temp` (target variable), `elev_std`, `urb_std`, `svf_std`, `day_idx`, `region` (first digit of the Planning Area code, as character);
  - `era5`: ERA5-Land pixel-days with `x`, `y`, `temp`, `elev_std`, `urb_std`, `svf_std`, `day_idx`;
  - `grid`: 500 m prediction grid with `x`, `y`, `elev_std`, `urb_std`, `svf_std`;
  - `mesh` (`inla.mesh.2d`) and `n_days`.
- `VALIDATION` (only for `02_validation.R tables`): `lro_spde` (`station_id`, `date`, `region`, `month`, `obs`, `pred`, `crps`) and `lro_lcz` (`station_id`, `date`, `region`, `err_lcz` = prediction - observation).
- `EXPOSURE` (for `03_exposure.R`): `aps` (sf, APS territories with `id`, `population`), `pop_1km` (sf, 1 km population grid, column `pop`), `lcz` (`terra::wrap()` of the LCZ map in LCZ4r format).

## Notes

- Fixed noise SDs (ERA5 0.7, official stations 0.3), priors, and the 10 of 18 windows in which GroupRho collapsed
  and was fixed at 0.55 (not estimated) are listed at the top of `01_fit_model.R`.
- The LCZ4r leave-region-out predictions are supplied in `VALIDATION$lro_lcz`; they are not recomputed here.
- Everything that uses sampling or resampling uses seed 42.
