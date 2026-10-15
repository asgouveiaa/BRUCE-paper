# 01_fit_model.R
# BRUCE: Bayesian spatiotemporal SPDE-INLA model for daily air temperature
# (daily median or maximum), fitted on monthly windows.
#
# Defines build_model() and fit_model() (also used by 02 and 03) and, when
# run directly, fits the requested windows and saves coefficients and
# hyperparameters.
#
#   Rscript 01_fit_model.R                    # all windows in WINDOWS
#   Rscript 01_fit_model.R median_2025-01     # selected window(s)
#
# Input : data/bruce_data.RData (see README.md); override with env BRUCE_DATA
# Output: results/fit_<window>.rds, results/fixed_effects.csv, results/hyperparameters.csv
# Time  : ~34 min per window; windows already fitted are skipped

suppressPackageStartupMessages({ library(INLA); library(dplyr) })

DATA_FILE <- Sys.getenv("BRUCE_DATA", "data/bruce_data.RData")
OUT_DIR   <- Sys.getenv("BRUCE_OUT", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Model settings (values used in the paper) ----------------------------------
SIGMA_ERA5  <- 0.7             # fixed noise SD, ERA5-Land (deg C)
SIGMA_OFIC  <- 0.3             # fixed noise SD, official stations (deg C)
RHO_FIXED   <- 0.55            # GroupRho fixed value in windows where it collapses to 1
RHO_COLLAPSE <- 0.95           # posterior median of GroupRho above this = collapsed
# windows where GroupRho collapsed and was fixed (10 of 18)
RHO_FIXED_WINDOWS <- c(
  paste0("max_",    c("2023-01", "2023-06", "2023-11", "2024-01", "2024-06", "2024-11", "2025-01")),
  paste0("median_", c("2024-01", "2024-06", "2024-11")))

fixed_priors <- list(
  mean = list(elev_std = 0, urb_std = 0, svf_std = 0, airport = 0, aura = 0, wunder = 0, default = 0),
  prec = list(elev_std = 1/4, urb_std = 1/4, svf_std = 1/4, airport = 1/9, aura = 1/9, wunder = 1/9, default = 0.001),
  mean.intercept = 25, prec.intercept = 1/100)

fixed_noise <- function(sigma) list(prec = list(initial = log(1 / sigma^2), fixed = TRUE))
pc_noise    <- list(prec = list(prior = "pc.prec", param = c(3, 0.01)))   # aura and wunder

# ---- Observation model: one likelihood per data source ---------------------------
# w: one element of WINDOWS; train: stations used for fitting (default: all).
# Uses AIRPORT_STATIONS from the data file.
build_model <- function(w, train = w$stations) {
  spde <- inla.spde2.pcmatern(w$mesh, alpha = 2, constr = TRUE,
                              prior.range = c(15000, 0.5), prior.sigma = c(8, 0.01))
  field <- inla.spde.make.index("s", n.spde = spde$n.spde, n.group = w$n_days)
  all_ids <- unique(w$stations$station_id)       # station effect indexed over all stations

  src <- split(train, factor(train$source, c("official", "aura", "wunder")))
  has <- sapply(src, nrow) > 0                   # official always present
  n_lik <- 1 + sum(has)                          # ERA5 + available station sources
  lik_col <- c(official = 2, aura = if (has["aura"]) 3, wunder = if (has["wunder"]) 2 + has["aura"] + 1)

  make_stack <- function(tag, d, col, is_station, aura = 0, wunder = 0, airport = TRUE) {
    y <- matrix(NA_real_, nrow(d), n_lik); y[, col] <- d$temp
    eff <- data.frame(b0 = 1, aura = aura, wunder = wunder,
                      airport = if (airport) as.integer(d$station_id %in% AIRPORT_STATIONS) else 0,
                      elev_std = d$elev_std, urb_std = d$urb_std, svf_std = d$svf_std, day_idx = d$day_idx,
                      station_idx = if (is_station) match(d$station_id, all_ids) else NA_integer_)
    A <- inla.spde.make.A(w$mesh, loc = cbind(d$x, d$y), group = d$day_idx, n.group = w$n_days)
    inla.stack(tag = tag, data = list(y = y), A = list(1, A), effects = list(eff, s = field))
  }

  stacks <- list(make_stack("era5", w$era5, 1, FALSE, airport = FALSE),
                 make_stack("official", src$official, 2, TRUE))
  noise  <- list(fixed_noise(SIGMA_ERA5), fixed_noise(SIGMA_OFIC))
  if (has["aura"]) {
    stacks <- c(stacks, list(make_stack("aura", src$aura, lik_col["aura"], TRUE, aura = 1, airport = FALSE)))
    noise <- c(noise, list(pc_noise))
  }
  if (has["wunder"]) {
    stacks <- c(stacks, list(make_stack("wunder", src$wunder, lik_col["wunder"], TRUE, wunder = 1, airport = FALSE)))
    noise <- c(noise, list(pc_noise))
  }
  list(stack = do.call(inla.stack, stacks), noise = noise, n_lik = n_lik, spde = spde,
       has_aura = has[["aura"]], has_wunder = has[["wunder"]],
       has_airport = any(src$official$station_id %in% AIRPORT_STATIONS))
}

# ---- Fit
# rho_fix = NULL: GroupRho estimated (PC prior); a number: GroupRho fixed at that value.
fit_model <- function(m, rho_fix = NULL, config = FALSE, cpo = FALSE) {
  spde <- m$spde                                  # the formula looks up `spde` in this environment
  rho_group <- if (is.null(rho_fix)) "hyper = list(rho = list(prior = 'pc.cor1', param = c(0.5, 0.75)))"
               else sprintf("hyper = list(rho = list(initial = %.6f, fixed = TRUE))", log((1 + rho_fix) / (1 - rho_fix)))
  form <- as.formula(paste(
    "y ~ 0 + b0", if (m$has_aura) "+ aura", if (m$has_wunder) "+ wunder", if (m$has_airport) "+ airport",
    "+ elev_std + urb_std + svf_std",
    "+ f(station_idx, model = 'iid', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))",
    "+ f(day_idx, model = 'ar1', hyper = list(rho = list(prior = 'pc.cor1', param = c(0.5, 0.75)),",
    "                                         prec = list(prior = 'pc.prec', param = c(3, 0.01))))",
    "+ f(s, model = spde, group = s.group, control.group = list(model = 'ar1',", rho_group, "))"))
  environment(form) <- environment()
  inla(form, family = rep("gaussian", m$n_lik), data = inla.stack.data(m$stack),
       control.predictor = list(compute = TRUE, A = inla.stack.A(m$stack)),
       control.family = lapply(m$noise, function(h) list(hyper = h)), control.fixed = fixed_priors,
       control.compute = list(cpo = cpo, config = config), control.inla = list(int.strategy = "ccd"))
}

# 
# Run directly: fit windows and save summaries
# 
if (sys.nframe() == 0) {
  load(DATA_FILE)
  todo <- commandArgs(TRUE); if (!length(todo)) todo <- names(WINDOWS)
  for (nm in todo) {
    f <- file.path(OUT_DIR, sprintf("fit_%s.rds", nm))
    if (file.exists(f)) { cat("skipping (exists):", f, "\n"); next }
    cat("\n===", nm, "===\n"); t0 <- proc.time()
    rho_fix <- if (nm %in% RHO_FIXED_WINDOWS) RHO_FIXED else NULL
    fit <- fit_model(build_model(WINDOWS[[nm]]), rho_fix = rho_fix, cpo = TRUE)
    saveRDS(list(window = nm, rho_fixed = rho_fix, fixed = fit$summary.fixed,
                 hyper = fit$summary.hyperpar, mlik = fit$mlik[1, 1]), f)
    cat(sprintf("done in %.1f min\n", (proc.time() - t0)["elapsed"] / 60)); print(fit$summary.fixed[, 1:5])
  }
  res <- lapply(list.files(OUT_DIR, "^fit_.*[0-9]\\.rds$", full.names = TRUE), readRDS)
  tab <- function(el) bind_rows(lapply(res, function(r) data.frame(window = r$window, term = rownames(r[[el]]), r[[el]], row.names = NULL)))
  write.csv(tab("fixed"), file.path(OUT_DIR, "fixed_effects.csv"), row.names = FALSE)
  write.csv(tab("hyper"), file.path(OUT_DIR, "hyperparameters.csv"), row.names = FALSE)
}
