#
# 02_validation.R
# Leave-region-out (LRO) validation: SPDE-INLA vs LCZ4r (RMSE with bootstrap
# CI over stations) and CRPS of SPDE-INLA (Jan, Jun, Nov 2025; daily median).
#
#   Rscript 02_validation.R tables    # tables from existing LRO predictions (no INLA needed)
#   Rscript 02_validation.R refit     # refit the 15 folds (region x month) and compute CRPS
#
# Design: regions = first digit of the Planning Area code (column `region`).
# For each region, all its official stations are held out; every other
# station (all networks) is used for training. Held-out airports are
# predicted with the `airport` coefficient of the same fit. If GroupRho
# collapses (median > 0.95), the fold is refitted with GroupRho fixed at the
# full-data estimate of the same window (results/fit_<window>.rds from
# script 01; 0.55 if absent). CRPS uses 200 posterior samples of the latent
# field (no station effect, observation noise or airport term).
# CIs: bootstrap over stations (seed 42, 2000 resamples).
# Tables mode reads VALIDATION$lro_spde / lro_lcz from the data file (README).
# 
suppressPackageStartupMessages(library(dplyr))
mode <- commandArgs(TRUE)[1]; if (is.na(mode)) mode <- "tables"
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))

# ---- refit 
if (mode == "refit") {
  here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  source(file.path(here, "01_fit_model.R")); library(scoringRules)
  load(DATA_FILE)

  lro_fold <- function(nm, reg) {
    f <- file.path(OUT_DIR, "lro", sprintf("%s_region%s.rds", nm, reg)); dir.create(dirname(f), showWarnings = FALSE)
    if (file.exists(f)) return(readRDS(f))
    w <- WINDOWS[[nm]]; st <- w$stations
    train <- st[st$region != reg, ]
    test  <- st[st$source == "official" & st$region == reg, ]
    m <- build_model(w, train); fit <- fit_model(m, config = TRUE)

    rho <- fit$summary.hyperpar[grep("GroupRho", rownames(fit$summary.hyperpar))[1], "0.5quant"]
    if (!is.na(rho) && rho > RHO_COLLAPSE) {
      ref <- file.path(OUT_DIR, sprintf("fit_%s.rds", nm)); rho_fix <- RHO_FIXED
      if (file.exists(ref)) { h <- readRDS(ref)$hyper; r <- h[grep("GroupRho", rownames(h))[1], "0.5quant"]
                              if (!is.na(r) && r <= RHO_COLLAPSE) rho_fix <- r }
      fit <- fit_model(m, rho_fix = rho_fix, config = TRUE)
    } else rho_fix <- NA_real_

    # point prediction (posterior means): intercept + covariates + day effect + spatial field
    fx <- fit$summary.fixed[, "mean"]; names(fx) <- rownames(fit$summary.fixed)
    beta_ap <- if ("airport" %in% names(fx)) fx[["airport"]] else 0
    A <- inla.spde.make.A(w$mesh, loc = cbind(test$x, test$y))
    field <- matrix(fit$summary.random$s$mean, nrow = w$mesh$n, ncol = w$n_days)
    k <- cbind(seq_len(nrow(test)), test$day_idx)
    pred0 <- fx[["b0"]] + fx[["elev_std"]] * test$elev_std + fx[["urb_std"]] * test$urb_std +
      fx[["svf_std"]] * test$svf_std + fit$summary.random$day_idx$mean[test$day_idx] + as.matrix(A %*% field)[k]

    # CRPS from posterior samples of the same latent quantity
    smp <- inla.posterior.sample(200, fit, seed = 42); rn <- rownames(smp[[1]]$latent)
    i_f <- match(paste0(c("b0", "elev_std", "urb_std", "svf_std"), ":1"), rn)
    i_d <- grep("^day_idx:", rn); i_s <- grep("^s:", rn)
    X <- cbind(1, test$elev_std, test$urb_std, test$svf_std)
    draws <- sapply(smp, function(s) { l <- s$latent[, 1]
      as.vector(X %*% l[i_f]) + l[i_d][test$day_idx] + as.matrix(A %*% matrix(l[i_s], w$mesh$n, w$n_days))[k] })

    res <- data.frame(window = nm, month = sub("^[a-z]+_", "", nm), region = reg,
                      station_id = test$station_id, date = as.character(test$date), obs = test$temp,
                      pred = pred0 + beta_ap * (test$station_id %in% AIRPORT_STATIONS),
                      rho_fixed = rho_fix, crps = crps_sample(y = test$temp, dat = draws))
    saveRDS(res, f); res
  }

  windows <- commandArgs(TRUE)[-1]; if (!length(windows)) windows <- c("median_2025-01", "median_2025-06", "median_2025-11")
  res <- bind_rows(lapply(windows, function(nm) bind_rows(lapply(as.character(1:5), function(r) lro_fold(nm, r)))))
  write.csv(res, file.path(OUT_DIR, "lro_spde.csv"), row.names = FALSE)
  quit(save = "no")
}

# ---- tables 
DATA_FILE <- Sys.getenv("BRUCE_DATA", "data/bruce_data.RData"); OUT_DIR <- Sys.getenv("BRUCE_OUT", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
e <- new.env(); if (file.exists(DATA_FILE)) load(DATA_FILE, envir = e)
if (exists("VALIDATION", e)) { spde <- e$VALIDATION$lro_spde; lcz <- e$VALIDATION$lro_lcz
} else { spde <- read.csv(file.path(OUT_DIR, "lro_spde.csv")); lcz <- NULL }
airports <- if (exists("AIRPORT_STATIONS", e)) e$AIRPORT_STATIONS else character()

key <- function(d) mutate(d, station_id = as.character(station_id), date = as.character(date), region = as.character(region))
d <- key(spde) %>% mutate(err_spde = pred - obs, airport = station_id %in% airports)
if (!is.null(lcz)) d <- left_join(d, key(lcz)[c("date", "region", "station_id", "err_lcz")], by = c("date", "region", "station_id"))

boot <- function(x, n = 2000) {                      # resample stations (with all their days)
  set.seed(42); idx <- split(seq_len(nrow(x)), x$station_id)
  b <- t(replicate(n, { i <- unlist(idx[sample(length(idx), replace = TRUE)], use.names = FALSE)
                        c(rmse(x$err_spde[i]), rmse(x$err_lcz[i])) }))
  q <- function(v) quantile(v, c(.025, .975))
  data.frame(n = nrow(x), stations = length(idx), rmse_spde = rmse(x$err_spde), rmse_lcz4r = rmse(x$err_lcz),
             spde_lo = q(b[, 1])[1], spde_hi = q(b[, 1])[2], lcz4r_lo = q(b[, 2])[1], lcz4r_hi = q(b[, 2])[2],
             diff = mean(b[, 2] - b[, 1]), diff_lo = q(b[, 2] - b[, 1])[1], diff_hi = q(b[, 2] - b[, 1])[2], row.names = NULL)
}
by_month <- function(x) bind_rows(lapply(split(x, x$month), function(s) cbind(month = s$month[1], boot(s))))

out <- list(spde_all = d %>% group_by(month) %>% summarise(n = n(), rmse_spde = rmse(err_spde)))
if ("err_lcz" %in% names(d)) {
  common <- filter(d, !is.na(err_lcz))
  out$common         <- by_month(common)                       # common station-days, with airports
  out$common_no_airp <- by_month(filter(common, !airport))     # common station-days, without airports
  out$by_region <- common %>% group_by(month, region) %>%
    summarise(n = n(), rmse_spde = rmse(err_spde), rmse_lcz4r = rmse(err_lcz), .groups = "drop")
}
if ("crps" %in% names(d)) out$crps <- d %>% group_by(month) %>% summarise(n = n(), crps = mean(crps, na.rm = TRUE))
for (n in names(out)) { cat("\n===", n, "===\n"); print(as.data.frame(out[[n]]), digits = 3, row.names = FALSE)
                        write.csv(out[[n]], file.path(OUT_DIR, sprintf("validation_%s.csv", n)), row.names = FALSE) }
