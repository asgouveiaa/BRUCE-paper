# ==================================================================
# 03_exposure.R
# Population exposure to heat on the peak day (2023-11-18, daily median):
# primary-care (APS) territories, with SPDE-INLA posterior uncertainty,
# compared with LCZ4r. Produces Figure 4 (panels a-c) and the tables.
#
#   Rscript 03_exposure.R
#
# Input : data/bruce_data.RData with WINDOWS, AIRPORT_STATIONS and EXPOSURE (README.md)
# Output: results/exposure_city.csv, exposure_aps.csv, figure4_exposure.png
# Time  : ~35 min (refit with config = TRUE; cached in results/ck_exposure/)
#
# Method: 1 km population grid distributed to 500 m pixels by area and
# rescaled so each APS territory sums to its registered population;
# 500 posterior samples of the latent field (no observation noise);
# territory "above" if P(population-weighted T > threshold) > 0.9, "below" if < 0.1.
# LCZ4r: point estimate plus an empirical interval from its same-day
# leave-one-out errors (airports excluded).
# ==================================================================
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(sf); library(terra); library(ggplot2); library(patchwork); library(LCZ4r)
})
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "01_fit_model.R"))              # build_model(), fit_model()
load(DATA_FILE)

WINDOW <- "median_2023-11"; DAY <- as.Date("2023-11-18")
THRESHOLDS <- c(30, 32, 34); MAIN <- 32; N_SAMPLES <- 500; CRS_M <- 31983
ck_dir <- file.path(OUT_DIR, "ck_exposure"); dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) { f <- file.path(ck_dir, name); if (file.exists(f)) readRDS(f) else { v <- expr; saveRDS(v, f); v } }

w <- WINDOWS[[WINDOW]]; mesh <- w$mesh; grid <- w$grid
day_g <- unique(w$stations$day_idx[as.Date(w$stations$date) == DAY]); stopifnot(length(day_g) == 1)
npix <- nrow(grid)
pix_pts <- st_as_sf(grid, coords = c("x", "y"), crs = CRS_M, remove = FALSE) %>% mutate(pix = row_number())
pixels  <- st_buffer(pix_pts, 250, endCapStyle = "SQUARE")                  # 500 m cells
A_grid  <- inla.spde.make.A(mesh, loc = as.matrix(grid[, c("x", "y")]))

# ---- 1. Population per pixel and territory x pixel weights -----------------------------
aps <- EXPOSURE$aps %>% st_transform(CRS_M) %>% st_make_valid() %>% mutate(aps_id = as.character(id))
pop_pix <- cached("pop_pixel.rds", {
  pop1km <- EXPOSURE$pop_1km %>% st_transform(CRS_M) %>% filter(pop > 0)
  coalesce(st_interpolate_aw(pop1km["pop"], pixels["pix"], extensive = TRUE, keep_NA = TRUE)$pop, 0)
})
inter <- cached("pixel_aps_intersection.rds",
  st_intersection(pixels["pix"], aps["aps_id"]) %>% mutate(area = as.numeric(st_area(geometry))) %>%
    st_drop_geometry() %>% select(pix, aps_id, area))
inter <- inter %>% mutate(pop = pop_pix[pix] * area / 500^2) %>% group_by(aps_id) %>%
  mutate(total = aps$population[match(aps_id[1], aps$aps_id)],       # grid gives the spatial split,
         pop = if (sum(pop) > 0) pop / sum(pop) * total else area / sum(area) * total) %>%  # registry the total
  ungroup()
pop <- numeric(npix); p <- inter %>% group_by(pix) %>% summarise(pop = sum(pop)); pop[p$pix] <- p$pop

f <- factor(inter$aps_id)
W <- sparseMatrix(i = as.integer(f), j = inter$pix, x = inter$pop, dims = c(nlevels(f), npix), dimnames = list(levels(f), NULL))
W <- W[rowSums(W) > 0, , drop = FALSE]
W_norm <- Diagonal(x = 1 / rowSums(W)) %*% W; dimnames(W_norm) <- dimnames(W)
cat(sprintf("Population used: %.0f | registry: %.0f\n", sum(pop), sum(aps$population)))

# ---- 2. Posterior samples of the latent field on the day --------------------------------
S <- cached(sprintf("samples_%s.rds", DAY), {
  fit <- fit_model(build_model(w), rho_fix = if (WINDOW %in% RHO_FIXED_WINDOWS) RHO_FIXED, config = TRUE)
  idx <- inla.spde.make.index("s", n.spde = mesh$n, n.group = w$n_days)
  pos <- which(idx$s.group == day_g); nodes <- idx$s[pos]
  pos_d <- which(fit$summary.random$day_idx$ID == day_g)
  smp <- inla.posterior.sample(N_SAMPLES, fit, seed = 42); rn <- rownames(smp[[1]]$latent)
  i_s <- match(paste0("s:", pos), rn); i_d <- match(paste0("day_idx:", pos_d), rn)
  i_f <- match(paste0(c("b0", "elev_std", "urb_std", "svf_std"), ":1"), rn)
  X <- cbind(1, grid$elev_std, grid$urb_std, grid$svf_std)
  sapply(smp, function(s) { l <- s$latent[, 1]; field <- numeric(mesh$n); field[nodes] <- l[i_s]
                            as.vector(X %*% l[i_f]) + l[i_d] + as.vector(A_grid %*% field) })   # pixels x samples
})

# ---- 3. LCZ4r on the same day: map + leave-one-out errors --------------------------------
st <- w$stations %>% filter(as.Date(date) == DAY, !is.na(temp)) %>%
  transmute(estacao_id = station_id, longitude, latitude, temp, date = DAY)    # LCZ4r column names
lcz <- terra::unwrap(EXPOSURE$lcz)
dd <- list(year = as.integer(format(DAY, "%Y")), month = as.integer(format(DAY, "%m")), day = as.integer(format(DAY, "%d")))

lcz_pix <- cached(sprintf("lcz_map_%s.rds", DAY), {
  m <- do.call(lcz_interp_map, c(list(lcz, data_frame = st, var = "temp", vg.model = "Ste", station_id = "estacao_id",
                                      sp.res = 500, tp.res = "day"), dd))
  v <- terra::extract(m, terra::project(terra::vect(pixels), terra::crs(m)), fun = mean, na.rm = TRUE)[, 2]
  v[is.nan(v)] <- NA; v                                                          # mean of 100 m cells per pixel
})
ev <- as.data.frame(cached(sprintf("lcz_loocv_%s.rds", DAY),
  do.call(lcz_interp_eval, c(list(lcz, data_frame = st, var = "temp", station_id = "estacao_id", LOOCV = TRUE, sp.res = 100,
                                  tp.res = "day", vg.model = "Ste", extract.method = "simple", LCZinterp = TRUE), dd))))
cat("lcz_interp_eval columns:", paste(names(ev), collapse = ", "), "\n")      # check the three names below
col <- function(pattern) names(ev)[grepl(pattern, names(ev), ignore.case = TRUE) & sapply(ev, is.numeric)][1]
c_obs <- col("obs|observ|measur|real"); c_pred <- col("pred|estim|interp|krig|fit")
c_id  <- names(ev)[grepl("station|estacao|id", names(ev), ignore.case = TRUE)][1]
err <- (ev[[c_pred]] - ev[[c_obs]])[!(ev[[c_id]] %in% AIRPORT_STATIONS)]; err <- err[is.finite(err)]
q_err <- quantile(err, c(0.025, 0.975))
cat(sprintf("LCZ4r LOOCV: n = %d, RMSE = %.2f, bias = %.2f\n", length(err), sqrt(mean(err^2)), mean(err)))
na_lcz <- is.na(lcz_pix)

# ---- 4. Tables ------------------------------------------------------------------------------
fmt <- function(x) format(round(x), big.mark = ",")
city <- bind_rows(lapply(THRESHOLDS, function(L) {
  exposed <- colSums((S > L) * pop)
  exposed_lcz <- sapply(err, function(e) sum(pop[!na_lcz][(lcz_pix[!na_lcz] - e) > L]))   # true = predicted - error
  data.frame(threshold = L,
             exposed_spde = sprintf("%s (%s-%s)", fmt(mean(exposed)), fmt(quantile(exposed, .025)), fmt(quantile(exposed, .975))),
             pct_spde = sprintf("%.1f%%", 100 * mean(exposed) / sum(pop)),
             exposed_lcz = fmt(sum(pop[which(lcz_pix > L)])),
             exposed_lcz_interval = sprintf("%s-%s", fmt(quantile(exposed_lcz, .025)), fmt(quantile(exposed_lcz, .975))),
             pop_without_lcz_estimate = fmt(sum(pop[na_lcz])))
}))
print(city)

T_aps <- as.matrix(W_norm %*% S)                                   # territories x samples
W_raw <- W[rownames(W_norm), , drop = FALSE]
w_ok  <- as.vector(W_raw %*% as.numeric(!na_lcz))
t_lcz <- as.vector(W_raw %*% ifelse(na_lcz, 0, lcz_pix)) / w_ok; t_lcz[w_ok == 0] <- NA
classify <- function(p) case_when(is.na(p) ~ "no estimate", p > 0.9 ~ "above", p < 0.1 ~ "below", TRUE ~ "uncertain")

terr <- tibble(aps = rownames(W_norm), pop = aps$population[match(rownames(W_norm), aps$aps_id)],
               temp_spde = rowMeans(T_aps), lo = apply(T_aps, 1, quantile, .025), hi = apply(T_aps, 1, quantile, .975),
               p_above = rowMeans(T_aps > MAIN), temp_lcz = t_lcz, lcz_coverage = w_ok / rowSums(W_raw),
               lcz_lo = t_lcz - q_err[[2]], lcz_hi = t_lcz - q_err[[1]],
               p_above_lcz = sapply(t_lcz, function(t) if (is.na(t)) NA_real_ else mean((t - err) > MAIN))) %>%
  mutate(class_spde = classify(p_above), class_lcz = classify(p_above_lcz)) %>% arrange(desc(p_above))
cat(sprintf("\nClassification at %d C:\n", MAIN)); print(table(SPDE = terr$class_spde, LCZ4r = terr$class_lcz))
write.csv(city, file.path(OUT_DIR, "exposure_city.csv"), row.names = FALSE)
write.csv(terr, file.path(OUT_DIR, "exposure_aps.csv"), row.names = FALSE)

# ---- 5. Figure 4 (a-c) ------------------------------------------------------------------------
g <- tibble(x = grid$x, y = grid$y, mean = rowMeans(S), p = rowMeans(S > MAIN), na_lcz = na_lcz)
theme_map <- theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
                                     plot.title = element_text(size = 9), legend.key.size = unit(0.35, "cm"))
pa <- ggplot(g, aes(x, y, fill = mean)) + geom_tile() + coord_equal() +
  scale_fill_gradient(low = "#fee090", high = "#a50026", name = "\u00b0C") +
  labs(title = sprintf("(a) SPDE-INLA posterior mean, daily median, %s", DAY), x = NULL, y = NULL) + theme_map
pb <- ggplot(g, aes(x, y, fill = p)) + geom_tile() + coord_equal() +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0, 1), name = "P") +
  geom_tile(data = filter(g, na_lcz), fill = NA, color = "black", linewidth = 0.2) +
  labs(title = sprintf("(b) P(T > %d\u00b0C); outlined: no LCZ4r estimate", MAIN), x = NULL, y = NULL) + theme_map
pc <- terr %>% arrange(temp_spde) %>% mutate(i = row_number()) %>% ggplot(aes(i)) +
  geom_hline(yintercept = MAIN, linetype = 2, color = "grey40") +
  geom_linerange(aes(ymin = lcz_lo, ymax = lcz_hi), color = "grey70", linewidth = 0.4, position = position_nudge(x = 0.35)) +
  geom_point(aes(y = temp_lcz), shape = 4, size = 0.8, color = "grey30", position = position_nudge(x = 0.35)) +
  geom_linerange(aes(ymin = lo, ymax = hi, color = p_above), linewidth = 0.4) +
  geom_point(aes(y = temp_spde, color = p_above), size = 0.5) +
  scale_color_distiller(palette = "Reds", direction = 1, limits = c(0, 1), name = "P(>thr)\nSPDE") +
  labs(title = "(c) Primary care territories: SPDE-INLA 95% CrI (coloured) vs LCZ4r estimate (x) with empirical 95% interval\nfrom same-day leave-one-out errors (grey; spatially constant width)",
       x = "Territories, ordered by SPDE-INLA estimate", y = "Population-weighted temperature (\u00b0C)") +
  theme_minimal() + theme(plot.title = element_text(size = 12), axis.text.x = element_blank())
ggsave(file.path(OUT_DIR, "figure4_exposure.png"), (pa | pb) / pc, width = 12, height = 8, dpi = 300)
