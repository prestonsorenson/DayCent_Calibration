# =============================================================================
# DayCent RM-level Bayesian calibration — ALL CROPS
#   smoke -> split -> Morris -> Sobol -> automatic selection -> N-response bound
#   for PRDX(1) -> LHS design -> GP emulators -> DEoptim (mode) -> DREAM(ZS)
#   -> posterior evaluation against CSU defaults on spatially/temporally
#   withheld data -> N-response curve
#
#   * one loop over crops; each crop writes to <crop>/<g$out_name>/
#   * calibrated parameters chosen by the Sobol rule (lower 95% CI of max
#     total-order index >= sobol_ST_thresh, at most max_calib), minus
#     never_calibrate
#   * ranges: absolute where the DayCent manual gives limits, otherwise
#     relative to each crop's CSU template value
#   * FINAL CONFIGURATION: PPDF(1) and PPDF(2) always calibrated within crop-specific
#     physiological ranges (g$temp_priors); variance-matched likelihood with
#     component emulation (g$lik = "vmc"); temporal effective sample size = number
#     of training years.
#   * (earlier test run) temperature response: PPDF(1) and PPDF_GAP (= PPDF(2) - PPDF(1),
#     8-25 C) screened separately; if either passes Sobol selection, both are
#     calibrated (paired rule). The gap parameterisation prevents narrow curves
#     (PPDF(2) close to PPDF(1)) that give zero production and crash DayCent
#     ("demand = 0.0 in nutrlm").
#   * HIMAX screened but never calibrated (fixed at CSU value): grain yield alone
#     cannot separate harvest index from production
#   * PRDX(1) upper bound: PRDX(1) at which yield at the recommended N rate
#     falls to nresp_min of yield at 2x N, INTERPOLATED on a fine grid
#     (skipped for pulses)
#   * priors: truncated normal centred on the CSU template (or uniform);
#     priors are NOT informed by the calibration data
#   * DEoptim on the emulated log-posterior locates the posterior mode as an
#     independent check of the DREAM(ZS) median (distance reported in
#     posterior-SD units); it does not inform priors or chain starts
#   * 100 DayCent runs at posterior draws check the emulator in the posterior region.
#     If posterior-region RMSE > refine$rmse_tol x hold-out RMSE or 95% coverage <
#     refine$cover_min for any emulated term, the emulator is refined automatically:
#     DayCent runs at posterior draws + a space-filling set in a box around the
#     posterior, global GPs refit on all runs, local GPs fit inside the box, DREAM(ZS)
#     rerun within the box, rechecked; up to refine$max_rounds. The point estimate is
#     the refined posterior median unless its DayCent run is > pool_delta_ll below the
#     best DayCent run, in which case the top-set medoid of all DayCent runs is used.
#   * likelihood: Gaussian on variance-matched spatial and temporal SSE (simulation
#     rescaled to observed SD before scoring) plus a lognormal penalty on the
#     variability ratio alpha (sd g$alpha_sd). Raw SSE ("sse") rewards damping
#     simulated variability toward r * obs SD, which pushes stress parameters
#     (e.g. WSCOEF(1)) to their bounds.
#     "vm"  emulates log SSE directly.
#     "vmc" emulates the components (spatial bias, atanh r spatial, atanh r temporal,
#           log alpha spatial/temporal) and computes the SSE exactly:
#             SSE_s = n_s [b^2 + 2 v_s (1 - r_s)],  SSE_t = 2 n_t v_t (1 - r_t)
#           (v_s, v_t: population variances of observed RM means and anomalies).
#           Near the optimum SSE_s is dominated by the squared bias, a sharp valley
#           a GP fits poorly on the log scale; the bias itself is near-linear.
#
# Run:  nohup Rscript all_crops_gp_dreamzs.R > log_all_crops.txt 2>&1 &
#       (outside conda: the conda PROJ database conflicts with sf)
# Every phase checkpoints; rerunning skips cached DayCent runs.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(parallel); library(sensitivity); library(lhs)
  library(DiceKriging); library(BayesianTools); library(coda); library(DEoptim)
})

# =============================================================================
# 1. GLOBAL SETTINGS
# =============================================================================

models_root <- "/home/preston/Projects/Daycent/Crop_Optimization/RM/Models"

crop_table <- tribble(
  ~folder,        ~crop_id, ~pulse,
  "Barley",       "BAR3",   FALSE,
  "Wheat",        "SW0",    FALSE,
  "Durum",        "DURUM",  FALSE,
  "Oats",         "OAT1",   FALSE,
  "Canola",       "CAN",    FALSE,
  "Peas",         "PEA",    TRUE,
  "Green_Lentil", "GLENT",  TRUE,
  "Red_Lentil",   "RLENT",  TRUE
)
crops_to_run <- crop_table$folder        # all crops; or e.g. c("Wheat", "Canola")

phases <- c("smoke", "split", "morris", "sobol", "select", "nbound",
            "design", "gp", "mcmc", "check", "refine", "post", "mgmt")

g <- list(
  sites_subdir = "qrts", exe = "DDcentEVI_rev491", schedule = "Crop_Opt.sch",
  out_prefix = "harvest", obs_file = "observed.csv",
  rm_coords_path = "/mnt/Elements/Saskatchewan/RM_Shapefile/RM_Centroids.shp",
  out_name = "gp_dreamzs_final_vmc",      # component-emulated variance-matched likelihood
  seed_from = "gp_dreamzs_final",         # reuse cached runs (incl. refinement runs) from the previous final run
  # likelihood
  lik      = "vmc",     # "sse" = raw SSE; "vm" = variance-matched log-SSE emulated directly;
  # "vmc" = variance-matched, components emulated and SSE computed exactly
  alpha_sd = 0.10,      # sd of log(alpha) penalty; alpha = sim SD / obs SD (spatial and temporal)
  n_cores = 100, seed = 20261013,
  # split
  n_blocks = 10, test_block_frac = 0.20, test_frac_tol = 0.05,
  test_year_start = 2016, min_years_rm = 4,
  # screening
  morris_r = 20, morris_levels = 6, morris_thresh = 0.10, max_fail_frac = 0.05,
  sobol_N = 500, sobol_nboot = 200,
  # selection
  sobol_ST_thresh = 0.05,        # applied to the LOWER 95% CI of max ST
  max_calib = 8,
  paired = list(),
  never_calibrate = c("FRTC(1)", "FRTC(3)", "HIMAX", "PPDF(3)", "PPDF(4)"),   # not identifiable from grain yield alone
  # Temperature response: PPDF(1) (optimum) and PPDF(2) (maximum) are calibrated for every crop,
  # bounded by physiologically justified ranges (deg C; cite sources in the paper).
  # The two ranges must be separated by at least min_temp_gap (PPDF(2) - PPDF(1)) so DayCent
  # never gets a degenerate curve ("demand = 0.0 in nutrlm"). NA = fill in before that crop runs.
  always_calibrate = c("PPDF(1)", "PPDF(2)"),
  min_temp_gap = 8,
  temp_priors = list(
    Barley       = list(`PPDF(1)` = c(15, 22), `PPDF(2)` = c(31, 40)),
    Wheat        = list(`PPDF(1)` = c(17, 22), `PPDF(2)` = c(32, 40)),
    Durum        = list(`PPDF(1)` = c(17, 22), `PPDF(2)` = c(32, 40)),
    Oats         = list(`PPDF(1)` = c(15, 22), `PPDF(2)` = c(31, 37)),
    Canola       = list(`PPDF(1)` = c(17, 22), `PPDF(2)` = c(30, 36)),
    Peas         = list(`PPDF(1)` = c(14, 20), `PPDF(2)` = c(28, 35)),
    Green_Lentil = list(`PPDF(1)` = c(16, 22), `PPDF(2)` = c(30, 38)),
    Red_Lentil   = list(`PPDF(1)` = c(16, 22), `PPDF(2)` = c(30, 38))),
  # always calibrated for a crop, whatever Sobol says (paired partners are added automatically);
  # they take priority within max_calib, lowest-ranked Sobol picks are dropped to make room
  force_calibrate = list(),
  # N-response bound for PRDX(1): fine grid near 1x, crossing interpolated
  nresp_min = 0.90, nresp_rms = 8,
  nresp_grid = c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2, 2.5),
  # priors
  prior_type = "tnorm",          # "tnorm" (centred on CSU template) or "uniform"
  prior_sd_frac = 0.25,          # sd = frac * (upper - lower)
  # DEoptim on the emulator (posterior-mode check)
  de_np_per_dim = 10, de_itermax = 300,
  # design / emulator / MCMC
  design_per_dim = 80, design_min = 300, n_val = 50, n_check = 100,
  neff_spatial_frac = 0.5, mcmc_iter = 60000, mcmc_chains = 3, burn_frac = 0.5,
  # emulator refinement near the posterior (runs only if the posterior-region check fails)
  refine = list(max_rounds = 3, n_post_pts = 100, n_box_pts = 50, box_q = c(0.005, 0.995),
                box_widen = 0.5, min_local_per_dim = 15, rmse_tol = 1.5, cover_min = 0.80,
                rhat_max = 1.2, pool_delta_ll = 2, pool_min_top = 20),
  # evaluation
  n_post = 60, n_boot = 1000, n_worst_years = 2,
  mgmt_mult = c(0.7, 1, 1.3, 2)
)

# =============================================================================
# 2. CANDIDATE PARAMETERS
#    type "abs": lower/upper used as given (DayCent manual limits or stability)
#    type "rel": lower/upper are multipliers of the crop's template value
# =============================================================================

cand <- tribble(
  ~label,        ~type, ~lower, ~upper, ~integer, ~group,
  "PRDX(1)",     "rel",  0.50,   2.50,  FALSE, "production",   # upper replaced by N bound
  "PPDF(1)",     "abs", 10.00,  25.00,  FALSE, "production",   # optimum; bounds replaced by g$temp_priors
  "PPDF(2)",     "abs",  0.00,  50.00,  FALSE, "production",   # bounds replaced by g$temp_priors
  "PPDF(3)",     "rel",  0.25,   2.00,  FALSE, "production",
  "PPDF(4)",     "rel",  0.25,   1.00,  FALSE, "production",
  "PLTMRF",      "rel",  0.25,   2.00,  FALSE, "production",
  "FULCAN",      "rel",  0.50,   1.50,  FALSE, "production",
  "WSCOEF(1)",   "abs",  0.20,   0.50,  FALSE, "water stress",
  "WSCOEF(2)",   "abs",  6.00,  30.00,  FALSE, "water stress",
  "CLAYPG",      "abs",  2,      6,     TRUE,  "water stress",
  "FRTC(1)",     "rel",  0.50,   1.50,  FALSE, "allocation",
  "FRTC(2)",     "rel",  0.50,   1.50,  FALSE, "allocation",
  "FRTC(3)",     "rel",  0.50,   2.00,  FALSE, "allocation",
  "FRTC(4)",     "abs",  0.05,   0.30,  FALSE, "allocation",
  "FRTC(5)",     "abs",  0.05,   0.30,  FALSE, "allocation",
  "RTDTMP",      "abs", -2.00,   5.00,  FALSE, "allocation",
  "MRTFRAC",     "abs",  0.00,   0.20,  FALSE, "allocation",
  "PRAMN(1,1)",  "rel",  0.50,   1.00,  FALSE, "nitrogen",
  "PRAMN(1,2)",  "rel",  0.50,   1.00,  FALSE, "nitrogen",
  "PRAMX(1,1)",  "rel",  1.00,   1.50,  FALSE, "nitrogen",
  "PRAMX(1,2)",  "rel",  1.00,   1.50,  FALSE, "nitrogen",
  "PRBMN(1,1)",  "rel",  0.50,   1.00,  FALSE, "nitrogen",
  "PRBMX(1,1)",  "rel",  1.00,   1.50,  FALSE, "nitrogen",
  "BIOMAX",      "rel",  0.50,   1.50,  FALSE, "nitrogen",
  "NO3PREF(1)",  "abs",  0.10,   0.60,  FALSE, "nitrogen",
  "VLOSSP",      "abs",  0.00,   0.10,  FALSE, "nitrogen",
  "EFRGRN(1)",   "abs",  0.50,   0.90,  FALSE, "nitrogen",
  "SNFXMX(1)",   "rel",  0.50,   1.50,  FALSE, "nitrogen",
  "HIMAX",       "rel",  0.80,   1.20,  FALSE, "harvest",
  "HIWSF",       "abs",  0.00,   1.00,  FALSE, "harvest",
  "HIMON(1)",    "abs",  1,      3,     TRUE,  "harvest",
  "FSDETH(1)",   "abs",  0.00,   0.20,  FALSE, "death",
  "RDRJ",        "abs",  0.00,   0.10,  FALSE, "death",
  "RDRM",        "abs",  0.00,   0.10,  FALSE, "death",
  "CKMRSPMX(1)", "rel",  0.50,   1.50,  FALSE, "respiration",
  "CKMRSPMX(2)", "rel",  0.50,   1.50,  FALSE, "respiration",
  "CKMRSPMX(3)", "rel",  0.50,   1.50,  FALSE, "respiration",
  "NPP2CS(1)",   "abs",  0.75,   1.50,  FALSE, "respiration"
)

# =============================================================================
# 3. crop.100 AND DAYCENT HELPERS
# =============================================================================

esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
crop_block <- function(lines, crop_id) {
  s <- grep(paste0("^\\s*", esc(crop_id), "(\\s|$)"), lines)
  if (!length(s)) stop("Crop block '", crop_id, "' not found")
  s <- s[1]; hdr <- grep("^[A-Z][A-Z0-9]+\\s+\\S", lines); nx <- hdr[hdr > s][1]
  c(s, ifelse(is.na(nx), length(lines), nx - 1))
}
find_line <- function(lines, blk, label) {
  idx <- grep(paste0("^\\s*[-+0-9.eE]+\\s+'?", esc(label), "'?(\\s|$)"),
              lines[blk[1]:blk[2]], perl = TRUE)
  if (!length(idx)) integer(0) else blk[1] + idx[1] - 1
}
get_value <- function(lines, blk, label) {
  i <- find_line(lines, blk, label)
  if (!length(i)) NA_real_ else as.numeric(sub("^\\s*([-+0-9.eE]+).*$", "\\1", lines[i]))
}
write_crop100 <- function(vals, tpl, crop_id, int_labels, path) {
  blk <- crop_block(tpl, crop_id)
  if ("PPDF_GAP" %in% names(vals)) {        # PPDF(2) = PPDF(1) + gap (gap >= 8 C avoids zero-production crashes)
    p1 <- if ("PPDF(1)" %in% names(vals)) vals[["PPDF(1)"]] else get_value(tpl, blk, "PPDF(1)")
    vals <- c(vals[names(vals) != "PPDF_GAP"], "PPDF(2)" = p1 + vals[["PPDF_GAP"]])
  }
  if ("TSHIFT" %in% names(vals)) {          # thermal shift: move optimum and maximum together
    p1 <- get_value(tpl, blk, "PPDF(1)"); p2 <- get_value(tpl, blk, "PPDF(2)")
    vals <- c(vals[names(vals) != "TSHIFT"],
              "PPDF(1)" = p1 + vals[["TSHIFT"]], "PPDF(2)" = p2 + vals[["TSHIFT"]])
  }
  for (lab in names(vals)) {
    i <- find_line(tpl, blk, lab)
    if (!length(i)) stop("Parameter ", lab, " not found")
    s <- if (lab %in% int_labels) sprintf("%d", as.integer(round(vals[[lab]]))) else sprintf("%.6f", vals[[lab]])
    tpl[i] <- sub("^(\\s*)\\S+", paste0("\\1", s), tpl[i])
  }
  writeLines(tpl, path); invisible(path)
}

run_rm <- function(rm, vals) {
  d <- file.path(W$sites, rm)
  ok <- tryCatch({ write_crop100(vals, W$tpl, W$crop, W$int, file.path(d, "crop.100")); TRUE },
                 error = function(e) FALSE)
  if (!ok) return(NULL)
  old <- list.files(d, pattern = paste0("^", W$prefix, ".*\\.(bin|lis|csv|out)$"), full.names = TRUE)
  if (length(old)) file.remove(old)
  ec <- system(sprintf("cd %s && ./%s -s %s -n %s > daycent.log 2>&1",
                       shQuote(d), W$exe, W$sch, W$prefix))
  f <- file.path(d, paste0(W$prefix, ".csv"))
  if (ec != 0 || !file.exists(f)) return(NULL)
  h <- tryCatch(read.csv(f), error = function(e) NULL)
  unlink(file.path(d, paste0(W$prefix, ".bin")))
  if (is.null(h) || nrow(h) == 0 || !all(c("time", "cgrain") %in% names(h))) return(NULL)
  y <- tapply(h$cgrain, floor(h$time), sum)
  data.frame(RMNO = rm, Year = as.integer(names(y)), sim = as.numeric(y))
}
fert_edit <- function(rm, n_mult) {
  d <- file.path(W$sites, rm); f <- file.path(d, "fert.100"); o <- file.path(d, "fert.100.calib_orig")
  if (!file.exists(o)) file.copy(f, o)
  x <- readLines(o, warn = FALSE); i1 <- grep("FERAMT\\(1\\)", x)
  if (length(i1) != 1) return(FALSE)
  v <- as.numeric(sub("^\\s*([-+0-9.eE]+).*$", "\\1", x[i1]))
  x[i1] <- sub("^\\s*[-+0-9.eE]+", sprintf("%.3f", v * n_mult), x[i1]); writeLines(x, f); TRUE
}
fert_restore <- function(rm) {
  d <- file.path(W$sites, rm); o <- file.path(d, "fert.100.calib_orig")
  if (file.exists(o)) file.copy(o, file.path(d, "fert.100"), overwrite = TRUE)
}
run_rm_n <- function(rm, vals, n_mult) {
  on.exit(fert_restore(rm)); if (!fert_edit(rm, n_mult)) return(NULL); run_rm(rm, vals)
}

cl <- NULL
start_cluster <- function() {
  if (!is.null(cl)) return(invisible())
  cl <<- makeCluster(g$n_cores, type = "PSOCK")
  clusterExport(cl, c("W", "esc", "crop_block", "find_line", "get_value", "write_crop100", "run_rm",
                      "fert_edit", "fert_restore", "run_rm_n"), envir = globalenv())
}
stop_cluster <- function() { if (!is.null(cl)) stopCluster(cl); cl <<- NULL }
run_set <- function(vals, rms, n_mult = 1) {
  start_cluster()
  res <- if (n_mult == 1) parLapplyLB(cl, rms, function(rm, v) run_rm(rm, v), v = vals)
  else parLapplyLB(cl, rms, function(rm, v, m) run_rm_n(rm, v, m), v = vals, m = n_mult)
  out <- bind_rows(res); attr(out, "failed") <- rms[vapply(res, is.null, logical(1))]; out
}
run_cached <- function(id, vals, rms, n_mult = 1) {
  f <- file.path(C$out_dir, "runs", paste0(id, ".rds"))
  if (file.exists(f)) {
    rec <- readRDS(f)
    if (isTRUE(all.equal(rec$vals[names(vals)], vals, tolerance = 1e-8)) &&
        identical(rec$n_mult %||% 1, n_mult)) return(rec)
  }
  sim <- run_set(vals, rms, n_mult)
  rec <- list(vals = vals, n_mult = n_mult, sim = sim, failed = attr(sim, "failed"))
  saveRDS(rec, f); rec
}
sim_ok <- function(rec, n_rms) !is.null(rec$sim) && nrow(rec$sim) > 0 &&
  all(c("RMNO", "Year") %in% names(rec$sim)) && length(rec$failed) <= g$max_fail_frac * n_rms

# =============================================================================
# 4. METRICS
# =============================================================================

join_obs <- function(sim, obs) inner_join(obs, sim, by = c("RMNO", "Year"))
# Squared error is minimised by shrinking simulated variability toward r * obs SD when
# correlation is imperfect, so raw SSE rewards damped spatial spread and interannual
# variability. The variance-matched SSE rescales the simulation to the observed SD
# (alpha = 1) before scoring: SSE_vm = n[(bias)^2 + 2 var_obs (1 - r)], which only
# improves with bias and correlation. alpha is constrained separately in the likelihood.
lik_names <- c("sse_s", "sse_t", "sse_s_vm", "sse_t_vm", "alpha_s", "alpha_t",
               "bias_s", "r_s", "r_t", "n_s", "n_t")
lik_stats <- function(sim, obs_train) {
  j <- join_obs(sim, obs_train)
  rm <- j %>% group_by(RMNO) %>% summarise(o = mean(obs), s = mean(sim), .groups = "drop")
  a  <- j %>% group_by(RMNO) %>% mutate(oa = obs - mean(obs), sa = sim - mean(sim)) %>% ungroup()
  vm <- function(s, o, centre) {
    al <- sd(s) / sd(o)
    if (!is.finite(al) || al <= 0) return(c(sse = NA_real_, alpha = NA_real_))
    s2 <- if (centre) mean(s) + (s - mean(s)) / al else s / al   # anomalies already centred per RM
    c(sse = sum((s2 - o)^2), alpha = al)
  }
  vs <- vm(rm$s, rm$o, TRUE); vt <- vm(a$sa, a$oa, FALSE)
  c(sse_s = sum((rm$s - rm$o)^2), sse_t = sum((a$sa - a$oa)^2),
    sse_s_vm = vs[["sse"]], sse_t_vm = vt[["sse"]], alpha_s = vs[["alpha"]], alpha_t = vt[["alpha"]],
    bias_s = mean(rm$s) - mean(rm$o), r_s = cor(rm$s, rm$o), r_t = cor(a$sa, a$oa),
    n_s = nrow(rm), n_t = nrow(a))
}
kge3 <- function(s, o) { r <- cor(s, o); a <- sd(s) / sd(o); b <- mean(s) / mean(o)
c(KGE = 1 - sqrt((r - 1)^2 + (a - 1)^2 + (b - 1)^2), r = r, alpha = a, beta = b) }
nse <- function(s, o) 1 - sum((s - o)^2) / sum((o - mean(o))^2)
ccc <- function(s, o) { ms <- mean(s); mo <- mean(o)
2 * mean((s - ms) * (o - mo)) / (mean((s - ms)^2) + mean((o - mo)^2) + (ms - mo)^2) }
spatial_metrics <- function(j) {
  m <- j %>% group_by(RMNO) %>% summarise(o = mean(obs), s = mean(sim), .groups = "drop")
  if (nrow(m) < 5) return(setNames(rep(NA_real_, 6), c("sp_KGE","sp_r","sp_alpha","sp_beta","sp_NSE","sp_CCC")))
  k <- kge3(m$s, m$o)
  c(sp_KGE = k[["KGE"]], sp_r = k[["r"]], sp_alpha = k[["alpha"]], sp_beta = k[["beta"]],
    sp_NSE = nse(m$s, m$o), sp_CCC = ccc(m$s, m$o))
}
temporal_metrics <- function(j, prefix = "tm_") {
  a <- j %>% group_by(RMNO) %>% filter(n() >= 3) %>%
    mutate(oa = obs - mean(obs), sa = sim - mean(sim)) %>% ungroup()
  nm <- paste0(prefix, c("r", "alpha", "KGE", "NSE", "CCC"))
  if (nrow(a) < 20) return(setNames(rep(NA_real_, 5), nm))
  r <- cor(a$sa, a$oa); al <- sd(a$sa) / sd(a$oa)
  setNames(c(r, al, 1 - sqrt((r - 1)^2 + (al - 1)^2),
             1 - sum((a$sa - a$oa)^2) / sum(a$oa^2), ccc(a$sa, a$oa)), nm)
}
all_metrics <- function(j) c(spatial_metrics(j), temporal_metrics(j),
                             temporal_metrics(filter(j, !Year %in% C$worst_years), "tmx_"))
subsets <- function(j) list(
  train_insample = filter(j, RMNO %in% C$train_rms, period == "train"),
  test_spatial   = filter(j, split == "test", period == "train"),
  test_temporal  = filter(j, RMNO %in% C$train_rms, period == "test"),
  test_both      = filter(j, split == "test", period == "test"))

# =============================================================================
# 5. ONE CROP
# =============================================================================

ck  <- function(x) file.path(C$out_dir, paste0(x, ".rds"))
has <- function(p) p %in% phases
status <- function(msg) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%H:%M"), C$folder, msg))
  cat(sprintf("%s\t%s\t%s\n", format(Sys.time()), C$folder, msg),
      file = file.path(models_root, paste0("all_crops_status_", g$out_name, ".tsv")), append = TRUE)
}

run_crop <- function(folder, crop_id, pulse) {
  stop_cluster()
  model_dir <- file.path(models_root, folder)
  C <<- list(folder = folder, crop_id = crop_id, pulse = pulse,
             sites = file.path(model_dir, g$sites_subdir),
             tpl_path = file.path(model_dir, "Copy_Folder", "crop.100"),
             out_dir = file.path(model_dir, g$out_name))
  dir.create(file.path(C$out_dir, "runs"), recursive = TRUE, showWarnings = FALSE)
  # reuse DayCent runs and observation-free phases from an earlier output folder
  if (!is.null(g$seed_from)) {
    src <- file.path(model_dir, g$seed_from)
    if (dir.exists(src)) {   # copy anything not already present (also after an interrupted earlier run)
      file.copy(list.files(file.path(src, "runs"), full.names = TRUE), file.path(C$out_dir, "runs"),
                overwrite = FALSE)
      seed_rds <- c("smoke", "split", "morris", "sobol", "select")   # N bound redone: base temperatures may change
      for (f in seed_rds)
        if (file.exists(file.path(src, paste0(f, ".rds"))) && !file.exists(file.path(C$out_dir, paste0(f, ".rds"))))
          file.copy(file.path(src, paste0(f, ".rds")), file.path(C$out_dir, paste0(f, ".rds")))
      for (f in c("morris_indices.csv", "morris_selection.csv", "sobol_indices.csv",
                  "sobol_selection.csv", "prdx_n_response.csv", "rm_split.csv"))
        if (file.exists(file.path(src, f))) file.copy(file.path(src, f), file.path(C$out_dir, f), overwrite = FALSE)
      status(sprintf("seeded from %s (morris.rds %s, sobol.rds %s)", g$seed_from,
                     ifelse(file.exists(ck("morris")), "present", "MISSING"),
                     ifelse(file.exists(ck("sobol")), "present", "MISSING")))
    }
  }
  
  # ---- parameter table for this crop ----------------------------------------
  tpl <- readLines(C$tpl_path); blk <- crop_block(tpl, crop_id)
  pt <- cand %>% mutate(template = map_dbl(label, ~ if (.x == "TSHIFT") 0 else if (.x == "PPDF_GAP")
    get_value(tpl, blk, "PPDF(2)") - get_value(tpl, blk, "PPDF(1)") else get_value(tpl, blk, .x))) %>%
    filter(!is.na(template))
  pt <- pt %>% mutate(
    lower = ifelse(type == "rel", lower * template, lower),
    upper = ifelse(type == "rel", upper * template, upper),
    lower = pmin(lower, upper), upper = pmax(lower, upper))
  widened <- pt %>% filter(template < lower | template > upper) %>% pull(label)
  pt <- pt %>% mutate(lower = pmin(lower, template), upper = pmax(upper, template)) %>%
    filter(upper > lower)
  if (length(widened)) status(paste("ranges widened to include template:", paste(widened, collapse = ", ")))
  if (!pulse) pt <- filter(pt, label != "SNFXMX(1)")
  # ---- physiological temperature ranges replace generic bounds --------------------
  tp <- g$temp_priors[[folder]]
  if (is.null(tp) || any(is.na(unlist(tp)))) stop(folder, ": fill in g$temp_priors (PPDF(1), PPDF(2) ranges)")
  if (tp[["PPDF(2)"]][1] - tp[["PPDF(1)"]][2] < g$min_temp_gap)
    stop(folder, ": PPDF(2) lower bound must exceed PPDF(1) upper bound by >= ", g$min_temp_gap, " C")
  pt <- pt %>% filter(label != "PPDF_GAP") %>%
    mutate(lower = case_when(label == "PPDF(1)" ~ tp[["PPDF(1)"]][1], label == "PPDF(2)" ~ tp[["PPDF(2)"]][1], TRUE ~ lower),
           upper = case_when(label == "PPDF(1)" ~ tp[["PPDF(1)"]][2], label == "PPDF(2)" ~ tp[["PPDF(2)"]][2], TRUE ~ upper),
           # prior centre: CSU template if physiologically plausible, otherwise the range midpoint
           centre = ifelse(template >= lower & template <= upper, template, (lower + upper) / 2))
  off <- pt %>% filter(label %in% g$always_calibrate, template != centre)
  if (nrow(off)) status(paste("template outside physiological range, prior centred on midpoint:",
                              paste(sprintf("%s (template %g, range %g-%g)", off$label, off$template, off$lower, off$upper), collapse = "; ")))
  C$fixed <<- numeric(0)
  C$pt <<- pt; C$int <<- pt$label[pt$integer]
  if (!"centre" %in% names(pt)) pt$centre <- pt$template
  C$base <<- setNames(pt$centre, pt$label)
  C$csu  <<- setNames(pt$template, pt$label)
  write_csv(pt, file.path(C$out_dir, "parameter_table.csv"))
  # fixed values are always written, so cached runs made with other values are never reused
  make_params <- function(x = NULL, csu = FALSE) { v <- c(if (csu) C$csu else C$base, C$fixed); if (!is.null(x)) v[names(x)] <- x
  v[C$int] <- round(v[C$int]); v }
  
  W <<- list(sites = C$sites, exe = g$exe, sch = g$schedule, prefix = g$out_prefix,
             crop = crop_id, tpl = tpl, int = C$int)
  
  # ---- observations -----------------------------------------------------------
  rms <- list.dirs(C$sites, full.names = FALSE, recursive = FALSE); rms <- rms[nchar(rms) > 0]
  obs_all <- map_dfr(rms, function(rm) {
    f <- file.path(C$sites, rm, g$obs_file); if (!file.exists(f)) return(NULL)
    d <- read.csv(f); tibble(RMNO = rm, Year = as.integer(d$Year), obs = as.numeric(d$cgrain))
  }) %>% filter(!is.na(obs))
  status(sprintf("%d RMs with obs, %d RM-years", n_distinct(obs_all$RMNO), nrow(obs_all)))
  
  # ---- smoke (CSU defaults, all RMs) ----------------------------------------
  if (has("smoke") && !file.exists(ck("smoke"))) {
    t0 <- Sys.time(); sim <- run_set(make_params(csu = TRUE), unique(obs_all$RMNO))   # true CSU defaults
    saveRDS(list(sim = sim, failed = attr(sim, "failed"),
                 secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), ck("smoke"))
    status(sprintf("smoke: %.1f s, %d failed RMs", readRDS(ck("smoke"))$secs, length(attr(sim, "failed"))))
  }
  
  # ---- split -------------------------------------------------------------------
  if (has("split") && !file.exists(ck("split"))) {
    r_ids <- unique(obs_all$RMNO)
    pts <- RM_PTS
    cand_f <- setdiff(names(pts), attr(pts, "sf_column"))
    hits <- sapply(cand_f, function(f) sum(as.character(suppressWarnings(as.integer(as.numeric(pts[[f]])))) %in% r_ids))
    fld <- cand_f[which.max(hits)]
    if (!all(sf::st_geometry_type(pts) == "POINT")) pts <- sf::st_centroid(pts)
    cc <- sf::st_coordinates(pts)
    xy <- tibble(RMNO = as.character(suppressWarnings(as.integer(as.numeric(pts[[fld]])))),
                 x = cc[, 1], y = cc[, 2]) %>% distinct(RMNO, .keep_all = TRUE) %>% filter(RMNO %in% r_ids)
    set.seed(g$seed)
    xy$block <- kmeans(scale(xy[, c("x", "y")]), centers = g$n_blocks, nstart = 50)$cluster
    xy <- left_join(xy, obs_all %>% group_by(RMNO) %>% summarise(m = mean(obs)), by = "RMNO")
    combos <- unlist(lapply(1:3, function(k) combn(sort(unique(xy$block)), k, simplify = FALSE)), recursive = FALSE)
    sc <- map_dfr(combos, function(cb) {
      te <- xy$m[xy$block %in% cb]; tr <- xy$m[!xy$block %in% cb]
      tibble(blocks = paste(cb, collapse = ","), frac = length(te) / nrow(xy),
             score = abs(mean(te) - mean(tr)) / sd(xy$m) + abs(sd(te) - sd(tr)) / sd(xy$m),
             inside = min(te) >= min(tr) & max(te) <= max(tr))
    }) %>% filter(abs(frac - g$test_block_frac) <= g$test_frac_tol) %>% arrange(desc(inside), score)
    tb <- as.integer(strsplit(sc$blocks[1], ",")[[1]])
    xy$split <- ifelse(xy$block %in% tb, "test", "train"); xy$m <- NULL
    oa <- obs_all %>% left_join(select(xy, RMNO, block, split), by = "RMNO") %>%
      mutate(period = ifelse(Year >= g$test_year_start, "test", "train"))
    tr_rms <- oa %>% filter(split == "train", period == "train") %>% count(RMNO) %>%
      filter(n >= g$min_years_rm) %>% pull(RMNO)
    write_csv(xy, file.path(C$out_dir, "rm_split.csv"))
    saveRDS(list(xy = xy, obs_all = oa, obs_train = filter(oa, RMNO %in% tr_rms, period == "train"),
                 train_rms = tr_rms, test_blocks = tb), ck("split"))
    status(sprintf("split: test blocks %s, %d train / %d test RMs",
                   paste(tb, collapse = ","), length(tr_rms), sum(xy$split == "test")))
  }
  sp <- readRDS(ck("split"))
  obs_all <- sp$obs_all; obs_train <- sp$obs_train
  C$train_rms <<- sp$train_rms
  C$worst_years <<- obs_all %>% group_by(Year) %>% summarise(m = mean(obs)) %>%
    slice_min(m, n = g$n_worst_years) %>% pull(Year)
  
  # Screening outputs (simulated only — observations are NOT used, so screening
  # stays independent of the calibration fit):
  #   1 mean yield              2 spread of RM means        3 mean within-RM SD
  #   4 spatial pattern change  = 1 - cor(RM means, CSU-default RM means)
  #   5 temporal pattern change = 1 - cor(anomalies, CSU-default anomalies)
  resp_names <- c("mean_yield", "spatial_spread", "temporal_variability",
                  "spatial_pattern", "temporal_pattern")
  ref <- join_obs(readRDS(ck("smoke"))$sim, obs_train) %>% select(RMNO, Year, ref = sim) %>%
    group_by(RMNO) %>% mutate(ref_a = ref - mean(ref)) %>% ungroup()
  ref_rm <- ref %>% group_by(RMNO) %>% summarise(ref_m = mean(ref), .groups = "drop")
  resp <- function(rec) {
    if (!sim_ok(rec, length(C$train_rms))) return(rep(NA_real_, 5))
    j <- join_obs(rec$sim, obs_train)
    rmm <- j %>% group_by(RMNO) %>% summarise(s = mean(sim), sd = sd(sim), .groups = "drop")
    sp <- inner_join(rmm, ref_rm, by = "RMNO")
    tp <- j %>% group_by(RMNO) %>% mutate(a = sim - mean(sim)) %>% ungroup() %>%
      inner_join(ref, by = c("RMNO", "Year"))
    pat <- function(x, y) if (sd(x) > 0 && sd(y) > 0) 1 - cor(x, y) else 0
    c(mean(j$sim), sd(rmm$s), mean(rmm$sd, na.rm = TRUE), pat(sp$s, sp$ref_m), pat(tp$a, tp$ref_a))
  }
  impute <- function(Y, lab) {
    nf <- sum(!complete.cases(Y))
    if (nf) status(sprintf("%s: %d of %d runs failed (%.1f%%), mean-imputed", lab, nf, nrow(Y), 100 * nf / nrow(Y)))
    for (k in seq_len(ncol(Y))) Y[is.na(Y[, k]), k] <- mean(Y[, k], na.rm = TRUE); Y
  }
  
  # ---- Morris -------------------------------------------------------------------
  if (has("morris") && !file.exists(ck("morris"))) {
    set.seed(g$seed)
    mo <- morris(model = NULL, factors = pt$label, r = g$morris_r,
                 design = list(type = "oat", levels = g$morris_levels, grid.jump = g$morris_levels / 2),
                 binf = pt$lower, bsup = pt$upper, scale = TRUE)
    Y <- t(sapply(seq_len(nrow(mo$X)), function(i)
      resp(run_cached(sprintf("morris_%04d", i), make_params(setNames(mo$X[i, ], colnames(mo$X))), C$train_rms))))
    fail_frac <- mean(!complete.cases(Y)); Y <- impute(Y, "Morris")
    ee <- map_dfr(seq_len(ncol(Y)), function(k) { m_k <- mo; tell(m_k, Y[, k])
    tibble(output = k, output_name = resp_names[k], label = colnames(m_k$ee), mu_star = colMeans(abs(m_k$ee)), sigma = apply(m_k$ee, 2, sd))
    }) %>% group_by(output) %>% mutate(rel = mu_star / max(mu_star)) %>% ungroup()
    sel <- ee %>% group_by(label) %>% summarise(max_rel = max(rel)) %>%
      left_join(select(pt, label, integer), by = "label") %>%
      mutate(passed = max_rel >= g$morris_thresh & !integer) %>% arrange(desc(max_rel))
    write_csv(ee, file.path(C$out_dir, "morris_indices.csv")); write_csv(sel, file.path(C$out_dir, "morris_selection.csv"))
    saveRDS(list(sel = sel, passed = sel$label[sel$passed], fail_frac = fail_frac), ck("morris"))
    status(sprintf("Morris: %d passed, %.1f%% runs failed", sum(sel$passed), 100 * fail_frac))
  }
  if (file.exists(ck("morris")) && readRDS(ck("morris"))$fail_frac > g$max_fail_frac)   # checked on every run
    stop(folder, ": Morris failure rate > 5% — check parameter ranges before continuing")
  
  # ---- Sobol -------------------------------------------------------------------
  if (has("sobol") && !file.exists(ck("sobol"))) {
    passed <- readRDS(ck("morris"))$passed
    ps <- pt %>% filter(label %in% passed) %>% slice(match(passed, label)); k <- length(passed)
    set.seed(g$seed + 10)
    draw <- function(n) { U <- randomLHS(n, k)
    X <- sweep(sweep(U, 2, ps$upper - ps$lower, "*"), 2, ps$lower, "+"); colnames(X) <- passed
    as.data.frame(X, check.names = FALSE) }
    so <- soboljansen(model = NULL, X1 = draw(g$sobol_N), X2 = draw(g$sobol_N), nboot = g$sobol_nboot)
    Xs <- as.matrix(so$X)
    Y <- t(sapply(seq_len(nrow(Xs)), function(i)
      resp(run_cached(sprintf("sobol_%05d", i), make_params(setNames(Xs[i, ], passed)), C$train_rms))))
    Y <- impute(Y, "Sobol")
    tab <- map_dfr(seq_len(ncol(Y)), function(o) { s_o <- so; tell(s_o, Y[, o])
    tibble(output = o, output_name = resp_names[o], label = passed, ST = s_o$T[, "original"],
           ST_lo = s_o$T[, "min. c.i."], ST_hi = s_o$T[, "max. c.i."]) })
    write_csv(tab, file.path(C$out_dir, "sobol_indices.csv"))
    saveRDS(list(tab = tab), ck("sobol"))
    status(sprintf("Sobol done (%d parameters, %d runs)", k, nrow(Xs)))
  }
  
  # ---- automatic selection ------------------------------------------------------
  if (has("select")) {
    tab <- readRDS(ck("sobol"))$tab
    st <- tab %>% group_by(label) %>% slice_max(ST_lo, n = 1, with_ties = FALSE) %>% ungroup() %>%
      arrange(desc(ST)) %>%
      mutate(eligible = ST_lo >= g$sobol_ST_thresh & !label %in% g$never_calibrate &
               !label %in% C$int & label %in% C$pt$label)   # e.g. PPDF_GAP from older screening is not a parameter now
    calib <- head(st$label[st$eligible], g$max_calib)
    forced <- intersect(unique(c(g$always_calibrate, g$force_calibrate[[folder]] %||% character(0))), C$pt$label)
    if (length(forced)) {
      calib <- c(forced, setdiff(calib, forced))
      partners <- unlist(lapply(g$paired, function(pr) if (any(pr %in% forced)) pr))
      keep <- unique(c(forced, intersect(partners, C$pt$label)))
      while (length(calib) > g$max_calib && length(setdiff(calib, keep))) calib <- setdiff(calib, tail(setdiff(calib, keep), 1))
      status(paste("forced into calibration:", paste(forced, collapse = ", ")))
    }
    # paired parameters: if one member is selected, add its partner; keep within
    # max_calib by dropping the lowest-ranked unpaired parameter(s)
    for (pr in g$paired) {
      pr <- pr[pr %in% C$pt$label]
      if (length(pr) == 2 && any(pr %in% calib) && !all(pr %in% calib)) {
        add <- setdiff(pr, calib)
        calib <- c(calib, add)
        while (length(calib) > g$max_calib) {
          drop_cand <- rev(setdiff(calib, unlist(g$paired)))
          if (!length(drop_cand)) break
          calib <- setdiff(calib, drop_cand[1])
        }
        status(paste("paired rule: added", paste(add, collapse = ", ")))
      }
    }
    st$selected <- st$label %in% calib
    st$forced <- st$label %in% forced
    st$added_by_pair <- st$label %in% setdiff(calib, head(st$label[st$eligible], g$max_calib))
    write_csv(st, file.path(C$out_dir, "sobol_selection.csv"))
    saveRDS(list(calib = calib, table = st), ck("select"))
    status(paste("selected:", paste(calib, collapse = ", ")))
  }
  calib <- readRDS(ck("select"))$calib
  # if the selected set changed, the cached design/emulator/posterior no longer apply
  if (file.exists(ck("design")) && !all(calib %in% names(readRDS(ck("design"))$des))) {
    status("calibrated set changed: rebuilding design (cached runs reused), GP, MCMC, check, post")
    unlink(c(ck("design"), ck("gp"), ck("mcmc_final"), ck("check"), ck("refine"),
             file.path(C$out_dir, c("evaluation_metrics.csv", "n_response.csv"))))
  }
  if (length(calib) < 1) { status("no parameter met the selection rule — defaults retained"); return(invisible()) }
  
  # ---- N-response bound for PRDX(1): interpolated crossing --------------------------
  if (has("nbound") && "PRDX(1)" %in% calib && !pulse && !file.exists(ck("nbound"))) {
    set.seed(g$seed + 20)
    nb_rms <- sp$xy %>% filter(RMNO %in% C$train_rms) %>% group_by(block) %>%
      slice_sample(n = 1) %>% ungroup() %>% slice_head(n = g$nresp_rms) %>% pull(RMNO)
    t0 <- C$base[["PRDX(1)"]]
    nb <- map_dfr(g$nresp_grid, function(m) {
      v <- make_params(c("PRDX(1)" = t0 * m))
      y1 <- run_cached(sprintf("nb_%03d_x1", round(100 * m)), v, nb_rms, 1)
      y2 <- run_cached(sprintf("nb_%03d_x2", round(100 * m)), v, nb_rms, 2)
      mean_y <- function(r) mean(join_obs(r$sim, obs_all)$sim)
      tibble(mult = m, prdx = t0 * m, y1 = mean_y(y1), y2 = mean_y(y2))
    }) %>% mutate(pct_at_1x = y1 / y2) %>% arrange(prdx)
    
    ub <- if (all(nb$pct_at_1x >= g$nresp_min)) {
      max(nb$prdx)                                   # never crosses within the grid
    } else if (nb$pct_at_1x[1] < g$nresp_min) {
      t0                                             # template already below threshold
    } else {
      k <- which(nb$pct_at_1x < g$nresp_min)[1]      # first grid point below threshold
      approx(nb$pct_at_1x[(k - 1):k], nb$prdx[(k - 1):k], xout = g$nresp_min)$y
    }
    nb$upper_bound <- ub
    write_csv(nb, file.path(C$out_dir, "prdx_n_response.csv"))
    saveRDS(list(table = nb, upper = ub), ck("nbound"))
    status(sprintf("PRDX(1) upper bound from N response: %.3f (template %.2f, %.2fx)", ub, t0, ub / t0))
  }
  if (file.exists(ck("nbound"))) {
    ub <- readRDS(ck("nbound"))$upper
    C$pt <<- C$pt %>% mutate(upper = ifelse(label == "PRDX(1)", ub, upper))
  }
  
  # ---- calibration bounds and prior ---------------------------------------------
  pc <- C$pt %>% filter(label %in% calib) %>% slice(match(calib, label))
  lo <- setNames(pc$lower, calib); hi <- setNames(pc$upper, calib); d <- length(calib)
  mu <- pmin(pmax(setNames(pc$centre, calib), lo), hi)   # template, or range midpoint if template implausible
  to_x <- function(U) sweep(sweep(as.matrix(U), 2, hi - lo, "*"), 2, lo, "+")
  to_u <- function(X) sweep(sweep(as.matrix(X), 2, lo, "-"), 2, hi - lo, "/")
  gpn <- paste0("p", seq_len(d))
  prior <- if (g$prior_type == "tnorm")
    createTruncatedNormalPrior(mean = mu, sd = g$prior_sd_frac * (hi - lo), lower = lo, upper = hi)
  else createUniformPrior(lower = lo, upper = hi)
  
  run_points <- function(Xm, tag) map_dfr(seq_len(nrow(Xm)), function(i) {
    rec <- run_cached(sprintf("%s_%04d", tag, i), make_params(setNames(Xm[i, ], calib)), C$train_rms)
    st <- if (!sim_ok(rec, length(C$train_rms))) setNames(rep(NA_real_, length(lik_names)), lik_names)
    else lik_stats(rec$sim, obs_train)
    as_tibble(as.list(c(setNames(Xm[i, ], calib), st)))
  })
  
  # ---- design -------------------------------------------------------------------
  if (has("design") && file.exists(ck("design")) && !all(lik_names %in% names(readRDS(ck("design"))$des))) {
    status("design table lacks current likelihood statistics: rebuilding from cached runs")
    unlink(c(ck("design"), ck("gp"), ck("mcmc_final"), ck("check"), ck("refine")))
  }
  if (has("design") && !file.exists(ck("design"))) {
    n_des <- max(g$design_min, g$design_per_dim * d); set.seed(g$seed + 1)
    U_des <- tryCatch(maximinLHS(n_des, d), error = function(e) randomLHS(n_des, d))
    U_val <- randomLHS(g$n_val, d)
    des <- run_points(to_x(U_des), "design"); val <- run_points(to_x(U_val), "val")
    saveRDS(list(des = des, val = val), ck("design"))
    status(sprintf("design: %d + %d runs, %d failed", n_des, g$n_val,
                   sum(is.na(des$sse_s)) + sum(is.na(val$sse_s))))
  }
  dz <- if (file.exists(ck("design"))) readRDS(ck("design"))
  
  # emulator targets
  gp_ok <- function(d) switch(g$lik,
                              sse = filter(d, !is.na(sse_s), !is.na(sse_t)),
                              vm  = filter(d, !is.na(sse_s_vm), !is.na(sse_t_vm), !is.na(alpha_s), !is.na(alpha_t)),
                              vmc = filter(d, !is.na(bias_s), !is.na(r_s), !is.na(r_t), !is.na(alpha_s), !is.na(alpha_t),
                                           abs(r_s) < 1, abs(r_t) < 1))
  gp_targets <- function(d) switch(g$lik,
                                   sse = list(s = log(d$sse_s), t = log(d$sse_t)),
                                   vm  = list(s = log(d$sse_s_vm), t = log(d$sse_t_vm), as = log(d$alpha_s), at = log(d$alpha_t)),
                                   vmc = list(b = d$bias_s, zs = atanh(d$r_s), zt = atanh(d$r_t),
                                              as = log(d$alpha_s), at = log(d$alpha_t)))
  gp_rows <- c(s = "spatial", t = "temporal", as = "alpha_spatial", at = "alpha_temporal",
               b = "bias_spatial", zs = "r_spatial", zt = "r_temporal")
  fit_gps <- function(des) {
    des <- gp_ok(des); U <- as.data.frame(to_u(des[, calib])); names(U) <- gpn
    map(gp_targets(des), ~ km(~., design = U, response = .x, covtype = "matern5_2",
                              nugget.estim = TRUE, control = list(trace = FALSE)))
  }
  gp_val <- function(gps, val) {
    val <- gp_ok(val); U <- as.data.frame(to_u(val[, calib])); names(U) <- gpn
    y <- gp_targets(val)
    f <- function(m, yy) { p <- predict(m, newdata = U, type = "UK", checkNames = FALSE)
    c(R2 = 1 - sum((yy - p$mean)^2) / sum((yy - mean(yy))^2), RMSE = sqrt(mean((yy - p$mean)^2)),
      cover95 = mean(abs(yy - p$mean) <= 1.96 * p$sd)) }
    out <- do.call(rbind, map(names(gps), ~ f(gps[[.x]], y[[.x]])))
    rownames(out) <- gp_rows[names(gps)]; out
  }
  
  # ---- likelihood and emulator helpers (shared by MCMC and refinement) ---------------
  if (!is.null(dz)) {
    des_ok <- filter(dz$des, !is.na(n_s)); n_s <- des_ok$n_s[1]; n_t <- des_ok$n_t[1]
  }
  ns_e <- g$neff_spatial_frac * (if (exists("n_s")) n_s else NA)
  nt_e <- n_distinct(obs_train$Year)    # years are the independent units (errors shared across RMs within a year)
  ll_vals <- function(ls, lt, las = 0, lat = 0) {
    out <- -0.5 * ns_e * (ls - log(n_s)) - 0.5 * nt_e * (lt - log(n_t))
    if (g$lik %in% c("vm", "vmc"))   # lognormal constraint: simulated variability should match observed
      out <- out - 0.5 * (las / g$alpha_sd)^2 - 0.5 * (lat / g$alpha_sd)^2
    out
  }
  # exact variance-matched SSE from the components; observed variances are constants
  rmo <- obs_train %>% group_by(RMNO) %>% summarise(o = mean(obs), .groups = "drop")
  ano <- obs_train %>% group_by(RMNO) %>% mutate(oa = obs - mean(obs)) %>% ungroup()
  v_s <- mean((rmo$o - mean(rmo$o))^2); v_t <- mean(ano$oa^2)
  sse_from <- function(b, zs, zt) list(
    s = log(n_s * (b^2 + 2 * v_s * (1 - tanh(zs)))),
    t = log(2 * n_t * v_t * (1 - tanh(zt))))
  ll_true <- function(runs) { y <- gp_targets(runs)
  if (g$lik == "vmc") { q <- sse_from(y$b, y$zs, y$zt); ll_vals(q$s, q$t, y$as, y$at) }
  else ll_vals(y$s, y$t, y$as %||% 0, y$at %||% 0) }
  # emulator = global GPs everywhere, local GPs inside the refinement box (if fitted)
  emu_pred <- function(emu, X, se = FALSE) {
    X <- matrix(X, ncol = d)
    pr <- function(gps, Xs) { U <- as.data.frame(to_u(Xs)); names(U) <- gpn
    map(gps, function(m) { p <- predict(m, newdata = U, type = "UK", checkNames = FALSE, se.compute = se)
    list(m = p$mean, s = p$sd) }) }
    out <- pr(emu$global, X)
    if (!is.null(emu$local)) {
      inb <- apply(X, 1, function(x) all(x >= emu$box$lo & x <= emu$box$hi))
      if (any(inb)) { pl <- pr(emu$local, X[inb, , drop = FALSE])
      for (k in names(out)) { out[[k]]$m[inb] <- pl[[k]]$m; if (se) out[[k]]$s[inb] <- pl[[k]]$s } }
    }
    out
  }
  ll_emu <- function(emu) function(x) { p <- emu_pred(emu, x)
  if (g$lik == "vmc") { q <- sse_from(p$b$m, p$zs$m, p$zt$m); ll_vals(q$s, q$t, p$as$m, p$at$m) }
  else ll_vals(p$s$m, p$t$m, p$as$m %||% 0, p$at$m %||% 0) }
  # check that the component formula reproduces the stored variance-matched SSE
  if (g$lik == "vmc" && !is.null(dz)) {
    dd <- gp_ok(dz$des); q <- sse_from(dd$bias_s, atanh(dd$r_s), atanh(dd$r_t))
    err <- max(abs(exp(q$s) / dd$sse_s_vm - 1), abs(exp(q$t) / dd$sse_t_vm - 1), na.rm = TRUE)
    status(sprintf("vmc: max relative error of component SSE vs stored SSE = %.2e", err))
    if (err > 0.02) status("WARNING: component SSE differs from stored SSE by > 2% (RMs dropped from some runs?)")
  }
  emu_val <- function(emu, runs) {
    runs <- gp_ok(runs); y <- gp_targets(runs); p <- emu_pred(emu, as.matrix(runs[, calib]), se = TRUE)
    out <- do.call(rbind, map(names(y), function(k)
      c(R2 = 1 - sum((y[[k]] - p[[k]]$m)^2) / sum((y[[k]] - mean(y[[k]]))^2),
        RMSE = sqrt(mean((y[[k]] - p[[k]]$m)^2)), cover95 = mean(abs(y[[k]] - p[[k]]$m) <= 1.96 * p[[k]]$s))))
    rownames(out) <- gp_rows[names(y)]; out
  }
  emu_ok <- function(vr_post, vr_hold)
    all(vr_post[, "RMSE"] <= g$refine$rmse_tol * vr_hold[rownames(vr_post), "RMSE"]) &&
    all(vr_post[, "cover95"] >= g$refine$cover_min)
  run_dream <- function(ll, pr, seed) {
    set.seed(seed)
    o <- runMCMC(createBayesianSetup(ll, prior = pr, names = calib), sampler = "DREAMzs",
                 settings = list(iterations = g$mcmc_iter, nrChains = g$mcmc_chains, message = FALSE))
    ml <- as.mcmc.list(lapply(getSample(o, coda = TRUE, parametersOnly = TRUE), function(ch) {
      n <- nrow(ch); mcmc(as.matrix(ch)[(floor(n * g$burn_frac) + 1):n, , drop = FALSE]) }))
    rhat <- gelman.diag(ml, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1]
    S <- do.call(rbind, lapply(ml, as.matrix)); colnames(S) <- calib
    list(S = S, rhat = setNames(rhat, calib))
  }
  
  # ---- GP -----------------------------------------------------------------------
  if (has("gp") && !file.exists(ck("gp"))) {
    gps <- fit_gps(dz$des); vr <- gp_val(gps, dz$val)
    saveRDS(list(gps = gps, val = vr), ck("gp"))
    write.csv(round(vr, 3), file.path(C$out_dir, "gp_validation.csv"))
    status(paste("GP R2:", paste(sprintf("%s %.3f", rownames(vr), vr[, "R2"]), collapse = ", ")))
    if (any(vr[, "R2"] < 0.85)) status("WARNING: emulator R2 < 0.85")
  }
  
  # ---- DEoptim (posterior mode) + DREAM(ZS) ----------------------------------------
  lik_cfg <- list(lik = g$lik, alpha_sd = g$alpha_sd, nt_e = "years", temp_priors = g$temp_priors[[folder]])
  if (file.exists(ck("mcmc_final")) && !identical(readRDS(ck("mcmc_final"))$lik_cfg, lik_cfg)) {
    status("likelihood or fixed-parameter settings changed: rerunning MCMC, check, post, mgmt")
    unlink(c(ck("mcmc_final"), ck("check"), ck("refine"), file.path(C$out_dir,
                                                                    c("evaluation_metrics.csv", "posterior_summary.csv", "n_response.csv",
                                                                      "sim_vs_obs_posterior_median.csv"))))
  }
  if (has("mcmc") && !file.exists(ck("mcmc_final"))) {
    gps <- readRDS(ck("gp"))$gps
    ll <- ll_emu(list(global = gps, local = NULL, box = NULL))
    
    # DEoptim on the emulated log-posterior (same likelihood and prior as DREAM(ZS))
    neg_post <- function(x) { lp <- prior$density(x); if (!is.finite(lp)) return(1e10); -(ll(x) + lp) }
    set.seed(g$seed + 30)
    de <- DEoptim(neg_post, lower = unname(lo), upper = unname(hi),
                  control = DEoptim.control(NP = g$de_np_per_dim * d, itermax = g$de_itermax,
                                            trace = FALSE))
    de_best <- setNames(de$optim$bestmem, calib)
    status(paste("DEoptim mode:", paste(sprintf("%s=%.4g", calib, de_best), collapse = " ")))
    
    # DREAM(ZS): independent chains started from the prior (DE result used only as a check)
    mc <- run_dream(ll, prior, g$seed + 2); S <- mc$S; rhat <- mc$rhat
    saveRDS(list(S = S, rhat = setNames(rhat, calib), de_best = de_best,
                 de_value = de$optim$bestval, lik_cfg = lik_cfg), ck("mcmc_final"))
    status(paste("R-hat:", paste(sprintf("%s=%.3f", calib, rhat), collapse = " ")))
    if (any(rhat > 1.1)) status("WARNING: R-hat > 1.1")
    # mode vs posterior: distance in posterior-SD units
    z <- (de_best - apply(S, 2, median)) / apply(S, 2, sd)
    status(paste("DE mode vs posterior median (SD units):",
                 paste(sprintf("%s=%.2f", calib, z), collapse = " ")))
    if (any(abs(z) > 2)) status("WARNING: DE mode > 2 posterior SD from median — check for multimodality")
  }
  if (!file.exists(ck("mcmc_final"))) return(invisible())
  fin <- readRDS(ck("mcmc_final")); S <- fin$S; med <- apply(S, 2, median)
  de_best <- fin$de_best %||% setNames(rep(NA_real_, length(calib)), calib)
  
  # ---- emulator check near the posterior ---------------------------------------------
  if (has("check") && !file.exists(ck("check"))) {
    set.seed(g$seed + 3)
    chk <- run_points(S[sample.int(nrow(S), g$n_check), , drop = FALSE], "check")
    vr <- emu_val(list(global = readRDS(ck("gp"))$gps, local = NULL, box = NULL), chk)
    saveRDS(list(runs = chk, val = vr), ck("check"))
    write.csv(round(vr, 3), file.path(C$out_dir, "gp_check_posterior_region.csv"))
    status(paste("posterior-region emulator RMSE:",
                 paste(sprintf("%s %.3f", rownames(vr), vr[, "RMSE"]), collapse = ", ")))
  }
  
  # ---- emulator refinement near the posterior (automatic, only if the check fails) ----
  if (has("refine") && file.exists(ck("check")) && !file.exists(ck("refine"))) {
    rf <- g$refine; gp0 <- readRDS(ck("gp")); hold <- gp0$val; chk0 <- readRDS(ck("check"))
    emu <- list(global = gp0$gps, local = NULL, box = NULL)
    S_cur <- S; rhat_cur <- fin$rhat; vr_post <- chk0$val
    pool <- bind_rows(dz$des, dz$val, chk0$runs)
    ratio <- function(v) max(v[, "RMSE"] / hold[rownames(v), "RMSE"])
    conv <- emu_ok(vr_post, hold)
    log_tab <- tibble(round = 0, n_pool = nrow(pool), rmse_ratio_max = ratio(vr_post),
                      cover_min = min(vr_post[, "cover95"]), rhat_max = max(rhat_cur), converged = conv, kept = TRUE)
    status(sprintf("refine: posterior-region check %s (max RMSE ratio %.2f, min coverage %.2f)",
                   if (conv) "passed, no refinement needed" else "failed", ratio(vr_post), min(vr_post[, "cover95"])))
    box_of <- function(S) {
      q1 <- apply(S, 2, quantile, rf$box_q[1]); q2 <- apply(S, 2, quantile, rf$box_q[2])
      w <- pmax(q2 - q1, 1e-6 * (hi - lo))
      list(lo = pmax(lo, q1 - rf$box_widen * w), hi = pmin(hi, q2 + rf$box_widen * w))
    }
    in_box <- function(n, b) sweep(sweep(matrix(runif(n * d), n, d), 2, b$hi - b$lo, "*"), 2, b$lo, "+")
    r <- 0
    while (!conv && r < rf$max_rounds) {
      r <- r + 1; box <- box_of(S_cur); set.seed(g$seed + 100 + r)
      Xn <- rbind(S_cur[sample.int(nrow(S_cur), rf$n_post_pts), , drop = FALSE],
                  tryCatch(maximinLHS(rf$n_box_pts, d), error = function(e) randomLHS(rf$n_box_pts, d)) %>%
                    sweep(2, box$hi - box$lo, "*") %>% sweep(2, box$lo, "+"))
      colnames(Xn) <- calib
      pool <- bind_rows(pool, run_points(Xn, sprintf("refine%d", r)))
      ok  <- gp_ok(pool)
      inb <- apply(as.matrix(ok[, calib]), 1, function(x) all(x >= box$lo & x <= box$hi))
      emu_new <- list(global = fit_gps(ok),
                      local  = if (sum(inb) >= rf$min_local_per_dim * d) fit_gps(ok[inb, ]) else NULL, box = box)
      bprior <- createPrior(density = function(x) prior$density(x),
                            sampler = function(n = 1) in_box(n, box), lower = box$lo, upper = box$hi)
      mc <- run_dream(ll_emu(emu_new), bprior, g$seed + 200 + r)
      if (max(mc$rhat) > rf$rhat_max) {
        status(sprintf("refine round %d: R-hat %.2f > %.2f, round discarded", r, max(mc$rhat), rf$rhat_max))
        log_tab <- bind_rows(log_tab, tibble(round = r, n_pool = nrow(pool), rmse_ratio_max = NA,
                                             cover_min = NA, rhat_max = max(mc$rhat), converged = FALSE, kept = FALSE))
        break
      }
      emu <- emu_new; S_cur <- mc$S; rhat_cur <- mc$rhat
      set.seed(g$seed + 300 + r)
      chk <- run_points(S_cur[sample.int(nrow(S_cur), rf$n_post_pts), , drop = FALSE], sprintf("rcheck%d", r))
      pool <- bind_rows(pool, chk)
      vr_post <- emu_val(emu, chk); conv <- emu_ok(vr_post, hold)
      log_tab <- bind_rows(log_tab, tibble(round = r, n_pool = nrow(pool), rmse_ratio_max = ratio(vr_post),
                                           cover_min = min(vr_post[, "cover95"]), rhat_max = max(rhat_cur),
                                           converged = conv, kept = TRUE))
      status(sprintf("refine round %d: %d runs in box (local GP %s), max RMSE ratio %.2f, min coverage %.2f, R-hat max %.3f%s",
                     r, sum(inb), if (is.null(emu$local)) "not fitted" else "fitted", ratio(vr_post),
                     min(vr_post[, "cover95"]), max(rhat_cur), if (conv) " - converged" else ""))
    }
    if (!conv) status("refine: emulator did not pass the posterior-region check; DayCent pool check decides the point estimate")
    # DayCent-verified point estimate: refined median unless a DayCent run there is clearly worse
    pl <- gp_ok(pool) %>% distinct(across(all_of(calib)), .keep_all = TRUE)
    pl$ll <- ll_true(pl); ll_best <- max(pl$ll)
    top <- filter(pl, ll >= ll_best - rf$pool_delta_ll)
    if (nrow(top) < rf$pool_min_top) top <- slice_max(pl, ll, n = rf$pool_min_top)
    Xt <- as.matrix(top[, calib]); z <- sweep(sweep(Xt, 2, apply(Xt, 2, median)), 2, pmax(apply(Xt, 2, sd), 1e-12), "/")
    medoid <- setNames(Xt[which.min(rowSums(z^2)), ], calib)
    med_cur <- apply(S_cur, 2, median)
    mrun <- gp_ok(run_points(matrix(med_cur, nrow = 1, dimnames = list(NULL, calib)), "refmed"))
    ll_med <- if (nrow(mrun)) ll_true(mrun) else -Inf
    use_pool <- ll_med < ll_best - rf$pool_delta_ll
    point <- if (use_pool) medoid else med_cur
    saveRDS(list(S = S_cur, rhat = rhat_cur, point = point, use_pool = use_pool, ll_median_run = ll_med,
                 ll_best = ll_best, converged = conv, rounds = r, log = log_tab, val = vr_post), ck("refine"))
    write_csv(log_tab, file.path(C$out_dir, "refine_log.csv"))
    write.csv(round(vr_post, 3), file.path(C$out_dir, "gp_check_posterior_region_refined.csv"))
    write_csv(arrange(pl, desc(ll)), file.path(C$out_dir, "pool_runs_by_loglik.csv"))
    unlink(file.path(C$out_dir, c("evaluation_metrics.csv", "posterior_summary.csv", "n_response.csv",
                                  "sim_vs_obs_posterior_median.csv")))   # posterior may have changed
    status(sprintf("refine: %d round(s); point estimate = %s (median-run loglik %.1f vs best %.1f)", r,
                   if (use_pool) "DayCent pool medoid" else "refined posterior median", ll_med, ll_best))
  }
  if (file.exists(ck("refine"))) {
    rfs <- readRDS(ck("refine")); S <- rfs$S; med <- rfs$point; fin$rhat <- rfs$rhat
  }
  
  # ---- posterior evaluation vs CSU defaults ------------------------------------------
  if (has("post") && !file.exists(file.path(C$out_dir, "evaluation_metrics.csv"))) {
    near <- 0.05 * (hi - lo)
    post_tab <- tibble(label = calib, prior_lower = lo, prior_upper = hi, template = mu,
                       de_mode = de_best[calib],
                       median = apply(S, 2, median), point_estimate = med,
                       q025 = apply(S, 2, quantile, 0.025),
                       q975 = apply(S, 2, quantile, 0.975), rhat = fin$rhat[calib],
                       frac_near_lower = colMeans(sweep(S, 2, lo + near, "<")),
                       frac_near_upper = colMeans(sweep(S, 2, hi - near, ">")))
    write_csv(post_tab, file.path(C$out_dir, "posterior_summary.csv"))
    
    all_r <- unique(obs_all$RMNO)
    set.seed(g$seed + 4)
    draws <- rbind(median = med, S[sample.int(nrow(S), g$n_post), , drop = FALSE])
    sims <- map(seq_len(nrow(draws)), function(k) {
      rec <- run_cached(sprintf("post_%03d", k), make_params(draws[k, ]), all_r)
      if (!sim_ok(rec, length(all_r))) NULL else join_obs(rec$sim, obs_all)
    })
    if (is.null(sims[[1]])) stop(folder, ": DayCent failed at the posterior median")
    n_drop <- sum(vapply(sims[-1], is.null, logical(1))); sims <- Filter(Negate(is.null), sims)
    def <- join_obs(readRDS(ck("smoke"))$sim, obs_all)
    
    boot <- function(jl, B) {
      unit <- if (n_distinct(jl[[1]]$block) >= 5) "block" else "RMNO"
      mats <- map(jl, ~ split(.x, .x[[unit]])); units <- names(mats[[1]])
      M <- do.call(rbind, map(seq_len(B), function(b) {
        dd <- mats[[sample.int(length(mats), 1)]]; bs <- sample(units, replace = TRUE)
        all_metrics(bind_rows(map(seq_along(bs), ~ mutate(dd[[bs[.x]]], RMNO = paste0(RMNO, "_", .x)))))
      }))
      tibble(metric = colnames(M), ci_lo = apply(M, 2, quantile, 0.025, na.rm = TRUE),
             ci_hi = apply(M, 2, quantile, 0.975, na.rm = TRUE))
    }
    set.seed(g$seed + 5)
    met <- bind_rows(
      imap_dfr(subsets(sims[[1]]), function(jm, nm) {
        jl <- map(sims[-1], ~ subsets(.x)[[nm]])
        tibble(version = "calibrated", subset = nm, metric = names(all_metrics(jm)),
               value = all_metrics(jm), n_rm = n_distinct(jm$RMNO), n_rm_years = nrow(jm)) %>%
          left_join(boot(jl, g$n_boot), by = "metric") }),
      imap_dfr(subsets(def), function(jm, nm)
        tibble(version = "csu_default", subset = nm, metric = names(all_metrics(jm)),
               value = all_metrics(jm), n_rm = n_distinct(jm$RMNO), n_rm_years = nrow(jm))))
    write_csv(met, file.path(C$out_dir, "evaluation_metrics.csv"))
    write_csv(sims[[1]], file.path(C$out_dir, "sim_vs_obs_posterior_median.csv"))
    write_crop100(make_params(med), tpl, crop_id, C$int,
                  file.path(C$out_dir, "crop_calibrated_posterior_median.100"))
    status(sprintf("post: %d draws dropped; test_spatial sp_KGE cal %.3f vs default %.3f", n_drop,
                   met$value[met$version == "calibrated" & met$subset == "test_spatial" & met$metric == "sp_KGE"],
                   met$value[met$version == "csu_default" & met$subset == "test_spatial" & met$metric == "sp_KGE"]))
  }
  
  # ---- N-response curve: posterior median and CSU defaults -----------------------------
  if (has("mgmt") && !file.exists(file.path(C$out_dir, "n_response.csv"))) {
    all_r <- unique(obs_all$RMNO)
    nr <- map_dfr(c(posterior = TRUE, csu_default = FALSE), function(post) {
      v <- if (post) make_params(med) else make_params()
      map_dfr(g$mgmt_mult, function(m) {
        rec <- run_cached(sprintf("n_%s_%03d", if (post) "post" else "def", round(100 * m)), v, all_r, m)
        tibble(version = if (post) "posterior" else "csu_default", n_mult = m,
               mean_sim = mean(join_obs(rec$sim, obs_all)$sim))
      })
    }) %>% group_by(version) %>%
      mutate(pct_of_max = 100 * mean_sim / max(mean_sim)) %>% ungroup()
    write_csv(nr, file.path(C$out_dir, "n_response.csv"))
    status(sprintf("N response at 1x: posterior %.1f%%, default %.1f%% of max",
                   nr$pct_of_max[nr$version == "posterior" & nr$n_mult == 1],
                   nr$pct_of_max[nr$version == "csu_default" & nr$n_mult == 1]))
  }
  stop_cluster()
  invisible()
}

# =============================================================================
# 6. LOOP OVER CROPS
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a
RM_PTS <- sf::st_read(g$rm_coords_path, quiet = TRUE)
if (is.na(sf::st_crs(RM_PTS))) stop("RM centroid shapefile has no CRS (.prj missing?): ", g$rm_coords_path)
if (isTRUE(sf::st_is_longlat(RM_PTS))) RM_PTS <- sf::st_transform(RM_PTS, crs = sf::st_crs(26913))
cat("RM centroids:", nrow(RM_PTS), "features, CRS:", sf::st_crs(RM_PTS)$input, "\n")
C <- list(); W <- list()
for (k in which(crop_table$folder %in% crops_to_run)) {
  ct <- crop_table[k, ]
  res <- tryCatch(run_crop(ct$folder, ct$crop_id, ct$pulse),
                  error = function(e) { C$folder <<- ct$folder; status(paste("STOPPED:", conditionMessage(e))); NULL })
  stop_cluster()
}
cat("\nAll crops processed. Status log:", file.path(models_root, paste0("all_crops_status_", g$out_name, ".tsv")), "\n")