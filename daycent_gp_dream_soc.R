# =============================================================================
# DayCent SOC Decomposition Parameter Calibration
# Gaussian Process Emulator + DREAM(ZS)
# Prairie Soil Carbon Balance (PSCB) Network — Saskatchewan
# =============================================================================
#
# WORKFLOW:
#   Phase 1 — Latin Hypercube Sampling: run DayCent at ~300 parameter sets
#             across the 78 matched sites to build training data for GP
#   Phase 2 — GP Emulator: train a Gaussian Process on training data
#             mapping parameter vector -> network mean SOC per calibration year
#   Phase 3 — DREAM(ZS): run Bayesian MCMC using GP predictions instead of
#             DayCent directly — orders of magnitude faster per iteration
#   Phase 4 — Validation: spot-check GP posterior against real DayCent runs
#   Phase 5 — Output: posterior samples, MAP estimates, tightened DEoptim bounds
#
# WHY NETWORK MEAN AS TARGET:
#   Individual-site observed SOC carries large random error; the network mean
#   averages that noise out, so calibration targets the mean 0-20 cm SOC stock
#   across the 78 sites at each observed year (single global parameter set).
#
# CALIBRATION OBJECTIVE:  CHANGE-BASED
#   Target = observed CHANGE in network-mean soil SOC relative to the baseline
#   year, i.e. the trajectory / rate of change rather than the absolute stock.
#   Because site.100 is pegged to observed baseline SOC and the same quarter-
#   sections are resampled over time, the change is computed as a PAIRED site-
#   level difference (soc_site(y) - soc_site(baseline)) then network-averaged.
#   Pairing cancels the large between-site variance, so the change target has a
#   much tighter SE than the stock target — the calibration is thus forced to
#   match the trajectory. One Gaussian likelihood term per change year.
#   (To revert to a stock objective, swap make_gp_likelihood back to the
#    obs_targets/stock version — obs_targets is still computed for reference.)
#
# UNITS:
#   Observed:  Mg C/ha (from each site's observed.csv)
#   Simulated: soil SOM pools som1c(2)+som2c(2)+som3c = cols 6,8,9 at day 365,
#              in g C/m2, converted to Mg/ha via /100 (surface pools excluded).
#
# PARAMETERS CALIBRATED (5) — reduced set from Morris + Sobol screening on the
# RE-INITIALIZED 78-site network (soil-pool SOC target; site.100 pegged to
# observed 1996). At the corrected near-steady state the CO2-partitioning terms
# co-dominate with temperature, and the old decaying-model set (which had TEFF(2)
# and lacked the CO2 partition terms) no longer applies:
#   TEFF(1)   temperature-on-decomp coef 1 *      (dominant; T ~0.36-0.45)
#   PMCO2(2)  CO2 frac, soil metabolic litter     (co-dominant; sets CUE at steady state)
#   P1CO2A(1) CO2 frac, surface active pool       (gates surface->soil C flux **)
#   DEC5(2)   soil slow-pool (som2c(2)) rate, yr-1 (moderate; faded from decaying model)
#   P1CO2B(2) clay slope, soil active CO2 frac    (the one active-pool CO2 term with
#                                                  independent signal; P1CO2A(2)/PS1CO2(2)
#                                                  were redundant and are fixed)
#
#   * TEFF is build-specific: DDcentEVI splits the DEFAC into a 2-param
#     temperature function (TEFF) + 2-param moisture function (WEFF), NOT the
#     public manual's arctangent teff(1-4). Functional form per DDcentEVI source
#     (calcdefac) / Gurung et al. (2020); confirm before the methods section.
#   ** P1CO2A(1) is a SURFACE parameter that modulates soil C INPUTS via the
#      surface->soil flux; calibrating it may absorb a residue/NPP input imbalance
#      (cf. the gentle per-zone SOC drift) rather than a true decomposition signal.
#
# PARAMETERS FIXED (at their fix.100 values — only the 5 above are perturbed):
#   Screening showed the rest below the noise floor on the re-initialized model:
#     - passive pool (DEC4, P3CO2, PS1S3(2), PS2S3): frozen on this timescale
#     - fast litter rates (DEC1-3), moisture (WEFF1/2), TEFF(2)
#     - redundant active-pool CO2 terms P1CO2A(2), PS1CO2(2) (confounded w/ P1CO2B(2))
#   MAXCLTEF (tillage) — sites are no-till from baseline; unidentifiable.
#   Crop parameters — calibrated separately via DEoptim against yield data.
#
# =============================================================================

library(BayesianTools)   # DREAM(ZS) sampler
library(DiceKriging)     # Gaussian Process emulator (km function)
library(lhs)             # Latin Hypercube Sampling
library(dplyr)
library(readr)
library(parallel)
library(purrr)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE_DIR    <- "/home/preston/Projects/Daycent/SOC_Optimization/soc/qrts"
OUTPUT_DIR  <- "/home/preston/Projects/Daycent/SOC_Optimization/dream_output"

# Each site folder contains its own copy of the DayCent binary. Set the name of
# that binary here; it is invoked from inside each folder as ./<EXE_NAME>.
EXE_NAME    <- "DDcentEVI_rev491"          # <-- verify this matches your binary

# DayCent -s argument. Some DDcent builds want the schedule prefix WITHOUT the
# .sch extension; verify against a manual run before launching the full design.
SCHEDULE_ARG <- "Annual_Crop.sch"          # <-- verify

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# LHS training design (5 parameters — 300 is generous; ~150 is ample)
N_LHS           <- 300     # DayCent runs for GP training
N_LHS_VALIDATE  <- 30      # held-out runs for GP validation

# DREAM settings
N_ITERATIONS    <- 30000   # per chain; reduce to 5000 for pilot run
N_CHAINS        <- 3
N_BURNIN        <- 10000

# Calibration objective (change-based)
CALIB_YEARS   <- c(1996, 2018)        # 1996->2018 change only (2005 dropped: re-analysis uncertainty)
BASELINE_YEAR <- 1996                 # changes referenced to this (the init anchor)

# =============================================================================
# Decomposition parameter names and bounds  (5-parameter re-init set; bounds match
# the ranges explored in the Morris/Sobol screening). R names use dot notation
# matching read.table output; fix.100 names use DayCent parenthesis syntax.
#
#   R name       fix.100 name  Default    Prior range    Description
#   -----------  ------------  ---------  -------------  ---------------------------
#   TEFF.1.      TEFF(1)       20.0797    5.0–30.0       Temperature-on-decomp coef 1 *
#   PMCO2.2.     PMCO2(2)      0.669333   0.30–0.75      CO2 frac, soil metabolic litter decomp
#   P1CO2A.1.    P1CO2A(1)     0.6000     0.30–0.90      CO2 frac intercept, surface active pool
#   DEC5.2.      DEC5(2)       0.07925    0.07–0.25      Max decomp rate, soil SLOW OM (som2c(2)), yr-1
#   P1CO2B.2.    P1CO2B(2)     0.6800     0.34–1.02      Clay slope, soil active-pool CO2 frac
#
#   * DDcentEVI-specific DEFAC parameterization; NOT the public manual's
#     arctangent teff(1-4). Definitions per DDcentEVI source / Gurung et al. (2020).
# =============================================================================

PARAM_NAMES <- c("TEFF.1.", "PMCO2.2.", "P1CO2A.1.", "DEC5.2.", "P1CO2B.2.")

# fix.100 parameter names (parenthesis notation matching the actual file)
FIX100_PARAM_NAMES <- c("TEFF(1)", "PMCO2(2)", "P1CO2A(1)", "DEC5(2)", "P1CO2B(2)")

# Default values from fix.100 (verified against uploaded file)
PARAM_DEFAULTS <- c(
  "TEFF.1."   = 20.079654,
  "PMCO2.2."  = 0.669333,
  "P1CO2A.1." = 0.6000,
  "DEC5.2."   = 0.07925,
  "P1CO2B.2." = 0.6800
)

# Prior bounds — uniform (same ranges the screening explored)
LOWER_BOUNDS <- c(
  "TEFF.1." = 5.0,  "PMCO2.2." = 0.30, "P1CO2A.1." = 0.30, "DEC5.2." = 0.07, "P1CO2B.2." = 0.34
)

UPPER_BOUNDS <- c(
  "TEFF.1." = 30.0, "PMCO2.2." = 0.75, "P1CO2A.1." = 0.90, "DEC5.2." = 0.25, "P1CO2B.2." = 1.02
)

# Helper: consistent column name for a calibration year (e.g. 2005 -> "y2005")
yr_col <- function(y) paste0("y", y)

# Helper: cores available (respects SLURM allocation on Nibi)
n_cores_available <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (is.na(n) || n < 1) n <- max(1, detectCores() - 1)
  n
}

# =============================================================================
# 1. PREPARE OBSERVED NETWORK MEANS FROM PER-SITE observed.csv
# =============================================================================
# Each site folder has observed.csv with columns Year, soc_obs (Mg/ha).
# The calibration target is the network MEAN across sites at each observed year;
# se_obs (SE of the mean) captures the large site-level random error.

prepare_obs_network <- function(base_dir) {

  site_ids <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

  obs <- purrr::map_dfr(site_ids, function(s) {
    f <- file.path(base_dir, s, "observed.csv")
    if (!file.exists(f)) {
      warning(sprintf("No observed.csv in site folder: %s", s))
      return(NULL)
    }
    readr::read_csv(f, show_col_types = FALSE) %>%
      dplyr::transmute(site_id = s,
                       year    = as.integer(Year),
                       soc_obs = as.numeric(soc_obs))   # Mg/ha
  })

  if (nrow(obs) == 0) stop("No observed.csv files found under BASE_DIR")

  # restrict to the calibration-year group
  obs <- obs %>% dplyr::filter(!is.na(soc_obs), year %in% CALIB_YEARS)
  calib_years <- sort(unique(obs$year))
  if (!(BASELINE_YEAR %in% calib_years))
    stop(sprintf("BASELINE_YEAR %d not present in observed data", BASELINE_YEAR))

  # stock summary (reference only — not the objective; still trains the GPs)
  obs_targets <- obs %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(mean_obs = mean(soc_obs),
                     se_obs   = sd(soc_obs) / sqrt(dplyr::n()),
                     n        = dplyr::n(), .groups = "drop")

  # PAIRED change targets relative to the baseline year (the objective)
  base_soc <- obs %>% dplyr::filter(year == BASELINE_YEAR) %>%
    dplyr::select(site_id, base = soc_obs)
  change_targets <- do.call(rbind, lapply(setdiff(calib_years, BASELINE_YEAR), function(y) {
    d <- obs %>% dplyr::filter(year == y) %>% dplyr::select(site_id, soc = soc_obs) %>%
      dplyr::inner_join(base_soc, by = "site_id") %>%
      dplyr::mutate(delta = soc - base)
    data.frame(year        = y,
               mean_change = mean(d$delta),
               se_change   = sd(d$delta) / sqrt(nrow(d)),
               n_paired    = nrow(d))
  }))

  message(sprintf("Sites with observed.csv: %d | calibration years: %s | baseline: %d",
                  length(unique(obs$site_id)), paste(calib_years, collapse = ", "), BASELINE_YEAR))
  message("=== Observed stocks (Mg C/ha, reference) ===");  print(obs_targets)
  message("=== Observed PAIRED change vs baseline (the calibration target) ===")
  print(change_targets)
  message("  (compare se_change here to se_obs above — pairing is why the change is the tighter target)")

  list(
    obs_targets    = obs_targets,
    change_targets = change_targets,
    baseline_year  = BASELINE_YEAR,
    all_matched    = unique(obs$site_id),
    calib_years    = calib_years
  )
}

# =============================================================================
# 2. DAYCENT INTERFACE — write params, run, extract SOC
# =============================================================================

# Modifies the 10 decomposition parameters in the SITE'S OWN fix.100 (in place).
# Idempotent: only the leading numeric token of each matched line is replaced,
# so re-running with a new parameter set overwrites cleanly.
write_decomp_params <- function(params, site_dir,
                                fix100_names = FIX100_PARAM_NAMES,
                                r_names      = PARAM_NAMES) {
  fix_path <- file.path(site_dir, "fix.100")
  lines    <- readLines(fix_path)

  # fix.100 format: VALUE  PARAM_NAME  # comment  (value comes FIRST on the line)
  replace_param <- function(lines, fix100_name, value) {
    escaped <- gsub("([()])", "\\\\\\1", fix100_name)
    # match the parameter NAME (after the value) followed by ws / tab / comment / EOL,
    # only within the non-comment portion of the line
    idx <- grep(paste0("^[^#]*", escaped, "(\\s|$|\t|#)"), lines, perl = TRUE)
    if (length(idx) == 0) {
      warning(sprintf("Parameter '%s' not found in fix.100", fix100_name)); return(lines)
    }
    if (length(idx) > 1) {
      warning(sprintf("Parameter '%s' matched %d lines — using first", fix100_name, length(idx)))
    }
    # replace the leading numeric token with the new value.
    # format="f" avoids scientific notation that DayCent's Fortran reader may reject.
    lines[idx[1]] <- sub("^(\\s*)([-0-9.eE+]+)",
                         sprintf("\\1%s", formatC(value, format = "f", digits = 7)),
                         lines[idx[1]])
    lines
  }

  name_map <- setNames(fix100_names, r_names)
  for (rname in names(params)) {
    fix_name <- name_map[[rname]]
    if (is.na(fix_name) || is.null(fix_name)) {
      warning(sprintf("No fix.100 mapping for parameter '%s'", rname)); next
    }
    lines <- replace_param(lines, fix_name, params[[rname]])
  }
  writeLines(lines, fix_path)
}

run_daycent_site <- function(site_id, params, base_dir, exe_name, calib_years) {

  site_dir <- file.path(base_dir, site_id)
  write_decomp_params(params, site_dir)

  cmd <- sprintf("cd '%s' && rm -f output.bin && './%s' -s %s -n output > /dev/null 2>&1",
                 site_dir, exe_name, SCHEDULE_ARG)
  if (system(cmd) != 0) { warning(sprintf("DayCent failed: %s", site_id)); return(NULL) }

  soilc_path <- file.path(site_dir, "soilc.out")
  if (!file.exists(soilc_path)) { warning(sprintf("No soilc.out: %s", site_id)); return(NULL) }

  soilc <- tryCatch(read.table(soilc_path, header = TRUE),
                    error = function(e) { warning(e$message); NULL })
  if (is.null(soilc)) return(NULL)
  if (ncol(soilc) < 9) { warning(sprintf("soilc.out has <9 columns: %s", site_id)); return(NULL) }

  # Total SOC = soil SOM pools only: som1c(2)+som2c(2)+som3c = cols 6,8,9 (g C/m2),
  # then g C/m2 -> Mg C/ha (/100). Surface pools (cols 5,7) and litter (3,4) excluded.
  # Positional indexing avoids read.table's make.names mangling of the som*(*) headers.
  data.frame(
    year         = as.integer(soilc$time),
    dayofyr      = soilc$dayofyr,
    soc_020_Mgha = rowSums(soilc[, c(6, 8, 9)]) / 100
  ) %>%
    filter(dayofyr == 365, year %in% calib_years) %>%
    select(year, soc_020_Mgha) %>%
    mutate(site_id = site_id)
}

# Run all matched sites and return network means at calibration years (Mg/ha)
run_network_mean <- function(params, matched_sites, base_dir, exe_name, calib_years) {

  results <- mclapply(
    matched_sites,
    function(s) run_daycent_site(s, params, base_dir, exe_name, calib_years),
    mc.cores = n_cores_available()
  )
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)

  bind_rows(results) %>%
    group_by(year) %>%
    summarise(mean_sim = mean(soc_020_Mgha, na.rm = TRUE),
              n_sim    = n(), .groups = "drop")
}

# =============================================================================
# 3. PHASE 1 — LATIN HYPERCUBE SAMPLING (training data for GP)
# =============================================================================

run_lhs_training <- function(n_lhs, n_validate, param_names, lower, upper,
                             matched_sites, base_dir, exe_name, calib_years,
                             output_dir) {

  message(sprintf("\n=== Phase 1: LHS Training Design (%d + %d runs) ===", n_lhs, n_validate))

  n_total <- n_lhs + n_validate
  set.seed(42)
  lhs_unit <- randomLHS(n_total, length(param_names))    # unit hypercube

  # Scale unit hypercube to parameter bounds
  lhs_params <- matrix(0, nrow = n_total, ncol = length(param_names))
  for (j in seq_along(param_names)) {
    lhs_params[, j] <- lower[j] + lhs_unit[, j] * (upper[j] - lower[j])
  }
  colnames(lhs_params) <- param_names

  n_cores <- n_cores_available()
  message(sprintf("Running %d DayCent network evaluations for GP training...", n_total))
  message(sprintf("Rough wall time: ~%.1f h  (%d evals x %d sites x ~18s / %d cores)",
                  n_total * length(matched_sites) * 18 / 3600 / n_cores,
                  n_total, length(matched_sites), n_cores))

  year_cols <- yr_col(calib_years)
  training_results <- vector("list", n_total)

  for (i in seq_len(n_total)) {
    if (i %% 10 == 0) message(sprintf("  LHS run %d / %d", i, n_total))

    params_i  <- lhs_params[i, ]
    sim_means <- tryCatch(
      run_network_mean(params_i, matched_sites, base_dir, exe_name, calib_years),
      error = function(e) NULL
    )
    if (is.null(sim_means)) next

    # reindex to all calibration years (NA where a year is missing)
    sim_vec <- setNames(rep(NA_real_, length(calib_years)), year_cols)
    sim_vec[yr_col(sim_means$year)] <- sim_means$mean_sim

    training_results[[i]] <- cbind(
      data.frame(run_id = i),
      as.data.frame(t(params_i)),
      as.data.frame(t(sim_vec))
    )
  }

  training_df <- bind_rows(training_results)
  # keep only rows with a simulated mean for every calibration year
  complete <- stats::complete.cases(training_df[, year_cols, drop = FALSE])
  training_df <- training_df[complete, , drop = FALSE]

  message(sprintf("Successful LHS runs: %d / %d", nrow(training_df), n_total))

  training_df$split <- ifelse(training_df$run_id <= n_lhs, "train", "validate")

  write_csv(training_df, file.path(output_dir, "lhs_training_data.csv"))
  message(sprintf("LHS data saved to %s/lhs_training_data.csv", output_dir))

  training_df
}

# =============================================================================
# 4. PHASE 2 — GAUSSIAN PROCESS EMULATOR (one GP per calibration year)
# =============================================================================

train_gp_emulator <- function(training_df, param_names, calib_years, output_dir) {

  message("\n=== Phase 2: Training GP Emulator ===")

  train <- training_df %>% filter(split == "train")
  valid <- training_df %>% filter(split == "validate")

  X_train <- as.data.frame(train[, param_names, drop = FALSE])
  X_valid <- as.data.frame(valid[, param_names, drop = FALSE])

  gp_models <- list()
  val_rows  <- list()

  png(file.path(output_dir, "gp_validation.png"),
      width = 500 * length(calib_years), height = 500)
  par(mfrow = c(1, length(calib_years)))

  for (y in calib_years) {
    col <- yr_col(y)
    message(sprintf("  Fitting GP for %d network mean SOC...", y))

    # Matern 5/2 kernel — smooth, robust for physical-model outputs
    gp <- km(
      formula      = ~1,
      design       = X_train,
      response     = train[[col]],
      covtype      = "matern5_2",
      optim.method = "BFGS",
      control      = list(trace = FALSE)
    )

    pred <- predict(gp, newdata = X_valid, type = "UK")
    rmse <- sqrt(mean((pred$mean - valid[[col]])^2))
    r2   <- cor(pred$mean, valid[[col]])^2
    message(sprintf("  %d GP: RMSE = %.3f Mg/ha | R2 = %.3f", y, rmse, r2))

    plot(valid[[col]], pred$mean,
         xlab = "DayCent (Mg/ha)", ylab = "GP prediction (Mg/ha)",
         main = sprintf("GP %d | RMSE=%.2f | R\u00b2=%.3f", y, rmse, r2),
         pch = 16, col = "steelblue")
    abline(0, 1, col = "red", lty = 2)

    if (r2 < 0.90) {
      warning(sprintf("GP R\u00b2 < 0.90 for %d — increase N_LHS or check for failed DayCent runs", y))
    }

    gp_models[[col]] <- gp
    val_rows[[col]]  <- data.frame(year = y, rmse_Mgha = rmse, r_squared = r2)
  }
  dev.off()

  validation_stats <- bind_rows(val_rows)
  write_csv(validation_stats, file.path(output_dir, "gp_validation_stats.csv"))
  message("\n=== GP validation ===")
  print(validation_stats)

  gp_models   # named list keyed by "y<year>"
}

# =============================================================================
# 5. PHASE 3 — DREAM(ZS) USING GP EMULATOR
# =============================================================================

make_gp_likelihood <- function(gp_models, change_targets, baseline_year, param_names) {
  # Returns a closure: params -> log likelihood (CHANGE-based objective).
  # For each change year y: Δsim(y) = gp_y(params) - gp_baseline(params), compared
  # to the observed PAIRED change. sigma combines the paired observation SE with
  # the GP predictive SD of both years (variances added in quadrature). One term
  # per change year; the baseline year is the anchor (sim(baseline) ~ observed by
  # construction of the re-initialization), so this targets the trajectory shape.
  base_col <- yr_col(baseline_year)

  function(params) {
    param_df <- as.data.frame(t(params))
    colnames(param_df) <- param_names

    gp_b <- gp_models[[base_col]]
    if (is.null(gp_b)) return(-Inf)
    pb <- tryCatch(predict(gp_b, newdata = param_df, type = "UK"), error = function(e) NULL)
    if (is.null(pb)) return(-Inf)

    ll <- 0
    for (i in seq_len(nrow(change_targets))) {
      gp <- gp_models[[yr_col(change_targets$year[i])]]
      if (is.null(gp)) return(-Inf)
      py <- tryCatch(predict(gp, newdata = param_df, type = "UK"), error = function(e) NULL)
      if (is.null(py)) return(-Inf)

      dsim  <- py$mean - pb$mean
      sigma <- sqrt(change_targets$se_change[i]^2 + py$sd^2 + pb$sd^2)
      ll    <- ll + dnorm(change_targets$mean_change[i], mean = dsim, sd = sigma, log = TRUE)
    }
    if (!is.finite(ll)) return(-Inf)
    ll
  }
}

run_dream <- function(gp_likelihood, param_names, lower, upper,
                      n_iterations, n_chains, n_burnin, run_parallel, output_dir) {

  message("\n=== Phase 3: DREAM(ZS) via GP Emulator ===")
  message(sprintf("  %d iterations x %d chains", n_iterations, n_chains))

  prior <- createUniformPrior(lower = lower, upper = upper)

  bayesian_setup <- createBayesianSetup(
    likelihood = gp_likelihood,
    prior      = prior,
    names      = param_names,
    parallel   = run_parallel
  )

  runMCMC(
    bayesianSetup = bayesian_setup,
    sampler       = "DEzs",
    settings      = list(
      iterations = n_iterations,
      nrChains   = n_chains,
      burnin     = n_burnin,
      message    = TRUE
    )
  )
}

# =============================================================================
# 6. PHASE 4 — VALIDATION: check GP posteriors against real DayCent
# =============================================================================

validate_posterior <- function(dream_result, param_names, matched_sites,
                               base_dir, exe_name, calib_years, obs_targets,
                               change_targets, baseline_year,
                               n_validate = 20, output_dir) {

  message(sprintf("\n=== Phase 4: Posterior Validation (%d DayCent runs) ===", n_validate))

  # burnin already applied in runMCMC; do not trim again
  post_samples <- getSample(dream_result)
  colnames(post_samples) <- param_names

  set.seed(99)
  val_idx    <- sample(nrow(post_samples), n_validate)
  val_params <- post_samples[val_idx, , drop = FALSE]

  year_cols   <- yr_col(calib_years)
  val_results <- vector("list", n_validate)

  for (i in seq_len(n_validate)) {
    if (i %% 5 == 0) message(sprintf("  Validation run %d / %d", i, n_validate))
    params_i <- val_params[i, ]

    sim_means <- tryCatch(
      run_network_mean(params_i, matched_sites, base_dir, exe_name, calib_years),
      error = function(e) NULL
    )
    if (is.null(sim_means)) next

    sim_vec <- setNames(rep(NA_real_, length(calib_years)), paste0("daycent_", calib_years))
    sim_vec[paste0("daycent_", sim_means$year)] <- sim_means$mean_sim

    val_results[[i]] <- cbind(
      data.frame(val_run = i),
      as.data.frame(t(params_i)),
      as.data.frame(t(sim_vec))
    )
  }

  val_df <- bind_rows(val_results)

  message("\n  === Posterior predictive: STOCKS (Mg/ha, reference) ===")
  for (y in calib_years) {
    obs_y <- obs_targets$mean_obs[obs_targets$year == y]
    dc    <- val_df[[paste0("daycent_", y)]]
    message(sprintf("  %d  observed: %.2f | DayCent posterior mean: %.2f (sd=%.2f)",
                    y, obs_y, mean(dc, na.rm = TRUE), sd(dc, na.rm = TRUE)))
  }

  message("\n  === Posterior predictive: CHANGE from baseline (Mg/ha, the objective) ===")
  base_dc <- val_df[[paste0("daycent_", baseline_year)]]
  for (i in seq_len(nrow(change_targets))) {
    y    <- change_targets$year[i]
    dsim <- val_df[[paste0("daycent_", y)]] - base_dc
    message(sprintf("  %d->%d  observed Δ: %+.2f (se %.2f) | DayCent posterior Δ: %+.2f (sd %.2f)",
                    baseline_year, y, change_targets$mean_change[i], change_targets$se_change[i],
                    mean(dsim, na.rm = TRUE), sd(dsim, na.rm = TRUE)))
  }

  write_csv(val_df, file.path(output_dir, "posterior_validation.csv"))
  val_df
}

# =============================================================================
# 7. DIAGNOSTICS AND OUTPUT
# =============================================================================

save_results <- function(dream_result, param_names, output_dir) {

  message("\n=== Convergence Diagnostics ===")
  gelman <- gelmanDiagnostics(dream_result)
  print(gelman)

  if (any(gelman$psrf[, 1] > 1.1, na.rm = TRUE)) {
    warning("R-hat > 1.1 for some parameters — consider more iterations")
  } else {
    message("Convergence good (all R-hat < 1.1)")
  }

  # burnin already applied in runMCMC
  post_samples <- getSample(dream_result)
  colnames(post_samples) <- param_names

  map_params <- MAP(dream_result)$parametersMAP

  summary_df <- data.frame(
    parameter = param_names,
    MAP       = map_params,
    post_mean = colMeans(post_samples),
    post_sd   = apply(post_samples, 2, sd),
    lower_95  = apply(post_samples, 2, quantile, 0.025),
    upper_95  = apply(post_samples, 2, quantile, 0.975)
  )
  message("\n=== Posterior Summary ===")
  print(summary_df)

  saveRDS(dream_result, file.path(output_dir, "dream_result.rds"))
  write_csv(as.data.frame(post_samples), file.path(output_dir, "posterior_samples.csv"))
  write_csv(summary_df, file.path(output_dir, "map_estimates.csv"))

  png(file.path(output_dir, "trace_plots.png"), width = 1200, height = 800)
  plot(dream_result); dev.off()

  png(file.path(output_dir, "marginal_posteriors.png"), width = 1200, height = 800)
  marginalPlot(dream_result); dev.off()

  png(file.path(output_dir, "correlation_posteriors.png"), width = 1000, height = 1000)
  correlationPlot(dream_result); dev.off()

  message(sprintf("\nAll outputs saved to: %s", output_dir))

  list(dream_result = dream_result, map_params = map_params,
       post_samples = post_samples, summary_df = summary_df)
}

# =============================================================================
# 8. HELPER: EXTRACT DEOPTIM BOUNDS FROM DREAM POSTERIORS
# =============================================================================

extract_deoptim_bounds <- function(posterior_samples_csv, ci_width = 0.95,
                                   output_path = NULL) {

  samples <- read_csv(posterior_samples_csv, show_col_types = FALSE)
  alpha   <- (1 - ci_width) / 2

  bounds <- data.frame(
    parameter = names(samples),
    lower     = apply(samples, 2, quantile, alpha),
    upper     = apply(samples, 2, quantile, 1 - alpha),
    MAP       = apply(samples, 2, function(x) { d <- density(x); d$x[which.max(d$y)] }),
    post_mean = colMeans(samples),
    post_sd   = apply(samples, 2, sd)
  )

  message(sprintf("\n=== DEoptim bounds from DREAM %d%% CI ===", ci_width * 100))
  print(bounds)

  if (!is.null(output_path)) {
    write_csv(bounds, output_path)
    message(sprintf("Bounds saved to: %s", output_path))
  }
  invisible(bounds)
}

# =============================================================================
# 9. MAIN WORKFLOW
# =============================================================================

main <- function() {

  message("=== DayCent SOC GP + DREAM Calibration ===")
  message(Sys.time())

  # --- 9a. Observed network means ---
  message("\n[1/5] Preparing observed network means from per-site observed.csv...")
  obs_data    <- prepare_obs_network(BASE_DIR)
  calib_years <- obs_data$calib_years

  # --- 9b. LHS training runs (skip if cached) ---
  lhs_path <- file.path(OUTPUT_DIR, "lhs_training_data.csv")
  if (file.exists(lhs_path)) {
    message(sprintf("\n[2/5] Loading existing LHS training data from %s", lhs_path))
    training_df <- read_csv(lhs_path, show_col_types = FALSE)
  } else {
    message("\n[2/5] Running LHS training design (Phase 1)...")
    training_df <- run_lhs_training(
      n_lhs         = N_LHS,
      n_validate    = N_LHS_VALIDATE,
      param_names   = PARAM_NAMES,
      lower         = LOWER_BOUNDS,
      upper         = UPPER_BOUNDS,
      matched_sites = obs_data$all_matched,
      base_dir      = BASE_DIR,
      exe_name      = EXE_NAME,
      calib_years   = calib_years,
      output_dir    = OUTPUT_DIR
    )
  }

  # --- 9c. Train GP emulator ---
  message("\n[3/5] Training GP emulator (Phase 2)...")
  gp_models <- train_gp_emulator(training_df, PARAM_NAMES, calib_years, OUTPUT_DIR)
  saveRDS(gp_models, file.path(OUTPUT_DIR, "gp_models.rds"))

  # --- 9d. DREAM via GP ---
  message("\n[4/5] Running DREAM(ZS) via GP emulator (Phase 3)...")
  gp_likelihood <- make_gp_likelihood(gp_models, obs_data$change_targets,
                                      obs_data$baseline_year, PARAM_NAMES)

  # Sanity check at fix.100 default values before the full run
  test_ll <- gp_likelihood(PARAM_DEFAULTS)
  message(sprintf("  Log likelihood at fix.100 defaults: %.3f", test_ll))
  if (!is.finite(test_ll)) stop("GP likelihood returned -Inf at defaults — check GP training")

  dream_result <- run_dream(
    gp_likelihood = gp_likelihood,
    param_names   = PARAM_NAMES,
    lower         = LOWER_BOUNDS,
    upper         = UPPER_BOUNDS,
    n_iterations  = N_ITERATIONS,
    n_chains      = N_CHAINS,
    n_burnin      = N_BURNIN,
    run_parallel  = FALSE,
    output_dir    = OUTPUT_DIR
  )

  # --- 9e. Save results and diagnostics ---
  message("\n[5/5] Saving results and diagnostics...")
  results <- save_results(dream_result, PARAM_NAMES, OUTPUT_DIR)

  # --- 9f. Posterior validation against real DayCent ---
  val_df <- validate_posterior(
    dream_result  = dream_result,
    param_names   = PARAM_NAMES,
    matched_sites = obs_data$all_matched,
    base_dir      = BASE_DIR,
    exe_name      = EXE_NAME,
    calib_years   = calib_years,
    obs_targets   = obs_data$obs_targets,
    change_targets = obs_data$change_targets,
    baseline_year  = obs_data$baseline_year,
    n_validate    = 20,
    output_dir    = OUTPUT_DIR
  )

  # --- 9g. Extract tightened DEoptim bounds ---
  extract_deoptim_bounds(
    posterior_samples_csv = file.path(OUTPUT_DIR, "posterior_samples.csv"),
    output_path           = file.path(OUTPUT_DIR, "deoptim_bounds.csv")
  )

  message("\n=== Complete ===")
  message(Sys.time())

  invisible(c(results, list(gp_models = gp_models, val_df = val_df,
                            obs_data = obs_data)))
}

# =============================================================================
# RUN
# =============================================================================

result <- main()

# To reload a completed run and extract bounds:
# dream_result <- readRDS(file.path(OUTPUT_DIR, "dream_result.rds"))
# gp_models    <- readRDS(file.path(OUTPUT_DIR, "gp_models.rds"))
# bounds       <- extract_deoptim_bounds(
#   file.path(OUTPUT_DIR, "posterior_samples.csv"),
#   output_path = file.path(OUTPUT_DIR, "deoptim_bounds.csv")
# )
