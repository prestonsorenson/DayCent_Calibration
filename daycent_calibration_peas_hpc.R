# =============================================================================
# DayCent PEA Field Peas Calibration — RM Level
# Saskatchewan — DDcentEVI_rev491
# =============================================================================
# Calibration targets are Rural Municipality (RM) means across quarter sections
# with >= 10 years of SCIC yield records, averaged over 1999-2019.
#
# Objective function mirrors the quarter-section calibration:
#   - RMSE on RM-mean yields
#   - Bias (normalised)
#   - Spatial correlation across RMs
#   - Anomaly correlation (RM-level year-to-year, site-normalised)
#   - Range exceedance penalty
#
# Anomaly weight retained: the annual RM series (1999-2019) carries temporal
# signal even though the primary calibration target is the RM mean. Anomaly
# correlation constrains the temporal response without overfitting to noisy
# single-year quarter-section records.
#
# OPTIMIZED (6 parameters) — same set as quarter-section calibration:
#   PRDX(1), PPDF(1), PPDF(2), WSCOEF(1), WSCOEF(2), CLAYPG
#
# FIXED at literature values (Saskatchewan HRSW):
#   HIMAX, HIWSF, DDBASE, BASETEMP, EMAX, TMPGERM,
#   FRTC(1), FRTC(3), PLTMRF, HIMON(1), CFRTCW(1)
# =============================================================================

library(tidyverse)
library(DEoptim)
library(parallel)
library(fs)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

config <- list(

  daycent_exe    = "/home/prs352/Daycent/Linux_Version_491/DDcentEVI_rev491",
  template_dir   = "/home/prs352/Daycent/Crop_Optimization/Models/Peas/Copy_Folder",
  sites_dir      = "/home/prs352/Daycent/Crop_Optimization/Models/Peas/qrts",
  results_dir    = "/home/prs352/Daycent/Crop_Optimization/Models/Peas/results",
  obs_data_path  = "/home/prs352/Daycent/Crop_Optimization/Models/input_data/peas_c_gmsq.csv",
  crop100_name   = "crop.100",
  schedule_name  = "Crop_Opt.sch",
  output_prefix  = "harvest",
  crop_id        = "PEA",

  n_cores        = parallel::detectCores() - 1,
  deoptim_NP     = 60,     # 10x n_params
  deoptim_iter   = 126,
  deoptim_CR     = 0.9,
  deoptim_F      = 0.8,

  # Objective weights
  # Anomaly weight reduced relative to spatial — RM-mean spatial gradient
  # is the primary calibration target; anomaly constrains temporal response.
  weights = list(
    bias    = 0.3,
    spatial = 1.0,
    anomaly = 0.3
  ),

  # Minimum matched RM-years to accept an objective evaluation
  min_rm_years = 1007
)

dir_create(config$results_dir)

# Clear stale calibration outputs
for (.f in c("calibration_log.csv", "best_parameters.csv")) {
  .fp <- file.path(config$results_dir, .f)
  if (file.exists(.fp)) file.remove(.fp)
}

# =============================================================================
# 2. FIXED PARAMETERS
# =============================================================================

fixed_params <- list(
  HIMAX    = 0.44,
  HIWSF    = 0.7,
  DDBASE   = 1100.0,
  BASETEMP = 0.0,
  EMAX     = 0.80,
  TMPGERM  = 2.0,
  FRTC1    = 0.45,
  FRTC3    = 100.0,
  PLTMRF   = 0.50,
  HIMON1   = 2,
  CFRTCW1  = 0.40
)

# =============================================================================
# 3. LOAD OBSERVED DATA
# =============================================================================
# Read observed.csv from each RM folder — annual RM-mean yields (cgrain, g C m-2)
# Only RMs with folders present are included.

rm_folders <- list.dirs(config$sites_dir, full.names = FALSE, recursive = FALSE)
rm_folders <- rm_folders[nchar(rm_folders) > 0]
cat(sprintf("Found %d RM folders\n", length(rm_folders)))

# Load single observed yield CSV — Year, RMNO, cgrain
obs_all <- read_csv(config$obs_data_path, show_col_types = FALSE) %>%
  mutate(
    RMNO   = as.character(RMNO),
    Year   = as.integer(Year),
    grain_c = as.numeric(grain_c)
  ) %>%
  filter(RMNO %in% rm_folders)   # restrict to RMs that have a model folder

cat(sprintf("Loaded %d RM-year observations across %d RMs\n",
            nrow(obs_all), n_distinct(obs_all$RMNO)))

# RM-mean (primary calibration target)
obs_rm_mean <- obs_all %>%
  group_by(RMNO) %>%
  summarise(obs_mean_grain_c = mean(grain_c, na.rm = TRUE), .groups = "drop")

# =============================================================================
# 4. WRITE crop.100 WITH CALIBRATION PARAMETERS
# =============================================================================

write_crop100 <- function(params, template_path, output_path, crop_name = "PEA") {

  lines <- readLines(template_path)

  crop_start  <- grep(paste0("^", crop_name, "\\b"), lines)
  if (length(crop_start) == 0) stop(paste("Crop", crop_name, "not found"))

  crop_headers <- grep("^[A-Z][A-Z0-9]+\\s+\\S", lines)
  next_header  <- crop_headers[crop_headers > crop_start][1]
  crop_end     <- ifelse(is.na(next_header), length(lines), next_header - 1)

  replace_param <- function(lines, block_start, block_end, param_pattern, new_value) {
    pattern <- paste0("^([0-9.-]+)\\s+(", param_pattern, ")(\\s.*)?$")
    idx     <- grep(pattern, lines[block_start:block_end], perl = TRUE)
    if (length(idx) == 0) {
      warning(paste("Parameter", param_pattern, "not found in", crop_name, "block"))
      return(lines)
    }
    abs_idx        <- block_start + idx - 1
    old_line       <- lines[abs_idx]
    old_value      <- regmatches(old_line, regexpr("^[0-9.-]+", old_line))
    new_value_str  <- formatC(new_value, format = "f", digits = 4)
    lines[abs_idx] <- sub(old_value, new_value_str, old_line, fixed = TRUE)
    return(lines)
  }

  if (!is.null(params$PRDX1))   lines <- replace_param(lines, crop_start, crop_end, "PRDX\\(1\\)",    params$PRDX1)
  if (!is.null(params$PPDF1))   lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(1\\)",    params$PPDF1)
  if (!is.null(params$PPDF2))   lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(2\\)",    params$PPDF2)
  if (!is.null(params$WSCOEF1)) lines <- replace_param(lines, crop_start, crop_end, "WSCOEF\\(1\\)",  params$WSCOEF1)
  if (!is.null(params$WSCOEF2)) lines <- replace_param(lines, crop_start, crop_end, "WSCOEF\\(2\\)",  params$WSCOEF2)
  if (!is.null(params$CLAYPG))  lines <- replace_param(lines, crop_start, crop_end, "CLAYPG",         round(params$CLAYPG))
  if (!is.null(params$HIMAX))   lines <- replace_param(lines, crop_start, crop_end, "HIMAX",          params$HIMAX)
  if (!is.null(params$HIWSF))   lines <- replace_param(lines, crop_start, crop_end, "HIWSF",          params$HIWSF)
  if (!is.null(params$DDBASE))  lines <- replace_param(lines, crop_start, crop_end, "DDBASE",         params$DDBASE)
  if (!is.null(params$BASETEMP))lines <- replace_param(lines, crop_start, crop_end, "BASETEMP(?!\\()",params$BASETEMP)
  if (!is.null(params$EMAX))    lines <- replace_param(lines, crop_start, crop_end, "EMAX",           params$EMAX)
  if (!is.null(params$TMPGERM)) lines <- replace_param(lines, crop_start, crop_end, "TMPGERM",        params$TMPGERM)
  if (!is.null(params$FRTC1))   lines <- replace_param(lines, crop_start, crop_end, "FRTC\\(1\\)",    params$FRTC1)
  if (!is.null(params$FRTC3))   lines <- replace_param(lines, crop_start, crop_end, "FRTC\\(3\\)",    params$FRTC3)
  if (!is.null(params$PLTMRF))  lines <- replace_param(lines, crop_start, crop_end, "PLTMRF",         params$PLTMRF)
  if (!is.null(params$HIMON1))  lines <- replace_param(lines, crop_start, crop_end, "HIMON\\(1\\)",   round(params$HIMON1))
  if (!is.null(params$CFRTCW1)) lines <- replace_param(lines, crop_start, crop_end, "CFRTCW\\(1\\)",  params$CFRTCW1)

  writeLines(lines, output_path)
  invisible(NULL)
}

# =============================================================================
# 5. RUN DAYCENT FOR ONE RM
# =============================================================================

run_daycent_rm <- function(rm_id, params, config) {

  site_dir     <- file.path(config$sites_dir, rm_id)
  crop100_path <- file.path(site_dir, config$crop100_name)
  output_csv   <- file.path(site_dir, "harvest.csv")

  # Write parameters to this RM's crop.100
  # Use template_dir as the source template so the base values are always clean
  template_path <- file.path(config$template_dir, config$crop100_name)
  tryCatch(
    write_crop100(params, template_path, crop100_path, config$crop_id),
    error = function(e) {
      warning(paste("write_crop100 failed for", rm_id, ":", e$message))
      return(NULL)
    }
  )

  # Remove stale outputs
  stale <- list.files(site_dir,
                      pattern    = paste0("^", config$output_prefix, ".*\\.(bin|lis|csv|out)$"),
                      full.names = TRUE)
  if (length(stale) > 0) file.remove(stale)

  # Run DayCent using absolute exe path
  cmd <- paste0(
    "cd ", shQuote(site_dir),
    " && ", shQuote(config$daycent_exe),
    " -s ", config$schedule_name,
    " -n ", config$output_prefix,
    " > daycent.log 2>&1"
  )

  exit_code <- system(cmd)

  if (exit_code != 0 || !file.exists(output_csv)) {
    warning(paste("DayCent failed:", rm_id, "— check daycent.log"))
    return(NULL)
  }

  # Remove .bin to prevent stale file errors on next evaluation
  bin_file <- file.path(site_dir, paste0(config$output_prefix, ".bin"))
  if (file.exists(bin_file)) file.remove(bin_file)

  return(output_csv)
}

# =============================================================================
# 6. EXTRACT SIMULATED YIELD FROM harvest.csv
# =============================================================================

extract_yield_rm <- function(output_csv, rm_id) {

  if (is.null(output_csv) || !file.exists(output_csv)) return(NULL)

  tryCatch({
    read_csv(output_csv, show_col_types = FALSE) %>%
      mutate(
        RMNO    = rm_id,
        Year    = as.integer(floor(as.numeric(time))),
        grain_c = as.numeric(cgrain)
      ) %>%
      select(RMNO, Year, grain_c)
  }, error = function(e) {
    warning(paste("Could not read harvest.csv for:", rm_id, "-", e$message))
    return(NULL)
  })
}

# =============================================================================
# 7. RUN ALL RMs IN PARALLEL
# =============================================================================

run_all_rms <- function(params, config) {

  clusterExport(cl, varlist = "params", envir = environment())

  results <- parLapply(
    cl,
    rm_folders,
    function(rm_id) {
      out_csv <- run_daycent_rm(rm_id, params, config)
      extract_yield_rm(out_csv, rm_id)
    }
  )

  bind_rows(Filter(Negate(is.null), results))
}

# =============================================================================
# 8. OBJECTIVE FUNCTION
# =============================================================================

objective_fn <- function(param_vec) {

  params <- c(
    list(
      PRDX1   = param_vec[1],
      PPDF1   = param_vec[2],
      PPDF2   = param_vec[3],
      WSCOEF1 = param_vec[4],
      WSCOEF2 = param_vec[5],
      CLAYPG  = round(param_vec[6])
    ),
    fixed_params
  )

  simulated <- tryCatch(
    run_all_rms(params, config),
    error = function(e) {
      warning(paste("run_all_rms failed:", e$message))
      return(NULL)
    }
  )

  if (is.null(simulated) || nrow(simulated) == 0) return(Inf)

  # Join simulated annual series to observed annual series
  joined <- simulated %>%
    inner_join(obs_all %>% rename(obs_grain_c = grain_c),
               by = c("RMNO", "Year"))

  if (nrow(joined) < config$min_rm_years) return(Inf)

  # -------------------------------------------------------------------
  # 1. RM-MEAN SPATIAL PERFORMANCE
  # -------------------------------------------------------------------

  rm_means <- joined %>%
    group_by(RMNO) %>%
    summarise(
      sim_mean = mean(grain_c,     na.rm = TRUE),
      obs_mean = mean(obs_grain_c, na.rm = TRUE),
      .groups  = "drop"
    )

  rmse <- sqrt(mean((rm_means$sim_mean - rm_means$obs_mean)^2))

  obs_grand_mean <- mean(rm_means$obs_mean)
  bias_norm      <- mean(rm_means$sim_mean - rm_means$obs_mean) / obs_grand_mean

  spatial_cor <- tryCatch(
    cor(rm_means$sim_mean, rm_means$obs_mean, method = "pearson"),
    error = function(e) 0
  )

  # -------------------------------------------------------------------
  # 2. ANOMALY CORRELATION (temporal, RM-normalised)
  # Site-mean and SD removed so every RM contributes equally regardless
  # of its position on the spatial gradient.
  # -------------------------------------------------------------------

  anomalies <- joined %>%
    group_by(RMNO) %>%
    filter(n() >= 4) %>%
    mutate(
      obs_anom = (obs_grain_c - mean(obs_grain_c)) / (sd(obs_grain_c) + 1e-6),
      sim_anom = (grain_c     - mean(grain_c))     / (sd(grain_c)     + 1e-6)
    ) %>%
    ungroup()

  anomaly_cor <- if (nrow(anomalies) >= 50) {
    tryCatch(
      cor(anomalies$sim_anom, anomalies$obs_anom, method = "pearson"),
      error = function(e) 0
    )
  } else 0

  # -------------------------------------------------------------------
  # 3. RANGE EXCEEDANCE PENALTY
  # -------------------------------------------------------------------

  obs_min   <- min(joined$obs_grain_c)
  obs_max   <- max(joined$obs_grain_c)
  obs_range <- obs_max - obs_min
  tolerance <- 0.05 * obs_range

  below_floor    <- pmax(0, (obs_min - tolerance) - joined$grain_c)
  above_ceil     <- pmax(0, joined$grain_c - (obs_max + tolerance))
  n_exceedance   <- sum(below_floor > 0 | above_ceil > 0)
  exceedance_mag <- sum(below_floor + above_ceil) / obs_range
  range_penalty  <- exceedance_mag * 10

  # -------------------------------------------------------------------
  # 4. COMBINED OBJECTIVE
  # -------------------------------------------------------------------

  objective <- rmse *
    (1 +
       config$weights$bias    * abs(bias_norm)    +
       config$weights$spatial * (1 - spatial_cor) +
       config$weights$anomaly * (1 - anomaly_cor)
    ) + range_penalty

  # Log
  log_entry <- data.frame(
    timestamp    = Sys.time(),
    PRDX1        = params$PRDX1,
    PPDF1        = params$PPDF1,
    PPDF2        = params$PPDF2,
    WSCOEF1      = params$WSCOEF1,
    WSCOEF2      = params$WSCOEF2,
    CLAYPG       = params$CLAYPG,
    rmse         = rmse,
    bias_norm    = bias_norm,
    anomaly_cor  = anomaly_cor,
    spatial_cor  = spatial_cor,
    range_penalty   = range_penalty,
    n_exceedance    = n_exceedance,
    objective       = objective,
    n_matched       = nrow(joined),
    n_anom_rms      = n_distinct(anomalies$RMNO)
  )
  write_csv(log_entry,
            file.path(config$results_dir, "calibration_log.csv"),
            append = file.exists(file.path(config$results_dir, "calibration_log.csv")))

  cat(sprintf(
    "[%s] RMSE: %.2f  Bias: %+.3f  Anom: %.3f  Spat: %.3f  Rng: %.1f(%d)  Obj: %.3f\n",
    format(Sys.time(), "%H:%M:%S"),
    rmse, bias_norm, anomaly_cor, spatial_cor, range_penalty, n_exceedance, objective
  ))

  return(objective)
}

# =============================================================================
# 9. PARAMETER BOUNDS
# =============================================================================
#            PRDX1  PPDF1  PPDF2  WSCOEF1  WSCOEF2  CLAYPG
lower <- c(  1.0,   13.0,  20.0,  0.20,    6.0,     1.5 )
upper <- c(  6.5,   20.0,  32.0,  0.50,   30.0,     4.5 )
param_names_opt <- c("PRDX1","PPDF1","PPDF2","WSCOEF1","WSCOEF2","CLAYPG")

# =============================================================================
# 10. START CLUSTER AND RUN DEOPTIM
# =============================================================================

cat("\nStarting persistent cluster...\n")
cl <- makeCluster(config$n_cores, type = "PSOCK")

clusterExport(cl, varlist = c(
  "run_daycent_rm", "write_crop100", "extract_yield_rm", "config", "rm_folders"
))
clusterEvalQ(cl, library(tidyverse))

cat(sprintf("Cluster ready: %d workers\n", config$n_cores))
cat(sprintf("Population size: %d  |  Max iterations: %d\n",
            config$deoptim_NP, config$deoptim_iter))
cat(sprintf("RM folders: %d\n\n", length(rm_folders)))

set.seed(42)

# Seed initial population around best solution from prior run (89 generations)
# PRDX1=4.984, PPDF1=18.291, PPDF2=26.341, WSCOEF1=0.260, WSCOEF2=7.174, CLAYPG=2
best_known <- c(4.984, 18.291, 26.341, 0.260, 7.174, 2.0)
noise_frac  <- 0.10
ranges      <- upper - lower

initialpop <- matrix(
  best_known, nrow = config$deoptim_NP, ncol = length(lower), byrow = TRUE
) + matrix(
  runif(config$deoptim_NP * length(lower), -1, 1) * rep(ranges * noise_frac, each = config$deoptim_NP),
  nrow = config$deoptim_NP
)
for (j in seq_along(lower)) {
  initialpop[, j] <- pmax(lower[j], pmin(upper[j], initialpop[, j]))
}
initialpop[1, ] <- best_known

deoptim_result <- tryCatch({
  DEoptim(
    fn      = objective_fn,
    lower   = lower,
    upper   = upper,
    control = DEoptim.control(
      NP           = config$deoptim_NP,
      itermax      = config$deoptim_iter,
      CR           = config$deoptim_CR,
      F            = config$deoptim_F,
      trace        = 10,
      initialpop   = initialpop,
      parallelType = 0
    )
  )
}, error = function(e) {
  stopCluster(cl)
  stop(e)
})

stopCluster(cl)
cat("Cluster stopped\n")

# =============================================================================
# 11. EXTRACT AND SAVE RESULTS
# =============================================================================

best_params <- c(
  list(
    PRDX1   = deoptim_result$optim$bestmem[1],
    PPDF1   = deoptim_result$optim$bestmem[2],
    PPDF2   = deoptim_result$optim$bestmem[3],
    WSCOEF1 = deoptim_result$optim$bestmem[4],
    WSCOEF2 = deoptim_result$optim$bestmem[5],
    CLAYPG  = round(deoptim_result$optim$bestmem[6])
  ),
  fixed_params
)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat(sprintf("Best objective value: %.4f\n\n", deoptim_result$optim$bestval))

# Warn if any parameter resting on a bound
.on_bound <- (abs(deoptim_result$optim$bestmem - lower) < 1e-6) |
             (abs(deoptim_result$optim$bestmem - upper) < 1e-6)
if (any(.on_bound)) {
  cat(sprintf("WARNING: parameter(s) on a bound: %s\n\n",
              paste(param_names_opt[.on_bound], collapse = ", ")))
}

cat(sprintf("  PRDX(1)   = %.4f\n", best_params$PRDX1))
cat(sprintf("  PPDF(1)   = %.4f\n", best_params$PPDF1))
cat(sprintf("  PPDF(2)   = %.4f\n", best_params$PPDF2))
cat(sprintf("  WSCOEF(1) = %.4f\n", best_params$WSCOEF1))
cat(sprintf("  WSCOEF(2) = %.4f\n", best_params$WSCOEF2))
cat(sprintf("  CLAYPG    = %d\n",   as.integer(best_params$CLAYPG))  )

# Save parameter table
best_params_df <- data.frame(
  parameter = c("PRDX1","PPDF1","PPDF2","WSCOEF1","WSCOEF2","CLAYPG",
                "HIMAX","HIWSF","DDBASE","BASETEMP",
                "EMAX","TMPGERM","FRTC1","FRTC3","PLTMRF","HIMON1","CFRTCW1"),
  status    = c(rep("optimized", 6), rep("fixed_literature", 11)),
  value     = c(best_params$PRDX1, best_params$PPDF1, best_params$PPDF2,
                best_params$WSCOEF1, best_params$WSCOEF2, best_params$CLAYPG,
                best_params$HIMAX, best_params$HIWSF, best_params$DDBASE,
                best_params$BASETEMP, best_params$EMAX, best_params$TMPGERM,
                best_params$FRTC1, best_params$FRTC3, best_params$PLTMRF,
                best_params$HIMON1, best_params$CFRTCW1)
)
write_csv(best_params_df, file.path(config$results_dir, "best_parameters.csv"))

# Write calibrated crop.100 — use master template from Copy_Folder as base
template_crop100 <- file.path(config$template_dir, config$crop100_name)
write_crop100(
  params        = best_params,
  template_path = template_crop100,
  output_path   = file.path(config$results_dir, "crop_calibrated.100"),
  crop_name     = config$crop_id
)
cat(sprintf("Calibrated crop.100 saved to: %s\n",
            file.path(config$results_dir, "crop_calibrated.100")))

# =============================================================================
# 12. FINAL EVALUATION WITH BEST PARAMETERS
# =============================================================================

cat("\nRunning final evaluation with best parameters...\n")

cl <- makeCluster(config$n_cores, type = "PSOCK")
clusterExport(cl, varlist = c(
  "run_daycent_rm", "write_crop100", "extract_yield_rm", "config",
  "rm_folders", "best_params"
))
clusterEvalQ(cl, library(tidyverse))

final_sim <- run_all_rms(best_params, config)
stopCluster(cl)

final_joined <- final_sim %>%
  inner_join(obs_all %>% rename(obs_grain_c = grain_c), by = c("RMNO", "Year"))

rm_final <- final_joined %>%
  group_by(RMNO) %>%
  summarise(
    sim_mean = mean(grain_c,     na.rm = TRUE),
    obs_mean = mean(obs_grain_c, na.rm = TRUE),
    bias     = sim_mean - obs_mean,
    .groups  = "drop"
  )

r_spatial <- cor(rm_final$sim_mean, rm_final$obs_mean)
rmse_f    <- sqrt(mean(rm_final$bias^2))
bias_f    <- mean(rm_final$bias)

cat(sprintf("\n--- Final RM-mean evaluation (n=%d RMs) ---\n", nrow(rm_final)))
cat(sprintf("Spatial r : %.3f\n", r_spatial))
cat(sprintf("RMSE      : %.1f g C m-2\n", rmse_f))
cat(sprintf("Mean bias : %.1f g C m-2\n", bias_f))

write_csv(rm_final, file.path(config$results_dir, "rm_sim_vs_obs_calibrated.csv"))

# Plot
p <- ggplot(rm_final, aes(x = obs_mean, y = sim_mean)) +
  geom_point(alpha = 0.5, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "lm", colour = "black", se = TRUE) +
  labs(
    title    = "DayCent RM-Level Field Peas Calibration: Simulated vs Observed",
    subtitle = sprintf("Spatial r=%.3f  RMSE=%.1f  Bias=%.1f g C m⁻²  (n=%d RMs)",
                       r_spatial, rmse_f, bias_f, nrow(rm_final)),
    x = "Observed mean grain C (g C m⁻²)",
    y = "Simulated mean grain C (g C m⁻²)"
  ) +
  theme_bw()

ggsave(file.path(config$results_dir, "rm_calibration_plot.png"),
       p, width = 7, height = 6, dpi = 150)

cat(sprintf("\nAll outputs written to %s\n", config$results_dir))
