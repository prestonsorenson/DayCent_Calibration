# =============================================================================
# DayCent Perennial Grassland / Pasture Calibration Pipeline
# Saskatchewan - DDcentEVI_rev491 (weather-driven, PRDX mode)
# =============================================================================
# Calibration target: MODIS EVI-derived peak aboveground biomass (g C m-2)
#   obs = SK_pasture_abg_c_peak_by_LLD.csv  (peak growing-season standing crop,
#         winsorized to 5th/95th percentile; observed max = 367.4 g C m-2)
#   sim = agcacc (DayCent accumulated aboveground C at end of growing season, g C m-2)
#
# REDUCED PARAMETER SET — same 7-parameter framework as annual crops, with the
# remaining parameters FIXED at literature values for Saskatchewan perennial
# grassland (reference species: smooth brome grass, Bromus inermis, and
# crested wheatgrass, Agropyron cristatum).
#
# WSCOEF NOTE: WSCOEF(1) and WSCOEF(2) are included in the optimized set.
# With CLAYPG=5 (deep perennial rooting zone to ~80 cm), maxrwcf will likely
# be higher than under CLAYPG=2 for annuals, so WSCOEF(1) bounds are set
# higher (0.60-0.90) to allow the inflection to sit above where maxrwcf
# typically lives. The optimizer will find the appropriate stress level.
#
# OPTIMIZED (6 parameters):
#   PRDX(1)   - production scalar / light use efficiency analog
#   PPDF(1)   - temperature optimum for production (°C)
#   PPDF(2)   - maximum temperature for production (°C)
#   CLAYPG    - rooting-zone layers, DISCRETE integer (see bounds vector below)
#   (PPDF(3)/PPDF(4) no longer calibrated — left at their crop.100 values)
#   WSCOEF(1) - water stress inflection
#   WSCOEF(2) - water stress curve steepness
#
# FIXED at literature values (Saskatchewan perennial grassland):
#   HIMAX     = 0.50   structural parameter; not a true harvest index for
#                      perennial; mid-range of the 0.4-0.7 structural bounds
#   HIWSF     = 0.85   structural; mild water-stress sensitivity
#   DDBASE    = 1500   GDD to end of growing season at BASETEMP=3C; smooth
#                      brome growing season in SK is ~DOY 100-260 (~160 days)
#   BASETEMP  = 3.0    C3 grass base temperature; smooth brome begins growth
#                      at ~3°C soil temperature (cooler than cereal base 0°C)
#   EMAX      = 0.70   ET scalar; deep-rooted perennial grasses have higher
#                      actual ET than shallow-rooted annuals but still
#                      constrained by dryland Prairie conditions
#   TMPGERM   = 3.0    green-up temperature; smooth brome resumes growth at
#                      ~3-4°C in spring
#                      wheatgrass root to ~1.0-1.5 m; CLAYPG=5 captures the
#                      upper portion of this rooting zone
#   FRTC(1)   = 0.55   initial root:shoot C allocation; perennial grasses
#                      allocate more C belowground than annuals
#   FRTC(3)   = 150    GDD for allocation shift; wider than cereals given
#                      longer perennial growing season
#   PLTMRF    = 0.50   production scaling
#   HIMON(1)  = 2      months of water-stress accumulation
#   CFRTCW(1) = 0.50   root allocation increase under drought; slightly higher
#                      than annuals reflecting perennial drought adaptation
# =============================================================================

# -----------------------------------------------------------------------------
# CALIBRATION TARGET — PERENNIAL GRASS BIOMASS (NOT GRAIN, NOT FERTILIZED)
# -----------------------------------------------------------------------------
# SKGRASS is calibrated against MODIS-GPP-derived aboveground biomass (agcacc),
# not grain, and SK grassland/pasture is essentially unfertilized — production
# is climate/water-limited, so the cereal "non-limiting N" regime does NOT apply.
# HIMAX/HIWSF act on grain HI and are immaterial to the agcacc target.
# NOTE: CLAYPG is calibrated over {4,5,6} here (deep grass rooting, ~0-80 cm at
# CLAYPG=5), NOT the {2,3,4} used for the shallow-rooted annuals.
# NOTE: SKGRASS lacks a WSCOEF field in the G1-G5 lineage; WSCOEF(1)/(2) in the
# free set may be inert (write_crop100 will warn if the line is absent).
# -----------------------------------------------------------------------------

library(tidyverse)
library(DEoptim)
library(parallel)
library(fs)

# =============================================================================
# 1. CONFIGURATION — edit these paths for your system
# =============================================================================

config <- list(
  
  # Paths
  daycent_exe    = "/home/prs352/Daycent/Linux_Version_491/DDcentEVI_rev491",
  template_dir   = "/home/prs352/Daycent/Crop_Optimization/Models/Perennial/Copy_Folder",
  sites_dir      = "/home/prs352/Daycent/Crop_Optimization/Models/Perennial/qrts",
  results_dir    = "/home/prs352/Daycent/Crop_Optimization/Models/Perennial/results",
  obs_data_path  = "/home/prs352/Daycent/Crop_Optimization/Crop_Yield/SK_pasture_abg_c_peak_by_LLD.csv",
  crop100_name   = "crop.100",
  schedule_name  = "Crop_Opt.sch",
  output_prefix  = "harvest",                    # DayCent -n argument
  crop_id        = "SKGRASS",  # Saskatchewan perennial mixed grass
  
  # Clustering
  n_clusters     = 125,    # representative sites across Saskatchewan
  
  # Optimization
  n_cores        = 125,   # use --max-connections=256 flag when launching Rscript
  deoptim_NP     = 60,    # population size — 10x n_params (7 parameters)
  deoptim_iter   = 126,
  deoptim_CR     = 0.9,
  deoptim_F      = 0.8,
  
  # Objective function weights — all penalties multiply RMSE.
  # ANOMALY disabled (weight 0): consistent with annual crop framework.
  # Spatial weight elevated — strong Brown/Dark Brown/Black zone gradient
  # in SK pasture productivity is the primary signal to capture.
  weights = list(
    bias    = 0.3,
    spatial = 1.2,   # elevated — pasture productivity has strong spatial gradient
    anomaly = 0.0    # DISABLED — consistent with annual crop framework
  )
)

dir_create(config$results_dir)

# Clear stale calibration outputs — prior logs were produced under the broken
# HIMAX / N-starved regime and must not be appended to or reused.
for (.f in c("calibration_log.csv", "best_parameters.csv")) {
  .fp <- file.path(config$results_dir, .f)
  if (file.exists(.fp)) file.remove(.fp)
}

# =============================================================================
# FIXED PARAMETERS — Saskatchewan perennial grassland literature values
# =============================================================================
# Reference species: smooth brome (Bromus inermis) and crested wheatgrass
# (Agropyron cristatum) — dominant managed perennial grasses in SK.
#
#   HIMAX     = 0.50  Structural parameter; not a true harvest index for
#                     perennial grasses. Mid-range of the 0.4-0.7 structural bounds.
#   HIWSF     = 0.85  Structural; mild water-stress sensitivity.
#   DDBASE    = 1500  GDD to end of growing season at BASETEMP=3C. Smooth
#                     brome green-up ~DOY 100, senescence ~DOY 260 in SK;
#                     at base 3C over ~160-day season ≈ 1400-1600 GDD.
#   BASETEMP  = 3.0   C3 grass base temperature; smooth brome and crested
#                     wheatgrass begin growth at ~3°C (cooler than cereal base 0°C).
#   EMAX      = 0.70  ET scalar; deep-rooted perennial grasses have higher potential
#                     ET than shallow-rooted annuals but dryland Prairie conditions
#                     still constrain actual ET substantially.
#   TMPGERM   = 3.0   Green-up temperature; smooth brome and crested wheatgrass
#                     resume spring growth at ~3-4°C.
#   CLAYPG    = 5     Rooting zone layers (~0-80 cm). Smooth brome and crested
#                     wheatgrass root to ~1.0-1.5 m; CLAYPG=5 captures the
#                     ecologically active upper root zone.
#   FRTC(1)   = 0.55  Initial root:shoot C allocation; perennial grasses allocate
#                     substantially more C belowground than annuals
#                     (Gill & Jackson 2000; root:shoot ratios 2-5 for Prairie grasses).
#   FRTC(3)   = 150   GDD for allocation shift; wider than cereals reflecting the
#                     longer perennial growing season.
#   PLTMRF    = 0.50  Production scaling; mid-range default.
#   HIMON(1)  = 2     Months of water-stress accumulation.
#   CFRTCW(1) = 0.50  Root allocation increase under water stress; slightly higher
#                     than annuals reflecting perennial drought-adaptation strategy.
fixed_params <- list(
  HIMAX    = 0.50,
  HIWSF    = 0.50,
  DDBASE   = 1500.0,
  BASETEMP = 3.0,
  EMAX     = 0.70,
  TMPGERM  = 3.0,
  FRTC1    = 0.55,
  FRTC3    = 150.0,
  PLTMRF   = 0.50,
  HIMON1   = 2,
  CFRTCW1  = 0.5
)

# =============================================================================
# 2. LOAD AND PREPARE OBSERVED DATA
# =============================================================================

obs_raw <- read_csv(config$obs_data_path, show_col_types = FALSE)
subsets=list.dirs('/home/prs352/Daycent/Crop_Optimization/Models/Perennial/qrts', full.names=FALSE)
obs_raw=obs_raw[obs_raw$LEGALLAND %in% subsets, ]

# One observation per site per year — rename column, remove duplicates,
# and cap implausible outliers. The peak-biomass input is already winsorized
# to the 5th/95th percentile (observed max = 367.4 g C m-2), so the cap is a
# safety guard set at the observed maximum.
cap_value <- 368

obs_summary <- obs_raw %>%
  rename(obs_grain_c = abg_c_peak_gmsq) %>%
  mutate(Year = as.integer(Year)) %>%
  select(LEGALLAND, Year, obs_grain_c) %>%
  distinct(LEGALLAND, Year, .keep_all = TRUE) %>%
  filter(obs_grain_c <= cap_value)

cat(sprintf("Removed %d implausible records (>%.0f g C m-2)\n",
            nrow(obs_raw) - nrow(obs_summary), cap_value))

cat(sprintf("Loaded %d observed site-year records from %d unique sites\n",
            nrow(obs_summary),
            n_distinct(obs_summary$LEGALLAND)))

# =============================================================================
# 3. LOAD REPRESENTATIVE SITES FROM EXISTING FOLDERS
# =============================================================================
# Sites were clustered and folders created externally.
# Reads folder names directly from sites_dir.

site_folders <- list.dirs(config$sites_dir, full.names = FALSE, recursive = FALSE)
site_folders <- site_folders[nchar(site_folders) > 0]  # drop empty strings

rep_sites <- data.frame(LEGALLAND = site_folders, stringsAsFactors = FALSE)

cat(sprintf("Found %d site folders in %s\n", nrow(rep_sites), config$sites_dir))

if (nrow(rep_sites) == 0) {
  stop(sprintf("No site folders found in: %s\nCheck that config$sites_dir is correct.",
               config$sites_dir))
}

# =============================================================================
# 4. SITE FOLDERS ALREADY EXIST — NO SETUP NEEDED
# =============================================================================

# =============================================================================
# 5. WRITE crop.100 WITH CALIBRATION PARAMETERS
# =============================================================================

write_crop100 <- function(params, template_path, output_path, crop_name = "SKGRASS") {
  
  lines <- readLines(template_path)
  
  # Locate crop block
  crop_start   <- grep(paste0("^", crop_name, "\\b"), lines)
  if (length(crop_start) == 0) stop(paste("Crop", crop_name, "not found"))
  
  crop_headers <- grep("^[A-Z][A-Z0-9]+\\s+\\S", lines)
  next_header  <- crop_headers[crop_headers > crop_start][1]
  crop_end     <- ifelse(is.na(next_header), length(lines), next_header - 1)
  
  # Replace a single parameter value within the crop block
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
  
  # Apply all calibration parameters
  if (!is.null(params$PRDX1))   lines <- replace_param(lines, crop_start, crop_end, "PRDX\\(1\\)",   params$PRDX1)
  if (!is.null(params$PPDF1))   lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(1\\)",   params$PPDF1)
  if (!is.null(params$PPDF2))   lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(2\\)",   params$PPDF2)
  if (!is.null(params$WSCOEF1)) lines <- replace_param(lines, crop_start, crop_end, "WSCOEF\\(1\\)", params$WSCOEF1)
  if (!is.null(params$WSCOEF2)) lines <- replace_param(lines, crop_start, crop_end, "WSCOEF\\(2\\)", params$WSCOEF2)
  if (!is.null(params$HIMAX))    lines <- replace_param(lines, crop_start, crop_end, "HIMAX",    params$HIMAX)
  if (!is.null(params$HIWSF))    lines <- replace_param(lines, crop_start, crop_end, "HIWSF",    params$HIWSF)
  if (!is.null(params$DDBASE))   lines <- replace_param(lines, crop_start, crop_end, "DDBASE",   params$DDBASE)
  if (!is.null(params$BASETEMP)) lines <- replace_param(lines, crop_start, crop_end, "BASETEMP(?!\\()", params$BASETEMP)
  if (!is.null(params$EMAX))     lines <- replace_param(lines, crop_start, crop_end, "EMAX",     params$EMAX)
  if (!is.null(params$TMPGERM))  lines <- replace_param(lines, crop_start, crop_end, "TMPGERM",  params$TMPGERM)
  if (!is.null(params$CLAYPG))   lines <- replace_param(lines, crop_start, crop_end, "CLAYPG",         round(params$CLAYPG))
  if (!is.null(params$FRTC1))    lines <- replace_param(lines, crop_start, crop_end, "FRTC\\(1\\)",      params$FRTC1)
  if (!is.null(params$FRTC3))    lines <- replace_param(lines, crop_start, crop_end, "FRTC\\(3\\)",      params$FRTC3)
  if (!is.null(params$PLTMRF))   lines <- replace_param(lines, crop_start, crop_end, "PLTMRF",         params$PLTMRF)
  if (!is.null(params$HIMON1))   lines <- replace_param(lines, crop_start, crop_end, "HIMON\\(1\\)",     round(params$HIMON1))
  if (!is.null(params$CFRTCW1))  lines <- replace_param(lines, crop_start, crop_end, "CFRTCW\\(1\\)",   params$CFRTCW1)
  if (!is.null(params$PPDF3))    lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(3\\)",      params$PPDF3)
  if (!is.null(params$PPDF4))    lines <- replace_param(lines, crop_start, crop_end, "PPDF\\(4\\)",      params$PPDF4)
  
  writeLines(lines, output_path)
  invisible(NULL)
}

# =============================================================================
# 6. RUN DAYCENT FOR ONE SITE
# =============================================================================

run_daycent_site <- function(site_id, params, config) {
  
  site_dir      <- file.path(config$sites_dir, site_id)
  crop100_path  <- file.path(site_dir, config$crop100_name)
  output_csv    <- file.path(site_dir, "harvest.csv")
  
  # Write parameters to this site's crop.100
  tryCatch(
    write_crop100(params, crop100_path, crop100_path, config$crop_id),
    error = function(e) {
      warning(paste("write_crop100 failed for", site_id, ":", e$message))
      return(NULL)
    }
  )
  
  # Remove ALL stale DayCent output files before each run
  # DayCent aborts with "binary file exists" if any output from a prior
  # run is present — delete all .bin, .lis, .csv outputs with the prefix
  stale_files <- list.files(
    site_dir,
    pattern    = paste0("^", config$output_prefix, "\\.(bin|lis|csv|out)$"),
    full.names = TRUE
  )
  # Also catch numbered variants e.g. "harvest 1.bin"
  stale_files <- c(stale_files, list.files(
    site_dir,
    pattern    = paste0("^", config$output_prefix, ".*\\.(bin|lis|csv|out)$"),
    full.names = TRUE
  ))
  if (length(stale_files) > 0) file.remove(unique(stale_files))
  
  # Run DayCent from site directory
  cmd <- paste0(
    "cd ", shQuote(site_dir),
    " && ", shQuote(config$daycent_exe),
    " -s ", config$schedule_name,
    " -n ", config$output_prefix,
    " > daycent.log 2>&1"
  )
  
  exit_code <- system(cmd)
  
  if (exit_code != 0 || !file.exists(output_csv)) {
    warning(paste("DayCent failed for site:", site_id,
                  "- check", file.path(site_dir, "daycent.log")))
    return(NULL)
  }
  
  # Delete harvest.bin after successful run — prevents "binary file exists"
  # error on next evaluation without needing to clean up beforehand
  bin_file <- file.path(site_dir, paste0(config$output_prefix, ".bin"))
  if (file.exists(bin_file)) file.remove(bin_file)
  
  return(output_csv)
}

# =============================================================================
# 7. EXTRACT SIMULATED YIELDS
# =============================================================================

extract_yield <- function(output_csv, site_id) {
  
  if (is.null(output_csv) || !file.exists(output_csv)) return(NULL)
  
  tryCatch({
    read_csv(output_csv, show_col_types = FALSE) %>%
      mutate(
        site_id = site_id,
        Year    = as.integer(floor(time))   # DayCent time is decimal year
      ) %>%
      rename(grain_c = agcacc) %>%   # agcacc = accumulated aboveground C — standing biomass at end of growing season
      select(site_id, Year, grain_c)
  }, error = function(e) {
    warning(paste("Could not read output for site:", site_id, "-", e$message))
    return(NULL)
  })
}

# =============================================================================
# 8. RUN ALL SITES IN PARALLEL
# =============================================================================
# Uses a persistent PSOCK cluster (cl) started in Section 11 before DEoptim.
# params exported to workers on each call via clusterExport.

run_all_sites <- function(params, rep_sites, config) {
  
  # Export current params to all workers
  clusterExport(cl, varlist = "params", envir = environment())
  
  results <- parLapply(
    cl,
    rep_sites$LEGALLAND,
    function(site_id) {
      out_csv <- run_daycent_site(site_id, params, config)
      extract_yield(out_csv, site_id)
    }
  )
  
  # Combine all site results, dropping failed runs
  bind_rows(Filter(Negate(is.null), results))
}

# =============================================================================
# 9. OBJECTIVE FUNCTION
# =============================================================================

objective_fn <- function(param_vec) {
  
  # Unpack the 6 optimized parameters and merge with fixed literature values.
  # CLAYPG is discrete: sampled continuously and rounded to an integer.
  # PPDF(3)/PPDF(4) are NOT set here, so write_crop100 leaves their crop.100 values.
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
  
  # Run all representative sites
  simulated <- tryCatch(
    run_all_sites(params, rep_sites, config),
    error = function(e) {
      warning(paste("run_all_sites failed:", e$message))
      return(NULL)
    }
  )
  
  if (is.null(simulated) || nrow(simulated) == 0) return(Inf)
  
  # Join simulated to observed on site + year
  joined <- simulated %>%
    rename(LEGALLAND = site_id) %>%
    inner_join(obs_summary, by = c("LEGALLAND", "Year"))
  
  if (nrow(joined) < 800) return(Inf)   # require at least 800 matched site-years (125 sites × 26 yrs = 3250 available)
  
  # ---------------------------------------------------------------------------
  # PERFORMANCE CRITERIA
  # ---------------------------------------------------------------------------
  
  # 1. RMSE — average prediction error across all site-years
  rmse <- sqrt(mean((joined$grain_c - joined$obs_grain_c)^2))
  
  # 2. BIAS — mean signed error (positive = model overestimates)
  #    Normalised by observed mean so it is unitless and comparable across crops
  obs_grand_mean <- mean(joined$obs_grain_c)
  bias_norm <- mean(joined$grain_c - joined$obs_grain_c) / obs_grand_mean
  
  # 3. ANOMALY CORRELATION
  #    For each site, subtract the site mean from both simulated and observed NPP,
  #    then normalise by the site standard deviation (z-score). This removes both the
  #    site-level mean productivity and the site-level variance, so every site
  #    contributes equally regardless of its position on the spatial gradient.
  #    Pooling all site-year standardised anomalies and correlating tests whether
  #    the model correctly ranks years as above or below average across the full
  #    network — drought/heat suppression years are below average and good years
  #    are above average in aggregate, independent of site mean or variance.
  anomalies <- joined %>%
    group_by(LEGALLAND) %>%
    filter(n() >= 4) %>%
    mutate(
      obs_anom = (obs_grain_c - mean(obs_grain_c)) / sd(obs_grain_c),
      sim_anom = (grain_c     - mean(grain_c))     / sd(grain_c)
    ) %>%
    ungroup()
  
  anomaly_cor <- if (nrow(anomalies) >= 50) {
    tryCatch(
      cor(anomalies$sim_anom, anomalies$obs_anom, method = "pearson"),
      error = function(e) 0
    )
  } else 0
  
  # 4. SPATIAL CORRELATION — simulated site means match observed site means
  site_means <- joined %>%
    group_by(LEGALLAND) %>%
    summarise(
      sim_mean = mean(grain_c),
      obs_mean = mean(obs_grain_c),
      .groups  = "drop"
    )
  
  spatial_cor <- tryCatch(
    cor(site_means$sim_mean, site_means$obs_mean, method = "pearson"),
    error = function(e) 0
  )
  
  # 5. RANGE EXCEEDANCE PENALTY
  #    Severe additive penalty when any simulated site-year NPP falls outside
  #    the global observed min/max. 5% tolerance buffer before penalty fires.
  #    Each unit of observed range exceeded contributes a x10 addition to the
  #    objective, firing independently of weights.
  obs_min   <- min(joined$obs_grain_c)
  obs_max   <- max(joined$obs_grain_c)
  obs_range <- obs_max - obs_min
  tolerance <- 0.05 * obs_range
  
  below_floor <- pmax(0, (obs_min - tolerance) - joined$grain_c)
  above_ceil  <- pmax(0, joined$grain_c - (obs_max + tolerance))
  
  n_exceedance   <- sum(below_floor > 0 | above_ceil > 0)
  exceedance_mag <- sum(below_floor + above_ceil) / obs_range
  
  range_penalty <- exceedance_mag * 10
  
  # ---------------------------------------------------------------------------
  # COMBINED OBJECTIVE
  # ---------------------------------------------------------------------------
  # Weights control the relative importance of each criterion.
  # All penalty terms are positive and push the objective higher (worse).
  #
  #   RMSE component:        base prediction error (g C m-2 yr-1 NPP)
  #   Bias penalty:          |bias_norm| — penalises systematic over/under
  #   Spatial penalty:       (1 - spatial_cor) — penalises wrong spatial NPP gradient
  #   Anomaly penalty:       (1 - anomaly_cor) — penalises failure to reproduce
  #                          site-normalised year-to-year anomalies across all sites
  #   Range penalty:         additive — severe penalty for any simulated NPP outside
  #                          observed min/max (5% tolerance); fires independently of weights
  #
  # Tune weights in config$weights to emphasise different criteria.
  
  objective <- rmse *
    (1 +
       config$weights$bias    * abs(bias_norm)      +
       config$weights$spatial * (1 - spatial_cor)   +
       config$weights$anomaly * (1 - anomaly_cor)
    ) + range_penalty
  
  # Log progress to file
  # Log progress to file (only optimized parameters vary)
  log_entry <- data.frame(
    timestamp   = Sys.time(),
    PRDX1       = params$PRDX1,
    PPDF1       = params$PPDF1,
    PPDF2       = params$PPDF2,
    WSCOEF1     = params$WSCOEF1,
    WSCOEF2     = params$WSCOEF2,
    CLAYPG      = params$CLAYPG,
    rmse        = rmse,
    bias_norm   = bias_norm,
    anomaly_cor  = anomaly_cor,
    spatial_cor  = spatial_cor,
    range_penalty   = range_penalty,
    n_exceedance    = n_exceedance,
    objective    = objective,
    n_matched    = nrow(joined),
    n_anom_sites = n_distinct(anomalies$LEGALLAND)
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
# 10. PARAMETER BOUNDS — 6 OPTIMIZED PARAMETERS
# =============================================================================
#            PRDX1  PPDF1  PPDF2  WSCOEF1  WSCOEF2  CLAYPG
# Literature justification (reference species: smooth brome, crested wheatgrass):
#   PRDX1:    2.0-5.0  — RUE for smooth brome and cool-season Prairie forages
#                        ~2.0-4.5 g DM/MJ; productive C3 grass with documented
#                        RUE above 2.0 under adequate moisture.
#   PPDF1:    14-20°C  — smooth brome photosynthesis optimum ~16-20°C; C3
#                        cool-season grass; slightly warmer than annual cereals.
#   PPDF2:    28-38°C  — smooth brome ceiling; C3 grasses retain some function
#                        to ~38°C; wider upper than annual crops reflecting
#                        perennial heat acclimation.
#   CLAYPG:   DISCRETE integer (see bounds vector) — rooting-zone layer count,
#                        sampled continuously and rounded; interacts with WSCOEF1.
#   WSCOEF1:  0.60-0.90 — WIDE band appropriate for CLAYPG=5. With deep rooting
#                        zone, maxrwcf is likely higher than at CLAYPG=2 for
#                        annuals (deep layers near FC). Inflection placed in the
#                        upper maxrwcf range so stress can activate during dry years.
#   WSCOEF2:  7.0-15.0 — wider than the narrow annual band; allows the optimizer
#                        to tune the sigmoid steepness more freely given the less-
#                        constrained stress regime of a deep-rooted perennial.
lower <- c( 1.0, 14.0, 28.0, 0.60, 7.0, 3.5 )
upper <- c( 5.0, 20.0, 38.0, 0.90, 15.0, 6.5 )
param_names_opt <- c("PRDX1","PPDF1","PPDF2","WSCOEF1","WSCOEF2","CLAYPG")

# =============================================================================
# 10b. OPTIONAL: SOBOL SENSITIVITY ANALYSIS
# =============================================================================
# Run this first (before DEoptim) to confirm which parameters actually drive
# yield variance in your system. Requires ~500-1000 DayCent evaluations.
# Comment out if you want to go straight to DEoptim.

run_sensitivity <- function() {
  library(sensitivity)
  
  n_sobol <- 100
  n_par   <- 6
  
  X1 <- data.frame(matrix(runif(n_par * n_sobol), ncol = n_par))
  X2 <- data.frame(matrix(runif(n_par * n_sobol), ncol = n_par))
  
  names(X1) <- names(X2) <- param_names_opt
  
  # Scale from [0,1] to actual parameter ranges
  scale_params <- function(X) {
    for (i in seq_along(param_names_opt)) {
      X[[i]] <- lower[i] + X[[i]] * (upper[i] - lower[i])
    }
    X
  }
  X1 <- scale_params(X1)
  X2 <- scale_params(X2)
  
  sob <- sobol2002(model = NULL, X1 = X1, X2 = X2, nboot = 50)
  
  # Evaluate model at all Sobol sample points
  y <- apply(sob$X, 1, objective_fn)
  tell(sob, y)
  
  # Save sensitivity results
  sens_results <- data.frame(
    parameter = param_names_opt,
    S1        = sob$S$original,   # first-order indices
    ST        = sob$T$original    # total-order indices
  ) %>% arrange(desc(ST))
  
  write_csv(sens_results,
            file.path(config$results_dir, "sensitivity_indices.csv"))
  
  cat("\nSobol Sensitivity Results:\n")
  print(sens_results)
  
  return(sens_results)
}

# Uncomment to run sensitivity analysis first:
# sens_results <- run_sensitivity()

# =============================================================================
# 11. RUN DEOPTIM CALIBRATION
# =============================================================================
# Cluster started ONCE here and reused across all DEoptim evaluations.
# This prevents the RAM spike caused by repeatedly creating/destroying
# worker processes on every objective function call.

cat("\nStarting persistent cluster...\n")
cl <- makeCluster(config$n_cores, type = "PSOCK")

clusterExport(cl, varlist = c(
  "run_daycent_site", "write_crop100", "extract_yield", "config"
))
clusterEvalQ(cl, library(tidyverse))

cat(sprintf("Cluster ready: %d workers\n", config$n_cores))
cat(sprintf("Population size: %d  |  Max iterations: %d\n",
            config$deoptim_NP, config$deoptim_iter))
cat(sprintf("Representative sites: %d\n\n", nrow(rep_sites)))

set.seed(42)

deoptim_result <- tryCatch({
  DEoptim(
    fn      = objective_fn,
    lower   = lower,
    upper   = upper,
    control = DEoptim.control(
      NP       = config$deoptim_NP,
      itermax  = config$deoptim_iter,
      CR       = config$deoptim_CR,
      F        = config$deoptim_F,
      trace    = 10,
      parallelType = 0  # parallelism handled via persistent cl
    )
  )
}, error = function(e) {
  stopCluster(cl)
  stop(e)
})

stopCluster(cl)
cat("Cluster stopped\n")

# =============================================================================
# 12. EXTRACT AND SAVE RESULTS
# =============================================================================

best_params <- c(
  list(
    PRDX1   = deoptim_result$optim$bestmem[1],
    PPDF1   = deoptim_result$optim$bestmem[2],
    PPDF2   = deoptim_result$optim$bestmem[3],
    WSCOEF1 = deoptim_result$optim$bestmem[4],
    WSCOEF2 = deoptim_result$optim$bestmem[5],
    CLAYPG  = round(deoptim_result$optim$bestmem[6])   # discrete integer
  ),
  fixed_params
)

cat("\n=== CALIBRATION COMPLETE ===\n")
cat(sprintf("Best objective value: %.4f\n\n", deoptim_result$optim$bestval))

# Diagnostic: flag any optimized parameter resting on a bound (re-examine bounds
# / N regime if so — a bound-pinned PRDX usually signals the N rate, not RUE).
.on_bound <- (abs(deoptim_result$optim$bestmem - lower) < 1e-6) |
             (abs(deoptim_result$optim$bestmem - upper) < 1e-6)
if (any(.on_bound)) {
  cat(sprintf("WARNING: parameter(s) on a bound: %s\n\n",
              paste(param_names_opt[.on_bound], collapse = ", ")))
}
cat("Optimized parameters:\n")
cat(sprintf("  PRDX(1)   = %.4f\n", best_params$PRDX1))
cat(sprintf("  PPDF(1)   = %.4f\n", best_params$PPDF1))
cat(sprintf("  PPDF(2)   = %.4f\n", best_params$PPDF2))
cat(sprintf("  WSCOEF(1) = %.4f  (bounded 0.60-0.90)\n", best_params$WSCOEF1))
cat(sprintf("  WSCOEF(2) = %.4f  (bounded 7.0-15.0)\n", best_params$WSCOEF2))
cat(sprintf("  CLAYPG    = %d       (discrete integer)\n", as.integer(best_params$CLAYPG)))
cat("\nFixed parameters (literature / crop.100 defaults, not optimized):\n")
cat(sprintf("  HIMAX     = %.2f   HIWSF     = %.2f\n", best_params$HIMAX, best_params$HIWSF))
cat(sprintf("  DDBASE    = %.0f  BASETEMP  = %.1f\n", best_params$DDBASE, best_params$BASETEMP))
cat(sprintf("  EMAX      = %.2f   TMPGERM   = %.1f\n", best_params$EMAX, best_params$TMPGERM))
cat(sprintf("  HIMON(1)  = %d      PPDF(3,4) = crop.100 defaults\n", as.integer(best_params$HIMON1)))
cat(sprintf("  FRTC(1)   = %.2f   FRTC(3)   = %.0f\n", best_params$FRTC1, best_params$FRTC3))
cat(sprintf("  PLTMRF    = %.2f   CFRTCW(1) = %.2f\n", best_params$PLTMRF, best_params$CFRTCW1))

# Save best parameters (optimized + fixed flagged)
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
write_csv(best_params_df,
          file.path(config$results_dir, "best_parameters.csv"))

# Write calibrated crop.100 to results folder
write_crop100(
  params        = best_params,
  template_path = file.path(config$template_dir, config$crop100_name),
  output_path   = file.path(config$results_dir, "crop_calibrated.100"),
  crop_name     = config$crop_id
)

cat(sprintf("\nCalibrated crop.100 saved to: %s\n",
            file.path(config$results_dir, "crop_calibrated.100")))

# =============================================================================
# 13. VALIDATION
# =============================================================================
# Run final simulation with best parameters and compute full diagnostics.
# Use sites NOT included in rep_sites for a true out-of-sample validation.

validate_calibration <- function(best_params, rep_sites, obs_summary, config) {
  
  cat("\nRunning validation...\n")
  
  # Get all sites that were NOT used in calibration
  all_sites <- obs_summary %>%
    distinct(LEGALLAND) %>%
    filter(dir_exists(file.path(config$sites_dir, LEGALLAND)),
           !LEGALLAND %in% rep_sites$LEGALLAND)
  
  # Sample up to 200 validation sites
  val_sites <- all_sites %>% slice_sample(n = min(200, nrow(all_sites)))
  
  cat(sprintf("Validating on %d held-out sites\n", nrow(val_sites)))
  
  simulated <- run_all_sites(best_params, val_sites, config)
  
  if (nrow(simulated) == 0) {
    cat("No validation output produced\n")
    return(invisible(NULL))
  }
  
  joined <- simulated %>%
    rename(LEGALLAND = site_id) %>%
    inner_join(obs_summary, by = c("LEGALLAND", "Year"))
  
  # Diagnostics
  rmse  <- sqrt(mean((joined$grain_c - joined$obs_grain_c)^2))
  bias  <- mean(joined$grain_c - joined$obs_grain_c)
  r2    <- cor(joined$grain_c, joined$obs_grain_c)^2
  
  cat(sprintf("\nValidation results (n=%d site-years):\n", nrow(joined)))
  cat(sprintf("  RMSE : %.2f g C m-2\n", rmse))
  cat(sprintf("  Bias : %.2f g C m-2 (%s)\n", bias,
              ifelse(bias > 0, "model overestimates", "model underestimates")))
  cat(sprintf("  R²   : %.3f\n", r2))
  
  # Save validation output
  write_csv(joined,
            file.path(config$results_dir, "validation_results.csv"))
  
  # Plot observed vs simulated
  p <- ggplot(joined, aes(x = obs_grain_c, y = grain_c)) +
    geom_point(alpha = 0.3, size = 0.8, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "lm", colour = "black", se = TRUE) +
    labs(
      title    = "SKGRASS Perennial Calibration Validation",
      subtitle = sprintf("RMSE=%.1f  Bias=%.1f  R²=%.3f  (n=%d)",
                         rmse, bias, r2, nrow(joined)),
      x        = "Observed aboveground C (g C m⁻²)",
      y        = "Simulated aboveground C (g C m⁻²)"
    ) +
    theme_bw()
  
  ggsave(file.path(config$results_dir, "validation_plot.png"),
         p, width = 7, height = 6, dpi = 150)
  
  cat(sprintf("Validation plot saved to: %s\n",
              file.path(config$results_dir, "validation_plot.png")))
  
  return(joined)
}

# Run validation
val_results <- validate_calibration(best_params, rep_sites, obs_summary, config)