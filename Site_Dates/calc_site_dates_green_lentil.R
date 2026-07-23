# =============================================================================
# Site-Specific Seeding and Harvest Date Calculator
# Saskatchewan Canadian Prairie Green Lentil (GLENT) — DDcentEVI 4.91
# =============================================================================
# Reads weather.wth from each site folder, computes site-specific seeding
# and harvest DOY, then updates Crop_Opt.sch in each folder.
#
# Schedule lines updated:
#   1 105 CULT KILL  → 1 {seed_doy - 5} CULT KILL
#   1 105 CULT NDRIL → 1 {seed_doy}     CULT NDRIL
#   1 105 FERT N10U  → 1 {seed_doy}     FERT N10U
#   1 105 CROP GLENT → 1 {seed_doy}     CROP GLENT
#   1 105 PLTM       → 1 {seed_doy}     PLTM
#   1 273 LAST       → 1 {harvest_doy}  LAST
#   1 273 HARV G     → 1 {harvest_doy}  HARV G
#
# Phenology values (TMPGERM, DDBASE, BASETEMP) are FIXED literature values
# matching the green lentil calibration crop.100 (GLENT block) — not optimized.
# BASETEMP=0 here matches the calibration and the Seeding_Date_App GDD targets.
# =============================================================================

library(tidyverse)

# =============================================================================
# 1. CONFIGURATION — update with calibrated values after optimization
# =============================================================================

config_dates <- list(
  
  sites_dir     = "/home/preston/Projects/Daycent/Crop_Optimization/RM/Models/Green_Lentil/qrts",
  results_dir   = "/home/preston/Projects/Daycent/Crop_Optimization/Models/Green_Lentil/results",
  schedule_name = "Crop_Opt.sch",
  weather_name  = "weather.wth",
  
  # GDD / phenology values — MUST match the green lentil calibration crop.100 (BASETEMP=0)
  # Consistency note: the DayCent LENT calibration fixes BASETEMP=0.0, so GDD
  # accumulation here uses base 0C to stay internally consistent with the
  # Seeding_Date_App. Green lentil matures in ~95-105 days in SK (slightly later
  # than red lentil); at base 0C ≈ 1450 GDD0. Matches daycent_calibration_green_lentil_hpc.R.
  #
  # TMPGERM: air-temperature proxy; lentil germinates at ~2-4C (cold tolerant,
  #          can be seeded early). Use 2.0 to match the calibration.
  # TMPKILL: -5C killing frost; lentil hypogeal germination, regrows from scale
  #          nodes, survives -4 to -6C (Saskatchewan Pulse Growers).
  # BASETEMP2: 27C ceiling — green lentil indeterminate, heat-sensitive during
  #            extended flowering.
  TMPGERM   = 2.0,    # min germination temperature (°C) — air temp proxy
  TMPKILL   = -5.0,   # killing frost temperature (°C) — survives -4 to -6C
  DDBASE    = 1450.0, # GDD to maturity at BASETEMP=0 — green lentil ~95-105 days
  BASETEMP  = 0.0,    # base temperature for GDD accumulation (°C) — matches calibration
  BASETEMP2 = 27.0,   # ceiling temperature (°C) — heat sensitive during flowering
  
  # Seeding window — uniform DOY 105-166 from the Seeding_Date_App standard
  earliest_seed = 105,   # DOY 105 = April 15
  latest_seed   = 166,   # DOY 166 = June 15
  
  # Harvest constraints — green lentil matures later than red lentil
  earliest_harvest = 210,  # DOY 210 = July 29
  latest_harvest   = 280   # DOY 280 = October 7 — hard cutoff (frost risk)
)

dir.create(config_dates$results_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. READ WEATHER.WTH FILE
# =============================================================================
# DayCent .wth format: Day  Month  Year  Tmax  Tmin  Precip
# Header line is skipped

read_weather_wth <- function(wth_path) {
  
  if (!file.exists(wth_path)) {
    warning(paste("Weather file not found:", wth_path))
    return(NULL)
  }
  
  # Read skipping the header line
  # Read positionally — handles both 7-column and 10-column (extra drivers) files.
  # Columns 1-7 are always: day, month, year, doy, tmax, tmin, precip
  tryCatch({
    raw <- read.table(wth_path, header = FALSE, fill = TRUE)
    
    df <- raw %>%
      transmute(
        day       = as.integer(V1),
        month     = as.integer(V2),
        year      = as.integer(V3),
        doy       = as.integer(V4),
        tmax_c    = as.numeric(V5),
        tmin_c    = as.numeric(V6),
        precip_mm = as.numeric(V7)
      ) %>%
      mutate(
        date  = as.Date(sprintf("%04d-%02d-%02d", year, month, day)),
        tmean = (tmax_c + tmin_c) / 2
      ) %>%
      filter(!is.na(date))
    
    return(df)
    
  }, error = function(e) {
    warning(paste("Could not read weather file:", wth_path, "-", e$message))
    return(NULL)
  })
}

# =============================================================================
# 3. CALCULATE SEEDING DATE FOR ONE YEAR
# =============================================================================
# First date after earliest_seed where 5-day running mean temp >= TMPGERM

calc_seeding_doy <- function(year_data, config_dates) {
  
  candidates <- year_data %>%
    filter(doy >= config_dates$earliest_seed,
           doy <= config_dates$latest_seed) %>%
    arrange(doy) %>%
    mutate(tmean_5day = zoo::rollmean(tmean, k = 5, fill = NA, align = "right"))
  
  seed_row <- candidates %>%
    filter(!is.na(tmean_5day), tmean_5day >= config_dates$TMPGERM) %>%
    slice(1)
  
  # If no suitable date found use latest_seed as fallback
  if (nrow(seed_row) == 0) return(as.integer(config_dates$latest_seed))
  return(as.integer(seed_row$doy))
}

# =============================================================================
# 4. CALCULATE HARVEST DATE FOR ONE YEAR
# =============================================================================
# Earliest of: cumulative GDD >= DDBASE, killing frost, or latest_harvest cap

calc_harvest_doy <- function(year_data, seed_doy, config_dates) {
  
  if (is.na(seed_doy)) return(as.integer(config_dates$latest_harvest))
  
  growing <- year_data %>%
    filter(doy >= seed_doy) %>%
    arrange(doy) %>%
    mutate(
      # GDD with base and ceiling temperatures — matches DayCent internal calc
      # GDD with base and ceiling: cap tmax at BASETEMP2, floor tmin at BASETEMP,
      # take mean, subtract base, floor at 0. Standard Miller et al. (2001) form.
      tmax_capped  = pmin(tmax_c, config_dates$BASETEMP2),
      tmin_floored = pmax(tmin_c, config_dates$BASETEMP),
      gdd_daily    = pmax(0, (tmax_capped + tmin_floored) / 2 - config_dates$BASETEMP),
      gdd_cumul    = cumsum(gdd_daily)
    )
  
  # Maturity — when cumulative GDD exceeds DDBASE
  maturity_row <- growing %>%
    filter(gdd_cumul >= config_dates$DDBASE,
           doy      >= config_dates$earliest_harvest) %>%
    slice(1)
  
  # Killing frost — first day Tmin drops below TMPKILL after earliest_harvest
  frost_row <- growing %>%
    filter(tmin_c <= config_dates$TMPKILL,
           doy    >= config_dates$earliest_harvest) %>%
    slice(1)
  
  # Take earliest of maturity, frost, or hard cap
  harvest_candidates <- c(
    if (nrow(maturity_row) > 0) as.integer(maturity_row$doy[1]) else NA_integer_,
    if (nrow(frost_row)    > 0) as.integer(frost_row$doy[1])    else NA_integer_,
    as.integer(config_dates$latest_harvest)
  )
  
  as.integer(min(harvest_candidates, na.rm = TRUE))
}

# =============================================================================
# 5. COMPUTE MEAN DATES ACROSS ALL YEARS FOR ONE SITE
# =============================================================================

calc_site_dates <- function(weather_df, config_dates) {
  
  years <- unique(weather_df$year)
  
  annual_dates <- map_dfr(years, function(yr) {
    yr_data     <- filter(weather_df, year == yr)
    seed_doy    <- calc_seeding_doy(yr_data, config_dates)
    harvest_doy <- calc_harvest_doy(yr_data, seed_doy, config_dates)
    data.frame(year = yr, seed_doy = seed_doy, harvest_doy = harvest_doy)
  })
  
  list(
    mean_seed_doy    = as.integer(round(mean(annual_dates$seed_doy,    na.rm = TRUE))),
    mean_harvest_doy = as.integer(round(mean(annual_dates$harvest_doy, na.rm = TRUE))),
    sd_seed          = round(sd(annual_dates$seed_doy,    na.rm = TRUE), 1),
    sd_harvest       = round(sd(annual_dates$harvest_doy, na.rm = TRUE), 1),
    annual           = annual_dates
  )
}

# =============================================================================
# 6. UPDATE SCHEDULE FILE
# =============================================================================
# Updates seeding and harvest DOY in Crop_Opt.sch
# Preserves all other lines exactly

update_schedule_file <- function(sch_path, seed_doy, harvest_doy) {
  
  if (!file.exists(sch_path)) {
    warning(paste("Schedule not found:", sch_path))
    return(invisible(NULL))
  }
  
  lines <- readLines(sch_path)
  
  # Helper to replace DOY on a matched line
  replace_doy <- function(lines, pattern, new_doy) {
    idx <- grep(pattern, lines)
    if (length(idx) == 0) return(lines)
    # Replace the leading integer (DOY) on matched lines
    lines[idx] <- sub("^(\\s*1\\s+)[0-9]+", 
                      paste0("\\1", new_doy), 
                      lines[idx])
    return(lines)
  }
  
  # Seeding events — all on original DOY 105
  # KILL cultivation 5 days before seeding
  lines <- replace_doy(lines, "CULT KILL",  seed_doy - 5)
  lines <- replace_doy(lines, "CULT NDRIL", seed_doy)
  lines <- replace_doy(lines, "FERT",       seed_doy)
  lines <- replace_doy(lines, "CROP GLENT", seed_doy)
  lines <- replace_doy(lines, "PLTM",       seed_doy)
  
  # Harvest events
  lines <- replace_doy(lines, "LAST",   harvest_doy)
  lines <- replace_doy(lines, "HARV G", harvest_doy)
  
  writeLines(lines, sch_path)
}

# =============================================================================
# 7. MAIN LOOP — PROCESS ALL SITE FOLDERS
# =============================================================================

site_folders <- list.dirs(config_dates$sites_dir,
                          full.names = FALSE,
                          recursive  = FALSE)
site_folders <- site_folders[nchar(site_folders) > 0]

cat(sprintf("Processing %d site folders...\n", length(site_folders)))

# Requires zoo for rolling mean
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
library(zoo)

results_list <- vector("list", length(site_folders))

for (i in seq_along(site_folders)) {
  
  site_id  <- site_folders[i]
  site_dir <- file.path(config_dates$sites_dir, site_id)
  wth_path <- file.path(site_dir, config_dates$weather_name)
  sch_path <- file.path(site_dir, config_dates$schedule_name)
  
  # Read weather
  weather_df <- read_weather_wth(wth_path)
  if (is.null(weather_df)) next
  
  # Calculate dates
  dates <- calc_site_dates(weather_df, config_dates)
  
  # Update schedule file
  update_schedule_file(sch_path, dates$mean_seed_doy, dates$mean_harvest_doy)
  
  # Store results
  results_list[[i]] <- data.frame(
    LEGALLAND        = site_id,
    mean_seed_doy    = dates$mean_seed_doy,
    mean_harvest_doy = dates$mean_harvest_doy,
    sd_seed          = dates$sd_seed,
    sd_harvest       = dates$sd_harvest,
    season_length    = dates$mean_harvest_doy - dates$mean_seed_doy
  )
  
  if (i %% 25 == 0 || i == length(site_folders)) {
    cat(sprintf("  %d / %d complete\n", i, length(site_folders)))
  }
}

# Combine and save results
site_dates <- bind_rows(Filter(Negate(is.null), results_list))

write_csv(site_dates,
          file.path(config_dates$results_dir, "site_dates.csv"))

# =============================================================================
# 8. SUMMARY
# =============================================================================

cat("\n=== Site Date Calculation Complete ===\n")
cat(sprintf("Sites processed:       %d\n",   nrow(site_dates)))
cat(sprintf("Mean seeding DOY:      %.1f (SD: %.1f) — approx %s\n",
            mean(site_dates$mean_seed_doy),
            sd(site_dates$mean_seed_doy),
            format(as.Date(round(mean(site_dates$mean_seed_doy)) - 1,
                           origin = "2000-01-01"), "%B %d")))
cat(sprintf("Mean harvest DOY:      %.1f (SD: %.1f) — approx %s\n",
            mean(site_dates$mean_harvest_doy),
            sd(site_dates$mean_harvest_doy),
            format(as.Date(round(mean(site_dates$mean_harvest_doy)) - 1,
                           origin = "2000-01-01"), "%B %d")))
cat(sprintf("Mean season length:    %.1f days\n",
            mean(site_dates$season_length)))
cat(sprintf("\nSeeding DOY range:     %d to %d\n",
            min(site_dates$mean_seed_doy), max(site_dates$mean_seed_doy)))
cat(sprintf("Harvest DOY range:     %d to %d\n",
            min(site_dates$mean_harvest_doy), max(site_dates$mean_harvest_doy)))
cat(sprintf("\nSite dates saved to:   %s\n",
            file.path(config_dates$results_dir, "site_dates.csv")))
cat("Schedule files updated in all site folders\n")

# =============================================================================
# 9. VERIFY ONE SITE
# =============================================================================
# Quick check — read back one updated schedule to confirm changes

cat("\nVerifying first site schedule:\n")
first_sch <- file.path(config_dates$sites_dir,
                       site_folders[1], config_dates$schedule_name)
cat(readLines(first_sch), sep = "\n")