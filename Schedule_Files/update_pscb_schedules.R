# =============================================================================
# update_pscb_schedules.R
# Update DayCent .sch files for each PSCB site folder based on the
# site-specific 4-crop rotation from pscb_site_rotations.csv.
#
# SCHEDULE FILE STRUCTURE (Harv_Test.sch):
#   - Header lines (years, site file, flags, initial conditions)
#   - Block header (repeats, output settings, weather file)
#   - 4-year rotation body: each year has management events
#   - Terminal line: -999 -999 X
#
# ROTATION LOGIC:
#   - 4 crops ranked by frequency from SCIC yield records
#   - Each crop gets one year slot in the rotation (years 1-4)
#   - Crop-specific fertilizer, planting DOY, harvest DOY applied
#   - Legumes (peas, lentils) get no N fertilizer
#   - All sites: no-till from 1996 (CULT NDRIL = no-till drill)
#
# FOLDER STRUCTURE:
#   base_dir/
#   ├── NE-03-19-09-2/
#   │   ├── Harv_Test.sch   <- updated by this script
#   │   └── ...
#   ├── NE-34-30-04-2/
#   └── ...
#
# MATCHING:
#   Rotation CSV uses Location column (e.g. "NE-03-19-09-2")
#   which matches folder names exactly.
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE_DIR     <- "/home/preston/Projects/Daycent/SOC_Optimization/soc/qrts"
ROTATION_CSV <- "/home/preston/Projects/Daycent/SOC_Optimization/soc/pscb_site_rotations.csv"
SCH_FILENAME <- "Annual_Crop.sch"
WTH_FILENAME <- "weather.wth"

DRY_RUN <- FALSE

# =============================================================================
# 1. CROP PARAMETERS
#
# Phenology parameters from calc_site_dates_*.R scripts — consistent with
# crop.100 calibration (BASETEMP=0 throughout).
# GDD method: Miller et al. (2001) with base/ceiling temperatures.
# Seeding: first DOY in window where 5-day rolling mean Tmean >= TMPGERM.
# Harvest: earliest of cumulative GDD >= DDBASE, killing frost, or latest cap.
# CULT KILL placed 5 days before seeding DOY.
# =============================================================================

# DayCent crop codes
CROP_CODES <- c(
  "Wheat - HRSW" = "SW0",
  "Wheat -Durum"  = "DURUM",
  "Can/Rapeseed"  = "CAN",
  "IP Canola"     = "CAN",
  "Barley"        = "BAR3",
  "Oats"          = "OAT1",
  "Field Peas"    = "PEA",
  "Lentils-Red"   = "RLENT",
  "Lentils-LGrn"  = "GLENT",
  "Lentils-Oth"   = "GLENT"
)

LEGUMES   <- c("PEA", "GLENT", "RLENT")
FERT_CODE <- list(
  "SW0"   = "N10U", "DURUM" = "N10U", "CAN"   = "N10U",
  "BAR3"  = "N10U", "OAT1"  = "N10U",
  "PEA"   = NULL,   "GLENT" = NULL,   "RLENT" = NULL
)

# Phenology config per crop code
CROP_CONFIG <- list(
  "SW0" = list(
    TMPGERM = 2.0, TMPKILL = -6.0, DDBASE = 1700, BASETEMP = 0, BASETEMP2 = 30,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 200, latest_harvest = 280
  ),
  "DURUM" = list(  # no dedicated script — use SW0 parameters as closest match
    TMPGERM = 2.0, TMPKILL = -6.0, DDBASE = 1700, BASETEMP = 0, BASETEMP2 = 30,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 200, latest_harvest = 280
  ),
  "CAN" = list(
    TMPGERM = 2.0, TMPKILL = -4.0, DDBASE = 1900, BASETEMP = 0, BASETEMP2 = 29,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 220, latest_harvest = 290
  ),
  "BAR3" = list(
    TMPGERM = 1.0, TMPKILL = -5.0, DDBASE = 1500, BASETEMP = 0, BASETEMP2 = 30,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 200, latest_harvest = 280
  ),
  "OAT1" = list(
    TMPGERM = 1.0, TMPKILL = -5.0, DDBASE = 1600, BASETEMP = 0, BASETEMP2 = 30,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 200, latest_harvest = 280
  ),
  "PEA" = list(
    TMPGERM = 2.0, TMPKILL = -4.0, DDBASE = 1100, BASETEMP = 0, BASETEMP2 = 25,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 195, latest_harvest = 275
  ),
  "GLENT" = list(
    TMPGERM = 2.0, TMPKILL = -5.0, DDBASE = 1450, BASETEMP = 0, BASETEMP2 = 27,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 210, latest_harvest = 280
  ),
  "RLENT" = list(
    TMPGERM = 2.0, TMPKILL = -5.0, DDBASE = 1350, BASETEMP = 0, BASETEMP2 = 27,
    earliest_seed = 105, latest_seed = 166, earliest_harvest = 195, latest_harvest = 270
  )
)

# =============================================================================
# 2. WEATHER FILE READER
# =============================================================================

read_weather_wth <- function(wth_path) {
  if (!file.exists(wth_path)) return(NULL)
  tryCatch({
    raw <- read.table(wth_path, header = FALSE, fill = TRUE)
    raw %>%
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
  }, error = function(e) NULL)
}

# =============================================================================
# 3. SEEDING DATE — 5-day rolling mean Tmean >= TMPGERM
# =============================================================================

calc_seeding_doy <- function(year_data, cfg) {
  candidates <- year_data %>%
    filter(doy >= cfg$earliest_seed, doy <= cfg$latest_seed) %>%
    arrange(doy) %>%
    mutate(tmean_5day = zoo::rollmean(tmean, k = 5, fill = NA, align = "right"))

  seed_row <- candidates %>%
    filter(!is.na(tmean_5day), tmean_5day >= cfg$TMPGERM) %>%
    slice(1)

  if (nrow(seed_row) == 0) return(as.integer(cfg$latest_seed))
  as.integer(seed_row$doy)
}

# =============================================================================
# 4. HARVEST DATE — GDD accumulation with base/ceiling, frost, hard cap
# =============================================================================

calc_harvest_doy <- function(year_data, seed_doy, cfg) {
  if (is.na(seed_doy)) return(as.integer(cfg$latest_harvest))

  growing <- year_data %>%
    filter(doy >= seed_doy) %>%
    arrange(doy) %>%
    mutate(
      tmax_capped  = pmin(tmax_c, cfg$BASETEMP2),
      tmin_floored = pmax(tmin_c, cfg$BASETEMP),
      gdd_daily    = pmax(0, (tmax_capped + tmin_floored) / 2 - cfg$BASETEMP),
      gdd_cumul    = cumsum(gdd_daily)
    )

  maturity_row <- growing %>%
    filter(gdd_cumul >= cfg$DDBASE, doy >= cfg$earliest_harvest) %>%
    slice(1)

  frost_row <- growing %>%
    filter(tmin_c <= cfg$TMPKILL, doy >= cfg$earliest_harvest) %>%
    slice(1)

  harvest_candidates <- c(
    if (nrow(maturity_row) > 0) as.integer(maturity_row$doy[1]) else NA_integer_,
    if (nrow(frost_row)    > 0) as.integer(frost_row$doy[1])    else NA_integer_,
    as.integer(cfg$latest_harvest)
  )
  as.integer(min(harvest_candidates, na.rm = TRUE))
}

# =============================================================================
# 5. SITE MEAN DATES ACROSS ALL YEARS
# =============================================================================

calc_site_dates <- function(weather_df, cfg) {
  years <- unique(weather_df$year)
  annual <- purrr::map_dfr(years, function(yr) {
    yr_data  <- filter(weather_df, year == yr)
    seed_doy <- calc_seeding_doy(yr_data, cfg)
    harv_doy <- calc_harvest_doy(yr_data, seed_doy, cfg)
    data.frame(year = yr, seed_doy = seed_doy, harvest_doy = harv_doy)
  })
  list(
    mean_seed    = as.integer(round(mean(annual$seed_doy,    na.rm = TRUE))),
    mean_harvest = as.integer(round(mean(annual$harvest_doy, na.rm = TRUE)))
  )
}

# =============================================================================
# 6. BUILD EVENT BLOCK FOR ONE CROP IN ONE YEAR SLOT
# =============================================================================

build_crop_block <- function(year_slot, crop_code, seed_doy, harvest_doy) {
  fert <- FERT_CODE[[crop_code]]
  lines <- c(
    sprintf(" %d %d CULT KILL",  year_slot, seed_doy - 5),
    sprintf(" %d %d CULT NDRIL", year_slot, seed_doy),
    if (!is.null(fert))
      sprintf(" %d %d FERT %s",  year_slot, seed_doy, fert),
    sprintf(" %d %d CROP %s",    year_slot, seed_doy, crop_code),
    sprintf(" %d %d PLTM   ",    year_slot, seed_doy),
    sprintf(" %d %d LAST",       year_slot, harvest_doy),
    sprintf(" %d %d HARV G10S",  year_slot, harvest_doy)
  )
  lines[!sapply(lines, is.null)]
}

# =============================================================================
# 7. BUILD COMPLETE ROTATION BLOCK (all 4 year slots)
# =============================================================================

build_rotation_block <- function(crops, site_doys) {
  # crops:     character vector of 4 DayCent crop codes
  # site_doys: named list, one entry per crop code with seed/harvest DOYs
  unlist(lapply(seq_along(crops), function(i) {
    doys <- site_doys[[crops[i]]]
    build_crop_block(
      year_slot    = i,
      crop_code    = crops[i],
      seed_doy     = doys$mean_seed,
      harvest_doy  = doys$mean_harvest
    )
  }))
}

# =============================================================================
# 8. READ SCHEDULE FILE — EXTRACT HEADER + BLOCK METADATA
# =============================================================================

parse_schedule <- function(sch_path) {

  lines <- readLines(sch_path, warn = FALSE)

  # Find the "Year Month Option" header line — rotation body starts after it
  block_header_idx <- which(trimws(lines) == "Year Month Option")
  if (length(block_header_idx) == 0) {
    stop(sprintf("Cannot find 'Year Month Option' line in %s", sch_path))
  }

  # Everything before and including "Year Month Option" is the file header
  file_header <- lines[seq_len(block_header_idx)]

  # After "Year Month Option": block metadata lines then rotation events
  # Block metadata: block#, last year, repeats, output start, output month,
  #                 output interval, weather choice, weather file = 8 lines
  block_meta_start <- block_header_idx + 1
  block_meta_end   <- block_header_idx + 8
  block_meta       <- lines[block_meta_start:block_meta_end]

  list(
    file_header      = file_header,
    block_meta       = block_meta,
    block_header_idx = block_header_idx
  )
}

# =============================================================================
# 9. WRITE UPDATED SCHEDULE FILE
# =============================================================================

write_schedule <- function(sch_path, parsed, rotation_events) {

  new_lines <- c(
    parsed$file_header,
    parsed$block_meta,
    rotation_events,
    "-999 -999 X"
  )

  writeLines(new_lines, sch_path)
}

# =============================================================================
# 10. MAIN LOOP — process each site folder
# =============================================================================

main <- function() {

  if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
  library(zoo)

  message("=== PSCB Schedule File Updater ===")
  if (DRY_RUN) message("DRY RUN — no files will be written\n")

  rotations <- read_csv(ROTATION_CSV, show_col_types = FALSE) %>%
    select(Location, crop_1, crop_2, crop_3, crop_4) %>%
    distinct(Location, .keep_all = TRUE)

  message(sprintf("Loaded %d unique site rotations", nrow(rotations)))

  site_folders <- list.dirs(BASE_DIR, full.names = TRUE, recursive = FALSE)
  message(sprintf("Found %d site folders in %s\n", length(site_folders), BASE_DIR))

  n_updated  <- 0
  n_skipped  <- 0
  n_no_match <- 0
  n_error    <- 0

  for (folder in site_folders) {

    site_name <- basename(folder)
    sch_path  <- file.path(folder, SCH_FILENAME)
    wth_path  <- file.path(folder, WTH_FILENAME)

    if (!file.exists(sch_path)) {
      message(sprintf("  SKIP (no .sch): %s", site_name))
      n_skipped <- n_skipped + 1
      next
    }

    rot <- rotations %>% filter(Location == site_name)
    if (nrow(rot) == 0) {
      message(sprintf("  NO MATCH in rotation CSV: %s", site_name))
      n_no_match <- n_no_match + 1
      next
    }

    # Map and deduplicate crop codes
    raw_crops  <- c(rot$crop_1, rot$crop_2, rot$crop_3, rot$crop_4)
    crop_codes <- CROP_CODES[raw_crops]
    unmapped   <- raw_crops[is.na(crop_codes)]
    if (length(unmapped) > 0) {
      message(sprintf("  ERROR — unmapped crops at %s: %s",
                      site_name, paste(unmapped, collapse = ", ")))
      n_error <- n_error + 1
      next
    }
    crop_codes <- unname(crop_codes)
    seen_codes <- c()
    crop_codes <- sapply(crop_codes, function(code) {
      if (code %in% seen_codes) {
        message(sprintf("    NOTE — duplicate '%s' replaced with OAT1", code))
        code <- "OAT1"
      }
      seen_codes <<- c(seen_codes, code)
      code
    })

    # Calculate site-specific seeding/harvest DOYs from weather file
    unique_crops <- unique(crop_codes)
    site_doys    <- list()

    wth <- if (file.exists(wth_path)) read_weather_wth(wth_path) else NULL

    for (cc in unique_crops) {
      cfg <- CROP_CONFIG[[cc]]
      if (!is.null(wth) && !is.null(cfg)) {
        site_doys[[cc]] <- calc_site_dates(wth, cfg)
      } else {
        # Fallback to window midpoints if weather unavailable
        site_doys[[cc]] <- list(mean_seed = 115, mean_harvest = 258)
        if (is.null(wth))
          message(sprintf("  WARN: no weather file for %s — using fallback DOY", site_name))
      }
    }

    # Parse schedule, build rotation, write
    parsed <- tryCatch(
      parse_schedule(sch_path),
      error = function(e) {
        message(sprintf("  ERROR parsing %s: %s", site_name, e$message))
        return(NULL)
      }
    )
    if (is.null(parsed)) { n_error <- n_error + 1; next }

    rotation_events <- build_rotation_block(crop_codes, site_doys)

    if (DRY_RUN) {
      summary <- paste(sapply(crop_codes, function(cc)
        sprintf("%s(s%d/h%d)", cc, site_doys[[cc]]$mean_seed,
                site_doys[[cc]]$mean_harvest)), collapse = " / ")
      message(sprintf("  WOULD UPDATE: %s  [%s]", site_name, summary))
    } else {
      tryCatch({
        write_schedule(sch_path, parsed, rotation_events)
        summary <- paste(sapply(crop_codes, function(cc)
          sprintf("%s(s%d/h%d)", cc, site_doys[[cc]]$mean_seed,
                  site_doys[[cc]]$mean_harvest)), collapse = " / ")
        message(sprintf("  UPDATED: %s  [%s]", site_name, summary))
        n_updated <- n_updated + 1
      }, error = function(e) {
        message(sprintf("  ERROR writing %s: %s", site_name, e$message))
        n_error <<- n_error + 1
      })
    }
  }

  message("\n=== Summary ===")
  message(sprintf("  Updated:       %d", n_updated))
  message(sprintf("  Skipped:       %d  (no .sch file)", n_skipped))
  message(sprintf("  No CSV match:  %d", n_no_match))
  message(sprintf("  Errors:        %d", n_error))
  if (DRY_RUN) message("\n  Set DRY_RUN <- FALSE to apply changes.")
}

# =============================================================================
# HELPER: verify one updated schedule file
# =============================================================================

verify_schedule <- function(site_name,
                             base_dir = BASE_DIR,
                             sch_file = SCH_FILENAME) {
  path <- file.path(base_dir, site_name, sch_file)
  if (!file.exists(path)) { cat("File not found:", path, "\n"); return(invisible(NULL)) }
  cat(readLines(path), sep = "\n")
}

# =============================================================================
# RUN
# =============================================================================

main()

# After running, spot-check a site:
verify_schedule("NE-03-19-09-2")

