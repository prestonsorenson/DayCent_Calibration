# =============================================================================
# ERA5-Land Precipitation Bias Correction
# Saskatchewan — using ECCC Climate Normals 1981-2010
# =============================================================================
# Method: Delta correction (ratio-based)
#   For each site and month:
#     correction_factor = station_normal_precip / era5_mean_precip
#     corrected_precip  = era5_precip * correction_factor
#
# This preserves ERA5's interannual variability while correcting
# the systematic monthly bias against observed station normals.
#
# Re-run safety: a pristine copy of every weather.wth is made ONCE into
# backup_dir. All reads (ERA5 means, daily correction, before/after plot)
# come from that backup; writes only ever go to sites_dir. This makes N
# runs identical to 1 run. NOTE: this treats whatever is in sites_dir at
# first run as pristine — if sites_dir may already hold corrected files
# from an earlier script version, delete the backup dir before running.
#
# Requirements:
#   install.packages("weathercan")
#   install.packages("sf")
# =============================================================================

library(tidyverse)
library(weathercan)
library(sf)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

config_bc <- list(
  sites_dir   = "/home/preston/Projects/Daycent/SOC_Optimization/soc/qrts",
  results_dir = "/home/preston/Projects/Daycent/SOC_Optimization/soc/results",
  
  # Site locations — needs x (LEGALLAND), lon, lat
  site_locs_path = "/home/preston/Projects/Daycent/SOC_Optimization/soc/pscb_locations.csv",
  
  # Search radius for station matching (km)
  station_search_radius = 100,
  
  # Normals period
  normals_years = "1981-2010"
)
# NOTE: dropped the unused `min_station_years` field — it was never referenced.
# ECCC normals carry per-element completeness codes (A/B/C/D); wiring a real
# minimum-years filter would mean parsing that code, which this script does
# not currently do. Left out rather than implying a guard that isn't there.

dir.create(config_bc$results_dir, recursive = TRUE, showWarnings = FALSE)

# Backup dir kept test-tagged so it can't collide with a non-test run.
backup_dir <- file.path(config_bc$sites_dir, "..", "qrts_test_pre_bias_correction")

# -- Shared weather.wth reader --------------------------------------------------
# Files are the extended DayCent format (10 cols: the 7 below + solrad, RH,
# windspeed). Read everything, name only the leading known columns, and leave
# any extras as V8, V9, ... so they survive downstream and get written back.
WTH_COLS <- c("day", "month", "year", "doy", "tmax_c", "tmin_c", "precip_cm")

read_wth <- function(path) {
  df <- read.table(path)
  names(df)[seq_along(WTH_COLS)] <- WTH_COLS
  df
}

# =============================================================================
# 2. LOAD SITE LOCATIONS
# =============================================================================

site_locs <- read_csv(config_bc$site_locs_path, show_col_types = FALSE)
stopifnot(all(c("x", "lon", "lat") %in% names(site_locs)))
cat(sprintf("Loaded %d site locations\n", nrow(site_locs)))

# =============================================================================
# 2b. BACK UP ORIGINAL WEATHER FILES (ONCE) — pristine read source
# =============================================================================

all_site_folders <- list.dirs(config_bc$sites_dir,
                              full.names = FALSE, recursive = FALSE)
all_site_folders <- all_site_folders[nchar(all_site_folders) > 0]
# Don't treat the backup dir itself as a site if it lands under sites_dir
all_site_folders <- setdiff(all_site_folders, basename(backup_dir))

if (!dir.exists(backup_dir)) {
  cat("Creating one-time backup of original weather files...\n")
  for (site_id in all_site_folders) {
    src <- file.path(config_bc$sites_dir, site_id, "weather.wth")
    if (!file.exists(src)) next
    backup_site <- file.path(backup_dir, site_id)
    dir.create(backup_site, recursive = TRUE, showWarnings = FALSE)
    file.copy(src, file.path(backup_site, "weather.wth"), overwrite = FALSE)
  }
  cat("Backup complete\n")
} else {
  cat("Backup dir exists — reading pristine weather from backup.\n")
}

# =============================================================================
# 3. DOWNLOAD (OR LOAD CACHED) ECCC CLIMATE NORMALS FOR SASKATCHEWAN
# =============================================================================

index_url     <- "https://dd.weather.gc.ca/today/climate/observations/normals/csv/1981-2010/SK"
normals_cache <- file.path(config_bc$results_dir, "station_normals.csv")

# Parse one station CSV file — returns monthly precip normals + metadata
parse_normals_csv <- function(url) {
  tryCatch({
    lines <- readLines(url, warn = FALSE, encoding = "latin1")
    line5 <- lines[5]
    
    station_name <- gsub('"', "", regmatches(line5,
                                             regexpr('^"([^"]+)"', line5))[[1]])
    
    climate_id <- gsub('"', "", regmatches(line5,
                                           regexpr('"([0-9]{7})"', line5))[[1]])
    
    dms_to_dd <- function(dms_str) {
      nums <- as.numeric(regmatches(dms_str,
                                    gregexpr("[0-9]+[.]?[0-9]*", dms_str))[[1]])
      if (length(nums) < 2) return(NA_real_)
      dd <- nums[1] + nums[2]/60 + ifelse(length(nums) >= 3, nums[3]/3600, 0)
      if (grepl("W|S", dms_str)) dd <- -dd
      return(dd)
    }
    
    parts <- strsplit(line5, '","')[[1]]
    parts <- gsub('^\"|\"$', "", parts)
    lat <- dms_to_dd(parts[3])
    lon <- dms_to_dd(parts[4])
    
    # Match the total-precipitation row specifically. "Precipitation (mm)" won't
    # match "Rainfall (mm)"/"Snowfall (cm)"; the "(mm)" immediately after the
    # word also excludes "Days with Precipitation >= 1 mm". If a file ever has
    # >1 match, this takes the first — worth a spot-check on a real file.
    precip_idx <- grep("Precipitation [(]mm[)]", lines)
    if (length(precip_idx) == 0) return(NULL)
    
    precip_line <- lines[precip_idx[1]]
    precip_vals <- strsplit(precip_line, ",")[[1]]
    precip_vals <- gsub('"', "", trimws(precip_vals))
    monthly_vals <- suppressWarnings(as.numeric(precip_vals[2:13]))
    
    data.frame(
      climate_id       = climate_id,
      station_name     = station_name,
      lat              = lat,
      lon              = lon,
      month            = 1:12,
      normal_precip_mm = monthly_vals,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

if (file.exists(normals_cache)) {
  cat("Loading cached station normals from:", normals_cache, "\n")
  # climate_id as character to preserve any leading structure / avoid coercion
  station_normals <- read_csv(normals_cache,
                              col_types = cols(climate_id = col_character(),
                                               .default = col_guess()))
} else {
  cat("Downloading Saskatchewan climate normals from MSC Datamart...\n")
  index_page <- readLines(index_url, warn = FALSE)
  
  csv_files <- regmatches(
    index_page,
    gregexpr('href="(climate_normals_SK_[^"]+\\.csv)"', index_page)
  )
  csv_files <- gsub('href="|"', "", unlist(csv_files))
  cat(sprintf("Found %d station CSV files\n", length(csv_files)))
  
  normals_list <- vector("list", length(csv_files))
  for (i in seq_along(csv_files)) {
    url <- paste0(index_url, "/", csv_files[i])
    normals_list[[i]] <- parse_normals_csv(url)
    if (i %% 50 == 0) cat(sprintf("  %d / %d files processed\n", i, length(csv_files)))
  }
  
  # Surface silent drops: which station files failed to parse / download
  failed <- csv_files[vapply(normals_list, is.null, logical(1))]
  if (length(failed) > 0) {
    cat(sprintf("  WARNING: %d/%d station file(s) returned no data:\n",
                length(failed), length(csv_files)))
    cat(paste0("    ", failed, "\n"))
  }
  
  station_normals <- bind_rows(Filter(Negate(is.null), normals_list)) %>%
    filter(!is.na(normal_precip_mm)) %>%
    mutate(normal_precip_cm = normal_precip_mm / 10) %>%
    distinct(climate_id, month, .keep_all = TRUE)
  
  write_csv(station_normals, normals_cache)
}

cat(sprintf("Normals: %d station-months across %d stations\n",
            nrow(station_normals), n_distinct(station_normals$climate_id)))

# =============================================================================
# 4. MATCH EACH SITE TO NEAREST STATION WITH NORMALS
# =============================================================================

cat("Matching sites to nearest stations...\n")

sites_sf <- st_as_sf(site_locs, coords = c("lon", "lat"), crs = 4326)

stations_sf <- station_normals %>%
  distinct(climate_id, station_name, lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

nearest_idx     <- st_nearest_feature(sites_sf, stations_sf)
nearest_dist_km <- as.numeric(st_distance(
  sites_sf, stations_sf[nearest_idx, ], by_element = TRUE
)) / 1000

site_station_match <- site_locs %>%
  mutate(
    climate_id   = stations_sf$climate_id[nearest_idx],
    station_name = stations_sf$station_name[nearest_idx],
    dist_km      = round(nearest_dist_km, 1)
  )

cat(sprintf("Mean distance to nearest station: %.1f km\n",
            mean(site_station_match$dist_km)))
cat(sprintf("Max distance to nearest station:  %.1f km\n",
            max(site_station_match$dist_km)))

site_station_match <- site_station_match %>%
  mutate(station_ok = dist_km <= config_bc$station_search_radius)

cat(sprintf("Sites within %d km of a station: %d / %d\n",
            config_bc$station_search_radius,
            sum(site_station_match$station_ok),
            nrow(site_station_match)))

write_csv(site_station_match,
          file.path(config_bc$results_dir, "site_station_match.csv"))

# =============================================================================
# 5. COMPUTE ERA5 MONTHLY MEANS PER SITE  (read pristine from backup)
# =============================================================================

cat("Computing ERA5 monthly means per site...\n")

era5_monthly_means <- map_dfr(all_site_folders, function(site_id) {
  wth_path <- file.path(backup_dir, site_id, "weather.wth")
  if (!file.exists(wth_path)) return(NULL)
  
  df <- read_wth(wth_path)
  
  # Monthly totals (not daily means) to match station normal units
  df %>%
    group_by(year, month) %>%
    summarise(monthly_total = sum(precip_cm), .groups = "drop") %>%
    group_by(month) %>%
    summarise(era5_mean_precip_cm = mean(monthly_total), .groups = "drop") %>%
    mutate(x = site_id)
})

cat(sprintf("Computed ERA5 monthly means for %d sites\n",
            n_distinct(era5_monthly_means$x)))

# =============================================================================
# 6. COMPUTE CORRECTION FACTORS
# =============================================================================

correction_factors <- site_station_match %>%
  filter(station_ok) %>%
  select(x, climate_id) %>%
  left_join(
    station_normals %>% select(climate_id, month, normal_precip_cm),
    by = "climate_id"
  ) %>%
  left_join(era5_monthly_means, by = c("x", "month")) %>%
  mutate(
    correction_factor = case_when(
      era5_mean_precip_cm > 0 ~ normal_precip_cm / era5_mean_precip_cm,
      TRUE                    ~ 1.0
    ),
    correction_factor = pmax(0.3, pmin(2.0, correction_factor))
  )

cat("\nCorrection factor summary:\n")
print(summary(correction_factors$correction_factor))

write_csv(correction_factors,
          file.path(config_bc$results_dir, "correction_factors.csv"))

# =============================================================================
# 7. APPLY BIAS CORRECTION  (read pristine from backup, write to sites_dir)
# =============================================================================

cat("\nApplying bias correction to weather files...\n")

sites_to_correct <- unique(correction_factors$x)
sites_to_correct <- sites_to_correct[
  file.exists(file.path(backup_dir, sites_to_correct, "weather.wth"))
]
cat(sprintf("Sites with weather files to correct: %d\n", length(sites_to_correct)))
n_corrected <- 0
n_skipped   <- 0

for (site_id in sites_to_correct) {
  
  src_wth <- file.path(backup_dir, site_id, "weather.wth")          # pristine
  dst_wth <- file.path(config_bc$sites_dir, site_id, "weather.wth") # target
  
  site_cf <- correction_factors %>%
    filter(x == site_id) %>%
    select(month, correction_factor)
  
  if (nrow(site_cf) < 12) {          # skip if missing months
    n_skipped <- n_skipped + 1
    next
  }
  
  # Preserve extra columns (solrad/RH/windspeed): join appends
  # correction_factor at the end, so dropping only that column keeps the
  # original column order and every extra intact.
  df <- read_wth(src_wth) %>%
    left_join(site_cf, by = "month") %>%
    mutate(precip_cm = pmax(0, round(precip_cm * correction_factor, 4))) %>%
    select(-correction_factor)
  
  write.table(df, dst_wth,
              row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = " ")
  
  n_corrected <- n_corrected + 1
}

cat(sprintf("Bias correction applied to %d weather files (%d skipped for <12 months)\n",
            n_corrected, n_skipped))

# =============================================================================
# 8. VERIFY CORRECTION
# =============================================================================

cat("\nVerifying correction — annual precip totals after correction:\n")

verify_sites <- head(sites_to_correct, 5)   # head() is safe if <5 sites

for (site_id in verify_sites) {
  wth_path <- file.path(config_bc$sites_dir, site_id, "weather.wth")
  df <- read_wth(wth_path)
  mean_annual <- df %>%
    group_by(year) %>%
    summarise(annual = sum(precip_cm), .groups = "drop") %>%
    summarise(mean_annual = mean(annual)) %>%
    pull(mean_annual)
  cat(sprintf("  %s: %.1f cm/year\n", site_id, mean_annual))
}

# =============================================================================
# 9. SUMMARY PLOT — before vs after
# =============================================================================

cat("\nGenerating before/after comparison...\n")

annual_mean <- function(dir_root, site_id, label) {
  wth_path <- file.path(dir_root, site_id, "weather.wth")
  if (!file.exists(wth_path)) return(NULL)
  read_wth(wth_path) %>%
    group_by(year) %>%
    summarise(annual_cm = sum(precip_cm), .groups = "drop") %>%
    summarise(mean_annual_cm = mean(annual_cm)) %>%
    mutate(x = site_id, type = label)
}

corrected_means <- map_dfr(sites_to_correct, ~annual_mean(config_bc$sites_dir, .x, "corrected"))
original_means  <- map_dfr(sites_to_correct, ~annual_mean(backup_dir,          .x, "original"))

comparison <- bind_rows(original_means, corrected_means)

p <- ggplot(comparison, aes(x = mean_annual_cm, fill = type)) +
  geom_histogram(alpha = 0.6, bins = 20, position = "identity") +
  scale_fill_manual(values = c("original" = "steelblue", "corrected" = "tomato")) +
  labs(
    title = "ERA5 Precipitation Before vs After Bias Correction",
    subtitle = "Saskatchewan sites — mean annual precipitation",
    x = "Mean Annual Precipitation (cm)",
    y = "Number of Sites",
    fill = NULL
  ) +
  theme_bw()

ggsave(file.path(config_bc$results_dir, "bias_correction_comparison.png"),
       p, width = 8, height = 5, dpi = 150)

cat(sprintf("\nComparison plot saved to: %s\n",
            file.path(config_bc$results_dir, "bias_correction_comparison.png")))
cat("Bias correction complete\n")
