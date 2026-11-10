# =============================================================================
# Set soil-zone-specific N rates (FERAMT(1)) and turn off automatic
# fertilization (AUFERT = 0) in every RM folder's fert.100.
#
# N rates: Saskatchewan Ministry of Agriculture (2019) Crop Planning Guide,
#   fertilizer assumptions by soil zone, converted lb N/ac -> g N m-2 (x 0.1121).
# Zones: RM centroid joined to Soil Zones of the Canadian Prairies.
#   Dark Gray / Gray (north) -> Black rate; durum has no Black rate -> Dark Brown.
#
# Run once per crop (set `crop` below). Originals are backed up as
# fert.100.bak_preNzone (only on the first run, never overwritten).
# =============================================================================

library(sf)
library(dplyr)

crop      <- "wheat"     # barley | wheat | durum | oats | canola | peas | green_lentil | red_lentil
crop_dir  <- "/home/preston/Projects/Daycent/Crop_Optimization/Models/Wheat/qrts"
dry_run   <- FALSE         # TRUE = report only; set FALSE to write files

rm_path   <- "/mnt/Elements/Saskatchewan/RM_Shapefile/RM_Centroids.shp"
zone_path <- "/mnt/Elements/Saskatchewan/Soil_Zones_of_the_Canadian_Prairies/Soil_Zones_of_the_Canadian_Prairies.shp"
out_csv   <- file.path(dirname(crop_dir), paste0("rm_zone_N_", crop, ".csv"))

# lb N/ac from the 2019 Crop Planning Guide (Brown, Dark Brown, Black),
# "Fertilizer" assumptions table for each soil zone. Matches the detrending
# reference period (yields adjusted to ~2019 levels).
lb <- list(
  barley       = c(85,  90,  99),    # feed barley
  wheat        = c(85,  98, 108),    # spring wheat
  durum        = c(94, 107, 107),    # no Black value; Dark Brown used
  oats         = c(61,  78,  99),
  canola       = c(102, 107, 111),
  peas         = c(8,   8,   10),    # edible yellow/green peas
  green_lentil = c(4,   4,    4),    # large green lentils
  red_lentil   = c(5,   5,    6)
)
rate <- setNames(lb[[crop]] * 0.1121, c("Brown", "Dark Brown", "Black"))

# --- RM centroids and soil zones ---------------------------------------------
zones <- st_read(zone_path, quiet = TRUE)
rms   <- st_read(rm_path,   quiet = TRUE) %>% st_transform(st_crs(zones))

# Find the attribute holding zone names (values mentioning Brown/Black)
cand <- names(zones)[sapply(zones, function(v)
  is.character(v) || is.factor(v))]
zfld <- cand[which.max(sapply(cand, function(f)
  sum(grepl("brown|black", zones[[f]], ignore.case = TRUE))))]
cat("Zone field:", zfld, "\nZone values:", paste(sort(unique(zones[[zfld]])), collapse = " | "), "\n")

norm_zone <- function(z) {
  z <- tolower(z)
  case_when(grepl("dark brown", z)            ~ "Dark Brown",
            grepl("brown", z)                 ~ "Brown",
            grepl("black|gray|grey", z)       ~ "Black",
            TRUE                              ~ NA_character_)
}

zones <- st_make_valid(zones)
j <- st_join(rms, zones[zfld], join = st_intersects, left = TRUE)
j <- j[!duplicated(j$RMNO), ]
# Centroids falling outside every polygon (e.g. far north): nearest polygon
miss <- is.na(j[[zfld]])
if (any(miss)) {
  nn <- st_nearest_feature(j[miss, ], zones)
  j[[zfld]][miss] <- zones[[zfld]][nn]
  cat("Nearest-polygon zone used for RMs:", paste(j$RMNO[miss], collapse = " "), "\n")
}

tab <- tibble(RMNO = as.character(as.integer(as.numeric(j$RMNO))),
              zone_raw = as.character(j[[zfld]]),
              zone = norm_zone(zone_raw)) %>%
  mutate(feramt = round(rate[zone], 3))

# --- Apply to RM folders -------------------------------------------------------
folders <- list.dirs(crop_dir, full.names = FALSE, recursive = FALSE)
folders <- folders[nchar(folders) > 0]
tab <- tab %>% filter(RMNO %in% folders)

unmatched <- setdiff(folders, tab$RMNO)
if (length(unmatched)) warning("No zone for folder(s): ", paste(unmatched, collapse = " "))
if (any(is.na(tab$zone))) {
  print(filter(tab, is.na(zone)))
  stop("Unrecognised zone names above — extend norm_zone().")
}

print(count(tab, zone, feramt))

if (!dry_run) {
  for (i in seq_len(nrow(tab))) {
    f <- file.path(crop_dir, tab$RMNO[i], "fert.100")
    if (!file.exists(f)) { warning("missing ", f); next }
    bak <- paste0(f, ".bak_preNzone")
    if (!file.exists(bak)) file.copy(f, bak)
    x <- readLines(f, warn = FALSE)
    i1 <- grep("FERAMT\\(1\\)", x); i2 <- grep("AUFERT", x)
    if (length(i1) != 1 || length(i2) != 1) { warning("unexpected layout: ", f); next }
    x[i1] <- sub("^\\s*[0-9.eE+-]+", sprintf("%.3f", tab$feramt[i]), x[i1])
    x[i2] <- sub("^\\s*[0-9.eE+-]+", "0", x[i2])
    writeLines(x, f, sep = "\n")
  }
  # Verify
  chk <- sapply(tab$RMNO, function(r) {
    x <- readLines(file.path(crop_dir, r, "fert.100"), warn = FALSE)
    c(feramt = as.numeric(sub("^\\s*([0-9.]+).*", "\\1", x[grep("FERAMT\\(1\\)", x)])),
      aufert = as.numeric(sub("^\\s*([0-9.]+).*", "\\1", x[grep("AUFERT", x)])))
  })
  cat("Written. AUFERT values:", paste(unique(chk["aufert", ]), collapse = ","),
      "| FERAMT(1) values:", paste(sort(unique(chk["feramt", ])), collapse = ","), "\n")
}

write.csv(tab, out_csv, row.names = FALSE)
cat("Zone table saved:", out_csv, "\n")
