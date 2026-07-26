# =============================================================================
# DayCent Decomposition Parameters — Morris Elementary-Effects Screening
# Stage 1 of a two-stage GSA:
#   Stage 1 (this script): Morris screen of all 28 fix.100 decomposition
#     parameters -> rank by mu* (influence), flag interactions/nonlinearity via
#     sigma, and cull to the influential subset. Cost = r*(k+1) DayCent network
#     runs (e.g. 20*29 = 580) — runs the REAL model, no emulator.
#   Stage 2 (sobol_direct): direct Sobol/Jansen variance decomposition on the
#     surviving m parameters only, non-survivors fixed at default. Cost =
#     N*(m+2) runs — affordable once m ~ 10, and needs no emulator.
#
# Why Morris first: mu* rankings track Sobol total-order indices closely for a
# keep-vs-fix decision, at a fraction of the cost, and without the high-
# dimensional GP that the 28-parameter emulator Sobol has to lean on.
#
# Same DayCent interface as the calibration/Sobol scripts:
#   per-site folders, per-folder binary + fix.100; total SOC = sum of soilc.out
#   soil SOM pools only: som1c(2)+som2c(2)+som3c = cols 6,8,9 (g C/m2) -> Mg/ha
#   (/100); day-365 rows; network mean over sites. Surface pools (5,7) excluded.
# =============================================================================

library(sensitivity)   # morris(), soboljansen()
library(lhs)           # LHS for the Sobol stage
library(dplyr)
library(readr)
library(purrr)
library(parallel)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

BASE_DIR     <- "/home/preston/Projects/Daycent/SOC_Optimization/soc/qrts"
OUTPUT_DIR   <- "/home/preston/Projects/Daycent/SOC_Optimization/morris_output"
EXE_NAME     <- "DDcentEVI_rev491"     # per-folder binary; verify
SCHEDULE_ARG <- "Annual_Crop.sch"      # DayCent -s argument; verify

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Morris design
R_MORRIS   <- 20     # trajectories kept  -> R_MORRIS*(k+1) = 580 network runs
R_POOL     <- 150    # pool for Campolongo-optimized trajectory selection (spread)
LEVELS     <- 6      # grid levels per factor
GRID_JUMP  <- 3      # jump (levels/2 gives well-spread elementary effects)

# Screening decision
KEEP_REL   <- 0.10   # keep factors with mu*_rel >= 10% of the most influential
# (heuristic — the mu*-sigma plot is the real guide)

# Stage-2 Sobol (on survivors)
AUTO_SOBOL_STAGE <- FALSE   # TRUE -> chain straight into direct Sobol on survivors
N_SOBOL_DIRECT   <- 300     # Sobol base sample for the reduced set
NBOOT            <- 200

set.seed(42)

# =============================================================================
# 1. PARAMETER SPECIFICATION  (28 decomposition parameters; defaults from fix.100)
# -----------------------------------------------------------------------------
# type: rate / fraction / coef.  lo/hi NA -> auto band default*(1 +/- REL_WIDTH).
# *** Bounds define the space Morris explores — edit to real prior uncertainty
#     wherever known. The 10 Gurung-set parameters are pre-filled. ***
# =============================================================================

REL_WIDTH <- 0.50

PARAM_SPEC <- tibble::tribble(
  ~fix_name,    ~type,       ~default,    ~lo,    ~hi,
  "DEC1(1)",    "rate",      3.9000,      NA,     NA,
  "DEC1(2)",    "rate",      4.9000,      NA,     NA,
  "DEC2(1)",    "rate",      14.800,      NA,     NA,
  "DEC2(2)",    "rate",      18.500,      NA,     NA,
  "DEC3(1)",    "rate",      6.0000,      NA,     NA,
  "DEC3(2)",    "rate",      7.3000,      NA,     NA,
  "DEC4",       "rate",      0.004946,    0.001,  0.005,
  "DEC5(1)",    "rate",      0.10,        NA,     NA,
  "DEC5(2)",    "rate",      0.07925,     0.07,   0.25,
  "P1CO2A(1)",  "fraction",  0.6000,      NA,     NA,
  "P1CO2A(2)",  "fraction",  0.1700,      NA,     NA,
  "P1CO2B(1)",  "coef",      0.0000,      0.0,    0.10,   # default 0 -> absolute band; EDIT
  "P1CO2B(2)",  "coef",      0.6800,      NA,     NA,
  "P2CO2(1)",   "fraction",  0.5500,      NA,     NA,
  "P2CO2(2)",   "fraction",  0.56789,     0.50,   0.85,
  "P3CO2",      "fraction",  0.5500,      0.30,   0.75,
  "TEFF(1)",    "coef",      20.079654,   5.0,    30.0,
  "TEFF(2)",    "coef",      0.160946,    0.10,   0.30,
  "WEFF(1)",    "coef",      30.000,      20.0,   45.0,
  "WEFF(2)",    "coef",      14.686709,   6.0,    16.0,
  "PS1CO2(1)",  "fraction",  0.4500,      NA,     NA,
  "PS1CO2(2)",  "fraction",  0.5500,      NA,     NA,
  "PS1S3(1)",   "coef",      0.0030,      NA,     NA,
  "PS1S3(2)",   "coef",      0.054486,    0.02,   0.06,
  "PS2S3(1)",   "coef",      0.0040,      NA,     NA,
  "PS2S3(2)",   "coef",      0.0090,      NA,     NA,
  "PMCO2(1)",   "fraction",  0.5500,      NA,     NA,
  "PMCO2(2)",   "fraction",  0.669333,    0.30,   0.75
)

PARAM_SPEC$id <- make.names(PARAM_SPEC$fix_name)
ID_TO_FIX <- setNames(PARAM_SPEC$fix_name, PARAM_SPEC$id)
IDS       <- PARAM_SPEC$id
K         <- length(IDS)
DEFAULTS  <- setNames(PARAM_SPEC$default, IDS)

auto_band <- function(default, type, w = REL_WIDTH) {
  if (type == "fraction") {
    if (default <= 0) return(c(0, w))
    return(c(max(0, default * (1 - w)), min(1, default * (1 + w))))
  }
  if (type == "rate") return(c(max(0, default * (1 - w)), default * (1 + w)))
  if (abs(default) < 1e-8) return(c(0, 0.10))
  sort(c(default * (1 - w), default * (1 + w)))
}
bounds <- t(mapply(function(def, ty, lo, hi)
  if (!is.na(lo) && !is.na(hi)) c(lo, hi) else auto_band(def, ty),
  PARAM_SPEC$default, PARAM_SPEC$type, PARAM_SPEC$lo, PARAM_SPEC$hi))
LOWER <- setNames(bounds[, 1], IDS)
UPPER <- setNames(bounds[, 2], IDS)

n_cores_available <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (is.na(n) || n < 1) n <- max(1, detectCores() - 1)
  n
}
yr_col <- function(y) paste0("y", y)

# =============================================================================
# 2. CALIBRATION YEARS (from per-site observed.csv)
# =============================================================================

get_calib_years <- function(base_dir) {
  site_ids <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
  yrs <- purrr::map(site_ids, function(s) {
    f <- file.path(base_dir, s, "observed.csv")
    if (!file.exists(f)) return(NULL)
    as.integer(readr::read_csv(f, show_col_types = FALSE)$Year)
  })
  sort(unique(unlist(yrs)))
}

# =============================================================================
# 3. DAYCENT INTERFACE
# =============================================================================

write_decomp_params <- function(params_fix, site_dir) {
  fix_path <- file.path(site_dir, "fix.100")
  lines    <- readLines(fix_path)
  replace_one <- function(lines, fix_name, value) {
    escaped <- gsub("([()])", "\\\\\\1", fix_name)
    idx <- grep(paste0("^[^#]*", escaped, "(\\s|$|\t|#)"), lines, perl = TRUE)
    if (length(idx) == 0) { warning(sprintf("'%s' not found in fix.100", fix_name)); return(lines) }
    lines[idx[1]] <- sub("^(\\s*)([-0-9.eE+]+)",
                         sprintf("\\1%s", formatC(value, format = "f", digits = 7)),
                         lines[idx[1]])
    lines
  }
  for (fx in names(params_fix)) lines <- replace_one(lines, fx, params_fix[[fx]])
  writeLines(lines, fix_path)
}

run_daycent_site <- function(site_id, params_fix, base_dir, exe_name, calib_years) {
  site_dir <- file.path(base_dir, site_id)
  write_decomp_params(params_fix, site_dir)
  cmd <- sprintf("cd '%s' && rm -f output.bin && './%s' -s %s -n output > /dev/null 2>&1",
                 site_dir, exe_name, SCHEDULE_ARG)
  if (system(cmd) != 0) return(NULL)
  soilc_path <- file.path(site_dir, "soilc.out")
  if (!file.exists(soilc_path)) return(NULL)
  soilc <- tryCatch(read.table(soilc_path, header = TRUE), error = function(e) NULL)
  if (is.null(soilc) || ncol(soilc) < 9) return(NULL)
  data.frame(
    year         = as.integer(soilc$time),
    dayofyr      = soilc$dayofyr,
    soc_020_Mgha = rowSums(soilc[, c(6, 8, 9)]) / 100   # soil SOM only: som1c(2)+som2c(2)+som3c
  ) %>%
    filter(dayofyr == 365, year %in% calib_years) %>%
    select(year, soc_020_Mgha)
}

run_network_mean <- function(params_id, matched_sites, base_dir, exe_name, calib_years) {
  params_fix <- setNames(as.numeric(params_id), ID_TO_FIX[names(params_id)])
  res <- mclapply(matched_sites,
                  function(s) run_daycent_site(s, params_fix, base_dir, exe_name, calib_years),
                  mc.cores = n_cores_available())
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) return(NULL)
  bind_rows(res) %>%
    group_by(year) %>%
    summarise(mean_sim = mean(soc_020_Mgha, na.rm = TRUE), .groups = "drop")
}

# Evaluate a FULL 28-column (id-named) design matrix -> per-year SOC + dSOC.
# One network run per row; parallel over sites within each run.
evaluate_design <- function(Xmat, matched_sites, calib_years, label = "design") {
  n <- nrow(Xmat)
  yr_cols <- yr_col(calib_years)
  out <- matrix(NA_real_, n, length(calib_years), dimnames = list(NULL, yr_cols))
  n_ok <- 0L
  for (i in seq_len(n)) {
    if (i %% 10 == 0) message(sprintf("  %s eval %d / %d (%d ok)", label, i, n, n_ok))
    p  <- setNames(as.numeric(Xmat[i, ]), colnames(Xmat))
    sm <- tryCatch(run_network_mean(p, matched_sites, BASE_DIR, EXE_NAME, calib_years),
                   error = function(e) NULL)
    if (!is.null(sm)) { out[i, yr_col(sm$year)] <- sm$mean_sim; n_ok <- n_ok + 1L }
  }
  if (n_ok == 0)
    stop("Every DayCent run failed (0 valid outputs). Likely wrong SCHEDULE_ARG / EXE_NAME, ",
         "or no day-365 output at the calibration years. Run preflight_check() to diagnose.")
  if (n_ok < 0.5 * n)
    warning(sprintf("%d of %d runs produced no output — Morris indices may be unreliable.", n - n_ok, n))
  df <- as.data.frame(out)
  if (length(calib_years) >= 2)
    df$dSOC <- df[[yr_col(max(calib_years))]] - df[[yr_col(min(calib_years))]]
  df
}

# =============================================================================
# 3b. PREFLIGHT — one default run to catch config errors before the full design
# =============================================================================

preflight_check <- function(matched_sites, calib_years) {
  message("=== Preflight: single default-parameter run on one site ===")
  s        <- matched_sites[1]
  site_dir <- file.path(BASE_DIR, s)
  sch      <- list.files(site_dir, pattern = "\\.sch$", ignore.case = TRUE)
  message(sprintf("  Site folder : %s", site_dir))
  message(sprintf("  .sch present : %s", if (length(sch)) paste(sch, collapse = ", ") else "NONE"))
  
  if (!file.exists(file.path(site_dir, EXE_NAME)))
    stop(sprintf("Binary '%s' not found in the site folder — check EXE_NAME.", EXE_NAME))
  
  sched_ok <- file.exists(file.path(site_dir, SCHEDULE_ARG)) ||
    file.exists(file.path(site_dir, paste0(SCHEDULE_ARG, ".sch")))
  if (!sched_ok)
    stop(sprintf("SCHEDULE_ARG '%s' not found in the site folder. Present .sch files: %s",
                 SCHEDULE_ARG, if (length(sch)) paste(sch, collapse = ", ") else "none"))
  
  write_decomp_params(setNames(as.numeric(DEFAULTS), ID_TO_FIX[IDS]), site_dir)
  cmd <- sprintf("cd '%s' && rm -f output.bin && './%s' -s %s -n output 2>&1", site_dir, EXE_NAME, SCHEDULE_ARG)
  log <- suppressWarnings(system(cmd, intern = TRUE))
  st  <- attr(log, "status")
  if (!is.null(st) && st != 0) {
    message(paste(tail(log, 15), collapse = "\n"))
    stop("DayCent exited non-zero on the test run — see the tail of its output above.")
  }
  
  sc <- file.path(site_dir, "soilc.out")
  if (!file.exists(sc)) stop("Test run produced no soilc.out — check schedule/output settings.")
  soilc <- read.table(sc, header = TRUE)
  hit   <- soilc[soilc$dayofyr == 365 & as.integer(soilc$time) %in% calib_years, , drop = FALSE]
  if (nrow(hit) == 0)
    stop(sprintf("soilc.out has no day-365 rows at years %s — check the schedule end year / output config.",
                 paste(calib_years, collapse = ", ")))
  message(sprintf("  OK: %d of %d calibration years present at day 365; sample SOC = %.1f Mg/ha",
                  length(unique(as.integer(hit$time))), length(calib_years),
                  rowSums(hit[1, c(6, 8, 9)]) / 100))
  invisible(TRUE)
}

# =============================================================================
# 4. STAGE 1 — MORRIS SCREEN
# =============================================================================

run_morris_screen <- function(matched_sites, calib_years, output_dir) {
  n_cores <- n_cores_available()
  n_runs  <- R_MORRIS * (K + 1)
  message(sprintf("=== Morris screen: %d network runs (r=%d, k=%d) ===", n_runs, R_MORRIS, K))
  message(sprintf("Rough wall time: ~%.1f h (assumes ~18 s/site, %d sites, %d cores)",
                  n_runs * length(matched_sites) * 18 / 3600 / n_cores,
                  length(matched_sites), n_cores))
  
  mo_path <- file.path(output_dir, "morris_design.rds")
  y_path  <- file.path(output_dir, "morris_outputs.csv")
  
  if (file.exists(mo_path) && file.exists(y_path)) {
    mo <- readRDS(mo_path)
    Y  <- read_csv(y_path, show_col_types = FALSE)
    if (all(vapply(Y, function(col) all(is.na(col)), logical(1)))) {
      message("Cached outputs are all NA (a previous run produced no DayCent output) — recomputing.")
      Y <- evaluate_design(mo$X, matched_sites, calib_years, label = "morris")
      write_csv(Y, y_path)
    } else {
      message("Loading cached Morris design + outputs.")
    }
  } else {
    mo <- morris(model = NULL, factors = IDS, r = c(R_MORRIS, R_POOL),
                 design = list(type = "oat", levels = LEVELS, grid.jump = GRID_JUMP),
                 binf = as.numeric(LOWER[IDS]), bsup = as.numeric(UPPER[IDS]))
    saveRDS(mo, mo_path)
    Y <- evaluate_design(mo$X, matched_sites, calib_years, label = "morris")
    write_csv(Y, y_path)
  }
  
  resp_cols <- names(Y)
  per_resp  <- list()
  for (rc in resp_cols) {
    mo_r <- mo
    mo_r <- tell(mo_r, Y[[rc]])
    ee   <- mo_r$ee                                  # r x k elementary effects
    mu_star <- apply(ee, 2, function(x) mean(abs(x), na.rm = TRUE))
    mu      <- apply(ee, 2, mean, na.rm = TRUE)
    sigma   <- apply(ee, 2, sd,   na.rm = TRUE)
    
    d <- data.frame(parameter = ID_TO_FIX[IDS], id = IDS, response = rc,
                    mu_star = mu_star, mu = mu, sigma = sigma,
                    mu_star_rel = mu_star / max(mu_star, na.rm = TRUE),
                    row.names = NULL)
    per_resp[[rc]] <- d[order(-d$mu_star), ]
    
    # classic mu*-sigma plot
    mx <- suppressWarnings(max(mu_star, na.rm = TRUE))
    if (is.finite(mx) && mx > 0) {
      png(file.path(output_dir, sprintf("morris_%s.png", rc)), width = 800, height = 700)
      plot(mu_star, sigma, pch = 19, col = "steelblue",
           xlab = expression(mu * "*  (mean |elementary effect|)"), ylab = expression(sigma),
           main = sprintf("Morris — %s", rc), xlim = c(0, mx * 1.15))
      text(mu_star, sigma, labels = ID_TO_FIX[IDS], pos = 4, cex = 0.7)
      abline(v = KEEP_REL * mx, col = "red", lty = 2)
      dev.off()
    } else {
      warning(sprintf("mu* all zero/NA for %s — skipping plot.", rc))
    }
  }
  
  results <- bind_rows(per_resp)
  write_csv(results, file.path(output_dir, "morris_per_response.csv"))
  
  # aggregate: max relative influence across responses -> keep/fix
  agg <- results %>%
    group_by(parameter, id) %>%
    summarise(max_mu_star_rel = max(mu_star_rel, na.rm = TRUE),
              max_mu_star     = max(mu_star, na.rm = TRUE),
              which_resp      = response[which.max(mu_star_rel)],
              max_sigma       = max(sigma, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(decision = ifelse(max_mu_star_rel >= KEEP_REL, "keep", "fix")) %>%
    arrange(desc(max_mu_star_rel))
  
  write_csv(agg, file.path(output_dir, "morris_ranking.csv"))
  message("\n=== Morris ranking (max relative mu* across metrics) ===")
  print(as.data.frame(agg %>% select(parameter, max_mu_star_rel, max_mu_star,
                                     which_resp, decision)), row.names = FALSE)
  
  survivor_ids <- agg$id[agg$decision == "keep"]
  list(results = results, ranking = agg, survivor_ids = survivor_ids)
}

# =============================================================================
# 5. STAGE 2 — DIRECT SOBOL ON SURVIVORS (non-survivors fixed at default)
# =============================================================================

sobol_direct <- function(survivor_ids, n_base = N_SOBOL_DIRECT, output_dir = OUTPUT_DIR) {
  m <- length(survivor_ids)
  if (m < 2) { message("Need >= 2 survivors for Sobol."); return(NULL) }
  matched_sites <- list.dirs(BASE_DIR, full.names = FALSE, recursive = FALSE)
  calib_years   <- get_calib_years(BASE_DIR)
  message(sprintf("=== Direct Sobol on %d survivors: %d*(%d+2) = %d network runs ===",
                  m, n_base, m, n_base * (m + 2)))
  
  lo <- LOWER[survivor_ids]; hi <- UPPER[survivor_ids]
  scale01 <- function(M) sapply(seq_len(m), function(j) lo[j] + M[, j] * (hi[j] - lo[j]))
  X1 <- scale01(randomLHS(n_base, m)); colnames(X1) <- survivor_ids
  X2 <- scale01(randomLHS(n_base, m)); colnames(X2) <- survivor_ids
  so <- soboljansen(model = NULL, X1 = as.data.frame(X1), X2 = as.data.frame(X2), nboot = NBOOT)
  
  # expand survivor design to full 28 columns, non-survivors at default
  full <- matrix(rep(DEFAULTS, each = nrow(so$X)), nrow = nrow(so$X),
                 dimnames = list(NULL, IDS))
  full[, survivor_ids] <- as.matrix(so$X)
  
  Y <- evaluate_design(full, matched_sites, calib_years, label = "sobol")
  
  out <- bind_rows(lapply(names(Y), function(rc) {
    s2 <- so; s2 <- tell(s2, Y[[rc]])
    data.frame(parameter = ID_TO_FIX[survivor_ids], response = rc,
               first_S = s2$S$original, total_T = s2$T$original,
               interaction = pmax(0, s2$T$original - s2$S$original), row.names = NULL)
  }))
  out <- out %>% arrange(response, desc(total_T))
  write_csv(out, file.path(output_dir, "sobol_direct_survivors.csv"))
  message("\n=== Direct Sobol on survivors ===")
  print(as.data.frame(out), row.names = FALSE)
  out
}

# =============================================================================
# 6. MAIN
# =============================================================================

main <- function() {
  message("=== DayCent Decomposition — Morris Screening ==="); message(Sys.time())
  matched_sites <- list.dirs(BASE_DIR, full.names = FALSE, recursive = FALSE)
  calib_years   <- get_calib_years(BASE_DIR)
  message(sprintf("Sites: %d | Years: %s | Factors: %d",
                  length(matched_sites), paste(calib_years, collapse = ", "), K))
  
  preflight_check(matched_sites, calib_years)
  screen  <- run_morris_screen(matched_sites, calib_years, OUTPUT_DIR)
  keepers <- screen$survivor_ids
  message(sprintf("\nSurvivors kept for optimization (%d): %s",
                  length(keepers), paste(ID_TO_FIX[keepers], collapse = ", ")))
  message(sprintf("Fixed at default (%d): %s",
                  K - length(keepers),
                  paste(setdiff(ID_TO_FIX[IDS], ID_TO_FIX[keepers]), collapse = ", ")))
  
  if (AUTO_SOBOL_STAGE) {
    sob <- sobol_direct(keepers, N_SOBOL_DIRECT, OUTPUT_DIR)
    invisible(list(morris = screen, sobol = sob))
  } else {
    message("\nReview the mu*-sigma plots and morris_ranking.csv, adjust the keep set")
    message("if needed, then run Stage 2, e.g.:")
    message("  sobol_direct(make.names(c('DEC5(2)','PMCO2(2)','P2CO2(2)','TEFF(1)')))")
    invisible(list(morris = screen))
  }
}

# =============================================================================
# RUN
# =============================================================================
result <- main()

# Then (after reviewing Morris output) run the trimmed Sobol on your chosen set:
sob <- sobol_direct(make.names(c("TEFF(1)","PMCO2(2)","P1CO2A(1)","P1CO2B(2)","DEC5(2)",
                                 "P1CO2A(2)","PMCO2(1)","PS1CO2(2)","TEFF(2)")),
                    n_base = 300)   # 300*(9+2) = 3,300 network runs

