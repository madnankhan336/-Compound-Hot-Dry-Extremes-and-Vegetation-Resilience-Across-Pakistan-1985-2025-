library(terra)
library(data.table)

terraclimate_dir <- "PATH/TO/TERRACLIMATE_YEARLY_STACKS"
ndvi_z_path <- "PATH/TO/ndvi_monthly_zscore_1985_2025.nc"
elevation_path <- "PATH/TO/ELEVATION_RASTER.tif"
template_raster_path <- "PATH/TO/TEMPLATE_RASTER.tif"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tc_files <- list.files(
  terraclimate_dir,
  pattern = "^pakistan_terraclimate_monthly_\\d{4}\\.tif$",
  full.names = TRUE
)

stopifnot(length(tc_files) > 0)

tc <- rast(tc_files)
ndvi_z <- rast(ndvi_z_path)
elev <- rast(elevation_path)
template <- rast(template_raster_path)

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  tc <- crop(tc, aoi)
  tc <- mask(tc, aoi)
  ndvi_z <- crop(ndvi_z, aoi)
  ndvi_z <- mask(ndvi_z, aoi)
  elev <- crop(elev, aoi)
  elev <- mask(elev, aoi)
  template <- crop(template, aoi)
  template <- mask(template, aoi)
}

if (!compareGeom(tc, template, stopOnError = FALSE)) {
  tc <- resample(tc, template, method = "bilinear")
}

if (!compareGeom(ndvi_z, template, stopOnError = FALSE)) {
  ndvi_z <- resample(ndvi_z, template, method = "bilinear")
}

if (!compareGeom(elev, template, stopOnError = FALSE)) {
  elev <- resample(elev, template, method = "bilinear")
}

meta <- strcapture(
  "^([a-z]+)_(\\d{4})_(\\d{2})$",
  names(tc),
  proto = list(var = character(), year = integer(), month = integer())
)

stopifnot(!any(is.na(meta$var)), !any(is.na(meta$year)), !any(is.na(meta$month)))

meta$date <- as.Date(sprintf("%04d-%02d-01", meta$year, meta$month))

get_var_stack <- function(x, meta, varname) {
  idx <- which(meta$var == varname)
  idx <- idx[order(meta$date[idx])]
  out <- x[[idx]]
  names(out) <- paste0(varname, "_", format(meta$date[idx], "%Y_%m"))
  time(out) <- meta$date[idx]
  out
}

tmmx <- get_var_stack(tc, meta, "tmmx")
pr <- get_var_stack(tc, meta, "pr")
vpd <- get_var_stack(tc, meta, "vpd")
soil <- get_var_stack(tc, meta, "soil")
wb <- get_var_stack(tc, meta, "wb")

dates <- as.Date(time(tmmx))
if (is.null(time(ndvi_z))) time(ndvi_z) <- dates

month_id <- as.integer(format(dates, "%m"))

qfun <- function(p) {
  force(p)
  function(x) quantile(x, probs = p, na.rm = TRUE, names = FALSE)
}

t95 <- tapp(tmmx, index = month_id, fun = qfun(0.95))
t90 <- tapp(tmmx, index = month_id, fun = qfun(0.90))
pr10 <- tapp(pr, index = month_id, fun = qfun(0.10))
pr05 <- tapp(pr, index = month_id, fun = qfun(0.05))
wb10 <- tapp(wb, index = month_id, fun = qfun(0.10))
vpd90 <- tapp(vpd, index = month_id, fun = qfun(0.90))
soil10 <- tapp(soil, index = month_id, fun = qfun(0.10))

zone <- classify(
  elev,
  rcl = matrix(c(
    -Inf,  500, 1,
     500, 1500, 2,
    1500, 2500, 3,
    2500, 3500, 4,
    3500,  Inf, 5
  ), ncol = 3, byrow = TRUE),
  include.lowest = TRUE,
  right = FALSE
)

zone_labels <- data.table(
  zone = 1:5,
  elev_band = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)

valid_mask <- !is.na(zone)
valid_cells <- which(values(valid_mask, mat = FALSE) == 1)
n_valid_pixels <- length(valid_cells)

calc_recovery <- function(ts, end_idx, max_window = 12) {
  if (is.na(end_idx) || end_idx >= length(ts)) return(NA_real_)
  last_idx <- min(length(ts), end_idx + max_window)
  if ((last_idx - end_idx) < 2) return(NA_real_)
  post <- ts[(end_idx + 1):last_idx]
  ok <- which(post[-length(post)] >= 0 & post[-1] >= 0)
  if (!length(ok)) return(NA_real_)
  ok[1]
}

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
q25 <- function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q75 <- function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)

make_flag_stack <- function(hot_stack, dry_mode) {
  out <- vector("list", length(dates))

  for (i in seq_along(dates)) {
    m <- month_id[i]

    hot_i <- if (hot_stack == "t95") tmmx[[i]] > t95[[m]] else tmmx[[i]] > t90[[m]]

    dry_i <- switch(
      dry_mode,
      pr10 = pr[[i]] < pr10[[m]],
      pr05 = pr[[i]] < pr05[[m]],
      wb10 = wb[[i]] < wb10[[m]],
      vpd90_soil10 = (vpd[[i]] > vpd90[[m]]) & (soil[[i]] < soil10[[m]])
    )

    out[[i]] <- hot_i & dry_i
  }

  r <- rast(out)
  names(r) <- paste0("HD_", format(dates, "%Y_%m"))
  time(r) <- dates
  r
}

make_pixel_metrics <- function(flag_stack, scenario_name) {
  bs <- blocks(flag_stack)
  nc <- ncol(flag_stack)

  readStart(flag_stack)
  readStart(ndvi_z)
  readStart(zone)

  on.exit(readStop(flag_stack), add = TRUE)
  on.exit(readStop(ndvi_z), add = TRUE)
  on.exit(readStop(zone), add = TRUE)

  out_list <- vector("list", bs$n)
  kk <- 0L

  for (b in seq_len(bs$n)) {
    flag_mat <- readValues(flag_stack, row = bs$row[b], nrows = bs$nrows[b], mat = TRUE)
    ndvi_mat <- readValues(ndvi_z, row = bs$row[b], nrows = bs$nrows[b], mat = TRUE)
    zone_vals <- readValues(zone, row = bs$row[b], nrows = bs$nrows[b], mat = FALSE)

    cell_start <- cellFromRowCol(flag_stack, bs$row[b], 1)
    cell_end <- cellFromRowCol(flag_stack, bs$row[b] + bs$nrows[b] - 1, nc)
    cells <- seq(cell_start, cell_end)

    block_rows <- vector("list", nrow(flag_mat))
    k <- 0L

    for (i in seq_len(nrow(flag_mat))) {
      z <- zone_vals[i]
      if (is.na(z)) next

      x <- flag_mat[i, ]
      ts <- ndvi_mat[i, ]

      if (all(is.na(x))) next

      x <- ifelse(is.na(x), 0L, as.integer(x > 0))

      if (!any(x == 1L)) next

      r <- rle(x)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1L
      keep <- which(r$values == 1L)

      if (!length(keep)) next

      st <- starts[keep]
      en <- ends[keep]
      du <- r$lengths[keep]

      resistance <- rep(NA_real_, length(keep))
      recovery <- rep(NA_real_, length(keep))

      for (j in seq_along(keep)) {
        seg <- ts[st[j]:en[j]]
        if (!all(is.na(seg))) resistance[j] <- min(seg, na.rm = TRUE)
        recovery[j] <- calc_recovery(ts, en[j], 12)
      }

      k <- k + 1L
      block_rows[[k]] <- data.table(
        scenario = scenario_name,
        cell = cells[i],
        zone = z,
        n_events = length(keep),
        event_months = sum(x == 1L),
        mean_duration = mean(du, na.rm = TRUE),
        mean_resistance = safe_mean(resistance),
        mean_recovery = safe_mean(recovery),
        suppression_prob = mean(resistance < 0, na.rm = TRUE)
      )
    }

    if (k > 0) {
      kk <- kk + 1L
      out_list[[kk]] <- rbindlist(block_rows[seq_len(k)])
    }
  }

  out <- rbindlist(out_list[seq_len(kk)], use.names = TRUE)
  out <- merge(out, zone_labels, by = "zone", all.x = TRUE)
  out
}

make_raster_from_table <- function(template, dt, value_col, out_name) {
  r <- rast(template)
  values(r) <- NA_real_
  values(r)[dt$cell] <- dt[[value_col]]
  names(r) <- out_name
  r
}

scenario_defs <- list(
  list(name = "HD95_PR10_baseline", hot = "t95", dry = "pr10"),
  list(name = "HD90_PR10", hot = "t90", dry = "pr10"),
  list(name = "HD95_PR05", hot = "t95", dry = "pr05"),
  list(name = "HD95_WB10", hot = "t95", dry = "wb10"),
  list(name = "HD95_VPD90_SM10", hot = "t95", dry = "vpd90_soil10")
)

national_out <- list()
elev_out <- list()

for (sc in scenario_defs) {
  flag_stack <- make_flag_stack(sc$hot, sc$dry)

  writeRaster(
    flag_stack,
    file.path(output_dir, paste0(sc$name, "_flag.tif")),
    overwrite = TRUE
  )

  px <- make_pixel_metrics(flag_stack, sc$name)

  fwrite(
    px,
    file.path(output_dir, paste0(sc$name, "_pixel_metrics.csv"))
  )

  r_count <- make_raster_from_table(template, px, "n_events", paste0(sc$name, "_event_count"))
  r_res <- make_raster_from_table(template, px, "mean_resistance", paste0(sc$name, "_mean_resistance"))
  r_rec <- make_raster_from_table(template, px, "mean_recovery", paste0(sc$name, "_mean_recovery"))

  writeRaster(r_count, file.path(output_dir, paste0(sc$name, "_event_count.tif")), overwrite = TRUE)
  writeRaster(r_res, file.path(output_dir, paste0(sc$name, "_mean_resistance.tif")), overwrite = TRUE)
  writeRaster(r_rec, file.path(output_dir, paste0(sc$name, "_mean_recovery.tif")), overwrite = TRUE)

  national_out[[sc$name]] <- px[, .(
    scenario = sc$name,
    pixels_with_event_pct = 100 * .N / n_valid_pixels,
    total_events = sum(n_events, na.rm = TRUE),
    median_resistance = median(mean_resistance, na.rm = TRUE),
    median_recovery_months = median(mean_recovery, na.rm = TRUE)
  )]

  elev_out[[sc$name]] <- px[, .(
    scenario = sc$name,
    pixels_n = .N,
    hd_events_per_pixel_median = median(n_events, na.rm = TRUE),
    hd_events_per_pixel_q25 = q25(n_events),
    hd_events_per_pixel_q75 = q75(n_events),
    resistance_median = median(mean_resistance, na.rm = TRUE),
    resistance_q25 = q25(mean_resistance),
    resistance_q75 = q75(mean_resistance),
    recovery_median = median(mean_recovery, na.rm = TRUE),
    recovery_q25 = q25(mean_recovery),
    recovery_q75 = q75(mean_recovery)
  ), by = .(zone, elev_band)][order(zone)]
}

national_summary <- rbindlist(national_out, use.names = TRUE)
elev_summary <- rbindlist(elev_out, use.names = TRUE)

fwrite(
  national_summary,
  file.path(output_dir, "Table_2_sensitivity_national_summary.csv")
)

fwrite(
  elev_summary,
  file.path(output_dir, "Table_S9_sensitivity_by_elevation.csv")
)

baseline <- fread(file.path(output_dir, "HD95_PR10_baseline_pixel_metrics.csv"))

for (sc in scenario_defs[-1]) {
  d <- fread(file.path(output_dir, paste0(sc$name, "_pixel_metrics.csv")))
  d <- merge(
    d[, .(cell, mean_resistance, mean_recovery, n_events)],
    baseline[, .(cell, mean_resistance_base = mean_resistance, mean_recovery_base = mean_recovery, n_events_base = n_events)],
    by = "cell",
    all = TRUE
  )

  d[, diff_resistance := mean_resistance - mean_resistance_base]
  d[, diff_recovery := mean_recovery - mean_recovery_base]
  d[, diff_event_count := n_events - n_events_base]

  fwrite(
    d,
    file.path(output_dir, paste0(sc$name, "_minus_baseline_pixel_differences.csv"))
  )

  r_diff_res <- make_raster_from_table(template, d[!is.na(cell)], "diff_resistance", paste0(sc$name, "_minus_baseline_resistance"))
  r_diff_rec <- make_raster_from_table(template, d[!is.na(cell)], "diff_recovery", paste0(sc$name, "_minus_baseline_recovery"))
  r_diff_evt <- make_raster_from_table(template, d[!is.na(cell)], "diff_event_count", paste0(sc$name, "_minus_baseline_event_count"))

  writeRaster(
    r_diff_res,
    file.path(output_dir, paste0(sc$name, "_minus_baseline_resistance.tif")),
    overwrite = TRUE
  )

  writeRaster(
    r_diff_rec,
    file.path(output_dir, paste0(sc$name, "_minus_baseline_recovery.tif")),
    overwrite = TRUE
  )

  writeRaster(
    r_diff_evt,
    file.path(output_dir, paste0(sc$name, "_minus_baseline_event_count.tif")),
    overwrite = TRUE
  )
}
