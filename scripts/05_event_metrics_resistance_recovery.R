library(terra)
library(data.table)

ndvi_z_path <- "PATH/TO/ndvi_monthly_zscore_1985_2025.nc"
events_dir <- "PATH/TO/EVENT_TABLES_DIRECTORY"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"

max_recovery_window <- 12

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

event_files <- c(
  HD = "hd_events_1985_2025.csv",
  HW = "hw_events_1985_2025.csv",
  HOT_ONLY = "hot_only_events_1985_2025.csv",
  DRY_ONLY = "dry_only_events_1985_2025.csv"
)

ndvi_z <- rast(ndvi_z_path)

dates <- time(ndvi_z)
if (!is.null(dates)) {
  dates <- as.Date(dates)
} else {
  dates <- as.Date(gsub("NDVI_Z_", "", names(ndvi_z)), format = "%Y_%m")
  time(ndvi_z) <- dates
}

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  ndvi_z <- crop(ndvi_z, aoi)
  ndvi_z <- mask(ndvi_z, aoi)
}

n_months <- nlyr(ndvi_z)
n_cols <- ncol(ndvi_z)

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

calc_recovery <- function(ts, end_idx, max_window = 12) {
  if (is.na(end_idx) || end_idx >= length(ts)) return(NA_real_)

  last_idx <- min(length(ts), end_idx + max_window)
  if ((last_idx - end_idx) < 2) return(NA_real_)

  post <- ts[(end_idx + 1):last_idx]
  ok <- which(post[-length(post)] >= 0 & post[-1] >= 0)

  if (!length(ok)) return(NA_real_)
  ok[1]
}

process_event_table <- function(events_path, event_label) {
  if (!file.exists(events_path)) return(NULL)

  ev <- fread(events_path)

  if (!nrow(ev)) return(NULL)

  ev <- ev[
    !is.na(cell) &
      !is.na(start_index) &
      !is.na(end_index) &
      cell >= 1 &
      cell <= ncell(ndvi_z) &
      start_index >= 1 &
      end_index <= n_months &
      start_index <= end_index
  ]

  if (!nrow(ev)) return(NULL)

  rc <- rowColFromCell(ndvi_z, ev$cell)

  ev[, `:=`(
    row = rc[, 1],
    col = rc[, 2],
    event_type = event_label
  )]

  bs <- blocks(ndvi_z)

  out <- vector("list", bs$n)
  kk <- 0L

  readStart(ndvi_z)
  on.exit(readStop(ndvi_z), add = TRUE)

  for (b in seq_len(bs$n)) {
    row_start <- bs$row[b]
    nrows_b <- bs$nrows[b]
    row_end <- row_start + nrows_b - 1L

    ev_b <- ev[row >= row_start & row <= row_end]

    if (!nrow(ev_b)) next

    mat <- readValues(ndvi_z, row = row_start, nrows = nrows_b, mat = TRUE)

    local_row <- ev_b$row - row_start + 1L
    local_idx <- (local_row - 1L) * n_cols + ev_b$col

    resistance <- rep(NA_real_, nrow(ev_b))
    recovery <- rep(NA_real_, nrow(ev_b))

    for (i in seq_len(nrow(ev_b))) {
      ts <- mat[local_idx[i], ]
      seg <- ts[ev_b$start_index[i]:ev_b$end_index[i]]

      if (!all(is.na(seg))) {
        resistance[i] <- min(seg, na.rm = TRUE)
      }

      recovery[i] <- calc_recovery(ts, ev_b$end_index[i], max_recovery_window)
    }

    ev_b[, `:=`(
      resistance = resistance,
      recovery_months = recovery,
      suppressed = fifelse(!is.na(resistance) & resistance < 0, 1L, 0L),
      recovered_12m = fifelse(!is.na(recovery_months), 1L, 0L)
    )]

    kk <- kk + 1L
    out[[kk]] <- ev_b[, .(
      cell,
      event_type,
      event_no,
      start_index,
      end_index,
      duration,
      start_date,
      end_date,
      resistance,
      recovery_months,
      suppressed,
      recovered_12m
    )]
  }

  out <- out[seq_len(kk)]

  if (!length(out)) return(NULL)

  event_metrics <- rbindlist(out)

  xy <- xyFromCell(ndvi_z, event_metrics$cell)

  event_metrics[, `:=`(
    x = xy[, 1],
    y = xy[, 2]
  )]

  pixel_summary <- event_metrics[, .(
    n_events = .N,
    mean_duration = mean(duration, na.rm = TRUE),
    median_duration = median(duration, na.rm = TRUE),
    mean_resistance = mean_na(resistance),
    median_resistance = median_na(resistance),
    mean_recovery = mean_na(recovery_months),
    median_recovery = median_na(recovery_months),
    suppression_prob = mean(suppressed, na.rm = TRUE),
    recovery_success_prob = mean(recovered_12m, na.rm = TRUE)
  ), by = .(cell, event_type)]

  xy2 <- xyFromCell(ndvi_z, pixel_summary$cell)

  pixel_summary[, `:=`(
    x = xy2[, 1],
    y = xy2[, 2]
  )]

  fwrite(
    event_metrics,
    file.path(output_dir, paste0(tolower(event_label), "_event_metrics.csv"))
  )

  fwrite(
    pixel_summary,
    file.path(output_dir, paste0(tolower(event_label), "_pixel_summary.csv"))
  )

  list(
    event_metrics = event_metrics,
    pixel_summary = pixel_summary
  )
}

results <- lapply(names(event_files), function(nm) {
  process_event_table(file.path(events_dir, event_files[[nm]]), nm)
})

names(results) <- names(event_files)

event_list <- lapply(results, function(x) if (!is.null(x)) x$event_metrics else NULL)
pixel_list <- lapply(results, function(x) if (!is.null(x)) x$pixel_summary else NULL)

event_list <- event_list[!vapply(event_list, is.null, logical(1))]
pixel_list <- pixel_list[!vapply(pixel_list, is.null, logical(1))]

if (length(event_list)) {
  all_event_metrics <- rbindlist(event_list, use.names = TRUE, fill = TRUE)
  fwrite(all_event_metrics, file.path(output_dir, "all_event_metrics.csv"))
}

if (length(pixel_list)) {
  all_pixel_summary <- rbindlist(pixel_list, use.names = TRUE, fill = TRUE)
  fwrite(all_pixel_summary, file.path(output_dir, "all_pixel_summary.csv"))
}

if (exists("all_pixel_summary")) {
  hd <- all_pixel_summary[event_type == "HD", .(cell, hd_resistance = mean_resistance)]
  hot_only <- all_pixel_summary[event_type == "HOT_ONLY", .(cell, hot_only_resistance = mean_resistance)]
  dry_only <- all_pixel_summary[event_type == "DRY_ONLY", .(cell, dry_only_resistance = mean_resistance)]

  amp <- merge(hd, hot_only, by = "cell", all = TRUE)
  amp <- merge(amp, dry_only, by = "cell", all = TRUE)

  amp[, strongest_single_driver := pmin(hot_only_resistance, dry_only_resistance, na.rm = TRUE)]
  amp[!is.finite(strongest_single_driver), strongest_single_driver := NA_real_]

  amp[, amplification := hd_resistance - strongest_single_driver]

  xy_amp <- xyFromCell(ndvi_z, amp$cell)

  amp[, `:=`(
    x = xy_amp[, 1],
    y = xy_amp[, 2]
  )]

  fwrite(amp, file.path(output_dir, "hd_amplification_pixel_summary.csv"))

  amp_rast <- rast(ndvi_z[[1]])
  values(amp_rast) <- NA_real_
  values(amp_rast)[amp$cell] <- amp$amplification
  names(amp_rast) <- "hd_amplification"

  writeRaster(
    amp_rast,
    file.path(output_dir, "hd_amplification.tif"),
    overwrite = TRUE
  )
}
