library(terra)
library(data.table)

pixel_summary_path <- "PATH/TO/OUTPUT_DIRECTORY/all_pixel_summary.csv"
amplification_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_amplification_pixel_summary.csv"
elevation_path <- "PATH/TO/ELEVATION_RASTER.tif"
template_raster_path <- "PATH/TO/TEMPLATE_RASTER.tif"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

px <- fread(pixel_summary_path)
amp <- fread(amplification_path)

elev <- rast(elevation_path)
template <- rast(template_raster_path)

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  elev <- crop(elev, aoi)
  elev <- mask(elev, aoi)
  template <- crop(template, aoi)
  template <- mask(template, aoi)
}

if (!compareGeom(elev, template, stopOnError = FALSE)) {
  elev <- resample(elev, template, method = "bilinear")
}

elev_class <- classify(
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

levels_df <- data.frame(
  value = 1:5,
  band = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)

as_dt <- function(x) {
  setDT(as.data.frame(x))
}

extract_zone <- function(cells) {
  ext(cells(elev_class, cells))
}

zone_rast <- elev_class
zone_vals <- values(zone_rast, mat = FALSE)

px[, zone := zone_vals[cell]]
amp[, zone := zone_vals[cell]]

px <- px[!is.na(zone)]
amp <- amp[!is.na(zone)]

px[, band := factor(levels_df$band[match(zone, levels_df$value)], levels = levels_df$band)]
amp[, band := factor(levels_df$band[match(zone, levels_df$value)], levels = levels_df$band)]

q25 <- function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q75 <- function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)

fmt_iqr <- function(med, lo, hi, digits = 2) {
  paste0(
    format(round(med, digits), nsmall = digits, trim = TRUE),
    " [",
    format(round(lo, digits), nsmall = digits, trim = TRUE),
    "–",
    format(round(hi, digits), nsmall = digits, trim = TRUE),
    "]"
  )
}

summ_event_type <- function(dt, event_name) {
  d <- dt[event_type == event_name]

  if (!nrow(d)) return(NULL)

  out <- d[, .(
    pixels_n = .N,
    event_freq_median = median(n_events, na.rm = TRUE),
    event_freq_q25 = q25(n_events),
    event_freq_q75 = q75(n_events),
    duration_median = median(mean_duration, na.rm = TRUE),
    duration_q25 = q25(mean_duration),
    duration_q75 = q75(mean_duration),
    resistance_median = median(mean_resistance, na.rm = TRUE),
    resistance_q25 = q25(mean_resistance),
    resistance_q75 = q75(mean_resistance),
    recovery_median = median(mean_recovery, na.rm = TRUE),
    recovery_q25 = q25(mean_recovery),
    recovery_q75 = q75(mean_recovery)
  ), by = .(zone, band)]

  out[, `:=`(
    frequency_iqr = fmt_iqr(event_freq_median, event_freq_q25, event_freq_q75, 0),
    duration_iqr = fmt_iqr(duration_median, duration_q25, duration_q75, 2),
    resistance_iqr = fmt_iqr(resistance_median, resistance_q25, resistance_q75, 3),
    recovery_iqr = fmt_iqr(recovery_median, recovery_q25, recovery_q75, 2),
    event_type = event_name
  )]

  out[order(zone)]
}

hd_summary <- summ_event_type(px, "HD")
hw_summary <- summ_event_type(px, "HW")
hot_only_summary <- summ_event_type(px, "HOT_ONLY")
dry_only_summary <- summ_event_type(px, "DRY_ONLY")

if (!is.null(hd_summary)) fwrite(hd_summary, file.path(output_dir, "table_hd_by_elevation.csv"))
if (!is.null(hw_summary)) fwrite(hw_summary, file.path(output_dir, "table_hw_by_elevation.csv"))
if (!is.null(hot_only_summary)) fwrite(hot_only_summary, file.path(output_dir, "table_hot_only_by_elevation.csv"))
if (!is.null(dry_only_summary)) fwrite(dry_only_summary, file.path(output_dir, "table_dry_only_by_elevation.csv"))

if (nrow(amp)) {
  amp_summary <- amp[, .(
    n_pixels = .N,
    mean = mean(amplification, na.rm = TRUE),
    sd = sd(amplification, na.rm = TRUE),
    median = median(amplification, na.rm = TRUE),
    q25 = q25(amplification),
    q75 = q75(amplification)
  ), by = .(zone, band)][order(zone)]

  fwrite(amp_summary, file.path(output_dir, "table_hd_amplification_by_elevation.csv"))
}

make_table_s1 <- function(hd_summary) {
  if (is.null(hd_summary)) return(NULL)

  out <- hd_summary[, .(
    elevation_zone = zone,
    elevation_band_m_asl = as.character(band),
    median_hot_dry_frequency_events = frequency_iqr,
    median_mean_duration_months = duration_iqr,
    median_resistance_unitless = resistance_iqr,
    median_recovery_months = recovery_iqr
  )]

  fwrite(out, file.path(output_dir, "Table_S1_HD_exposure_response_by_elevation.csv"))
  out
}

make_table_s2 <- function(hw_summary) {
  if (is.null(hw_summary)) return(NULL)

  out <- hw_summary[, .(
    elevation_zone = zone,
    elevation_band_m_asl = as.character(band),
    median_hw_frequency_events = frequency_iqr,
    median_mean_duration_months = duration_iqr,
    median_resistance_zscore = resistance_iqr,
    median_recovery_months = recovery_iqr
  )]

  fwrite(out, file.path(output_dir, "Table_S2_HW_exposure_response_by_elevation.csv"))
  out
}

make_table_s3 <- function(amp_summary) {
  if (is.null(amp_summary)) return(NULL)

  out <- amp_summary[, .(
    elevation_zone = zone,
    elevation_band_m_asl = as.character(band),
    n_pixels,
    mean = round(mean, 3),
    sd = round(sd, 3),
    median = round(median, 3),
    q25 = round(q25, 3),
    q75 = round(q75, 3)
  )]

  fwrite(out, file.path(output_dir, "Table_S3_HD_amplification_by_elevation.csv"))
  out
}

tab_s1 <- make_table_s1(hd_summary)
tab_s2 <- make_table_s2(hw_summary)
tab_s3 <- if (exists("amp_summary")) make_table_s3(amp_summary) else NULL

if (!is.null(hd_summary) || !is.null(hw_summary)) {
  all_elev_summary <- rbindlist(
    list(hd_summary, hw_summary, hot_only_summary, dry_only_summary),
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(all_elev_summary, file.path(output_dir, "all_event_types_by_elevation.csv"))
}
