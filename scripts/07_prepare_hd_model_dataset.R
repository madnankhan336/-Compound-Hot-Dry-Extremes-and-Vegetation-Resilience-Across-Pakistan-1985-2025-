library(terra)
library(data.table)

terraclimate_dir <- "PATH/TO/TERRACLIMATE_YEARLY_STACKS"
ndvi_z_path <- "PATH/TO/ndvi_monthly_zscore_1985_2025.nc"
event_metrics_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_event_metrics.csv"
template_raster_path <- "PATH/TO/TEMPLATE_RASTER.tif"
elevation_path <- "PATH/TO/ELEVATION_RASTER.tif"
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
template <- rast(template_raster_path)
elev <- rast(elevation_path)
events <- fread(event_metrics_path)

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  tc <- crop(tc, aoi)
  tc <- mask(tc, aoi)
  ndvi_z <- crop(ndvi_z, aoi)
  ndvi_z <- mask(ndvi_z, aoi)
  template <- crop(template, aoi)
  template <- mask(template, aoi)
  elev <- crop(elev, aoi)
  elev <- mask(elev, aoi)
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
pet <- get_var_stack(tc, meta, "pet")
aet <- get_var_stack(tc, meta, "aet")
wb <- get_var_stack(tc, meta, "wb")

dates <- as.Date(time(tmmx))

if (is.null(time(ndvi_z))) {
  time(ndvi_z) <- dates
}

events <- events[
  !is.na(cell) &
    !is.na(start_index) &
    !is.na(end_index) &
    !is.na(resistance) &
    start_index >= 1 &
    end_index <= length(dates) &
    start_index <= end_index
]

events <- unique(events)

u_cells <- sort(unique(events$cell))

extract_matrix <- function(r, cells) {
  x <- extract(r, cells, cells = TRUE)
  x <- as.data.table(x)
  setnames(x, names(x)[1], "ID")
  x[, ID := NULL]
  as.matrix(x)
}

tmmx_mat <- extract_matrix(tmmx, u_cells)
pr_mat   <- extract_matrix(pr, u_cells)
vpd_mat  <- extract_matrix(vpd, u_cells)
soil_mat <- extract_matrix(soil, u_cells)
pet_mat  <- extract_matrix(pet, u_cells)
aet_mat  <- extract_matrix(aet, u_cells)
wb_mat   <- extract_matrix(wb, u_cells)
ndvi_mat <- extract_matrix(ndvi_z, u_cells)

xy <- xyFromCell(template, u_cells)
elev_vals <- values(elev, mat = FALSE)[u_cells]

cell_ref <- data.table(
  cell = u_cells,
  row_id = seq_along(u_cells),
  x = xy[, 1],
  y = xy[, 2],
  elev_m = elev_vals
)

setkey(cell_ref, cell)
setkey(events, cell)

events <- cell_ref[events]

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

n <- nrow(events)

tmax_mean <- rep(NA_real_, n)
pr_mean   <- rep(NA_real_, n)
vpd_mean  <- rep(NA_real_, n)
soil_min  <- rep(NA_real_, n)
pet_mean  <- rep(NA_real_, n)
aet_mean  <- rep(NA_real_, n)
wb_mean   <- rep(NA_real_, n)
ndvi_pre  <- rep(NA_real_, n)

for (i in seq_len(n)) {
  r <- events$row_id[i]
  s <- events$start_index[i]
  e <- events$end_index[i]

  idx <- s:e

  tmax_mean[i] <- safe_mean(tmmx_mat[r, idx])
  pr_mean[i]   <- safe_mean(pr_mat[r, idx])
  vpd_mean[i]  <- safe_mean(vpd_mat[r, idx])
  soil_min[i]  <- safe_min(soil_mat[r, idx])
  pet_mean[i]  <- safe_mean(pet_mat[r, idx])
  aet_mean[i]  <- safe_mean(aet_mat[r, idx])
  wb_mean[i]   <- safe_mean(wb_mat[r, idx])

  pre_start <- max(1, s - 6)
  pre_end <- s - 1

  if (pre_end >= pre_start) {
    ndvi_pre[i] <- safe_mean(ndvi_mat[r, pre_start:pre_end])
  }
}

events[, `:=`(
  tmax_mean = tmax_mean,
  pr_mean = pr_mean,
  vpd_mean = vpd_mean,
  soil_min = soil_min,
  pet_mean = pet_mean,
  aet_mean = aet_mean,
  wb_mean = wb_mean,
  ndvi_pre = ndvi_pre
)]

events[, event_start_date := as.Date(start_date)]
events[, event_end_date := as.Date(end_date)]
events[, duration := as.integer(duration)]

events[, elev_band := fifelse(
  elev_m < 500, "<500",
  fifelse(
    elev_m < 1500, "500–1500",
    fifelse(
      elev_m < 2500, "1500–2500",
      fifelse(
        elev_m < 3500, "2500–3500",
        "≥3500"
      )
    )
  )
)]

events[, pixel_id := as.integer(factor(cell))]
events[, event_id := .I]

model_dt <- events[, .(
  event_id,
  pixel_id,
  cell,
  event_start_date,
  event_end_date,
  start_index,
  end_index,
  duration,
  resistance,
  recovery_months,
  suppressed,
  recovered_12m,
  ndvi_pre,
  tmax_mean,
  pr_mean,
  vpd_mean,
  soil_min,
  pet_mean,
  aet_mean,
  wb_mean,
  elev_m,
  elev_band,
  x,
  y
)]

fwrite(
  model_dt,
  file.path(output_dir, "hd_model_dataset_raw.csv")
)

model_dt_complete <- model_dt[
  !is.na(resistance) &
    !is.na(ndvi_pre) &
    !is.na(tmax_mean) &
    !is.na(pr_mean) &
    !is.na(vpd_mean) &
    !is.na(soil_min) &
    !is.na(wb_mean) &
    !is.na(elev_m) &
    !is.na(x) &
    !is.na(y)
]

fwrite(
  model_dt_complete,
  file.path(output_dir, "hd_model_dataset_complete.csv")
)
