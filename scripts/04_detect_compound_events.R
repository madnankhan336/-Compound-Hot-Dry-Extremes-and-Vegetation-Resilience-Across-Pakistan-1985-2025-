library(terra)
library(data.table)

terraclimate_dir <- "PATH/TO/TERRACLIMATE_YEARLY_STACKS"
ndvi_z_path <- "PATH/TO/ndvi_monthly_zscore_1985_2025.nc"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"

hot_prob <- 0.95
dry_prob <- 0.10
wet_prob <- 0.90

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tc_files <- list.files(
  terraclimate_dir,
  pattern = "^pakistan_terraclimate_monthly_\\d{4}\\.tif$",
  full.names = TRUE
)

stopifnot(length(tc_files) > 0)

tc <- rast(tc_files)

info <- strcapture(
  "^([a-z]+)_(\\d{4})_(\\d{2})$",
  names(tc),
  proto = list(var = character(), year = integer(), month = integer())
)

stopifnot(!any(is.na(info$year)), !any(is.na(info$month)))

info$date <- as.Date(sprintf("%04d-%02d-01", info$year, info$month))

get_var_stack <- function(x, meta, varname) {
  idx <- which(meta$var == varname)
  idx <- idx[order(meta$date[idx])]
  out <- x[[idx]]
  names(out) <- paste0(varname, "_", format(meta$date[idx], "%Y_%m"))
  time(out) <- meta$date[idx]
  out
}

tmmx <- get_var_stack(tc, info, "tmmx")
pr <- get_var_stack(tc, info, "pr")

dates <- time(tmmx)
stopifnot(length(dates) == nlyr(tmmx), nlyr(tmmx) == nlyr(pr))

ndvi_z <- rast(ndvi_z_path)

if (!is.null(time(ndvi_z))) {
  ndvi_dates <- as.Date(time(ndvi_z))
} else {
  ndvi_dates <- dates
  time(ndvi_z) <- dates
}

stopifnot(length(ndvi_dates) == length(dates))
stopifnot(all(ndvi_dates == dates))

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  tmmx <- crop(tmmx, aoi)
  tmmx <- mask(tmmx, aoi)
  pr <- crop(pr, aoi)
  pr <- mask(pr, aoi)
  ndvi_z <- crop(ndvi_z, aoi)
  ndvi_z <- mask(ndvi_z, aoi)
}

if (!compareGeom(tmmx, pr, stopOnError = FALSE)) {
  pr <- resample(pr, tmmx, method = "bilinear")
}

if (!compareGeom(tmmx, ndvi_z, stopOnError = FALSE)) {
  ndvi_z <- resample(ndvi_z, tmmx, method = "bilinear")
}

month_id <- as.integer(format(dates, "%m"))

t_hot <- tapp(tmmx, index = month_id, fun = function(x) quantile(x, probs = hot_prob, na.rm = TRUE))
p_dry <- tapp(pr, index = month_id, fun = function(x) quantile(x, probs = dry_prob, na.rm = TRUE))
p_wet <- tapp(pr, index = month_id, fun = function(x) quantile(x, probs = wet_prob, na.rm = TRUE))

names(t_hot) <- paste0("t_hot_", sprintf("%02d", 1:12))
names(p_dry) <- paste0("p_dry_", sprintf("%02d", 1:12))
names(p_wet) <- paste0("p_wet_", sprintf("%02d", 1:12))

hot_layers <- vector("list", length(dates))
dry_layers <- vector("list", length(dates))
wet_layers <- vector("list", length(dates))
hd_layers <- vector("list", length(dates))
hw_layers <- vector("list", length(dates))
hot_only_layers <- vector("list", length(dates))
dry_only_layers <- vector("list", length(dates))

for (i in seq_along(dates)) {
  m <- month_id[i]

  hot_i <- tmmx[[i]] > t_hot[[m]]
  dry_i <- pr[[i]] < p_dry[[m]]
  wet_i <- pr[[i]] > p_wet[[m]]

  hd_i <- hot_i & dry_i
  hw_i <- hot_i & wet_i
  hot_only_i <- hot_i & !dry_i & !wet_i
  dry_only_i <- dry_i & !hot_i

  hot_layers[[i]] <- hot_i
  dry_layers[[i]] <- dry_i
  wet_layers[[i]] <- wet_i
  hd_layers[[i]] <- hd_i
  hw_layers[[i]] <- hw_i
  hot_only_layers[[i]] <- hot_only_i
  dry_only_layers[[i]] <- dry_only_i
}

hot_flag <- rast(hot_layers)
dry_flag <- rast(dry_layers)
wet_flag <- rast(wet_layers)
hd_flag <- rast(hd_layers)
hw_flag <- rast(hw_layers)
hot_only_flag <- rast(hot_only_layers)
dry_only_flag <- rast(dry_only_layers)

names(hot_flag) <- paste0("hot_", format(dates, "%Y_%m"))
names(dry_flag) <- paste0("dry_", format(dates, "%Y_%m"))
names(wet_flag) <- paste0("wet_", format(dates, "%Y_%m"))
names(hd_flag) <- paste0("hd_", format(dates, "%Y_%m"))
names(hw_flag) <- paste0("hw_", format(dates, "%Y_%m"))
names(hot_only_flag) <- paste0("hot_only_", format(dates, "%Y_%m"))
names(dry_only_flag) <- paste0("dry_only_", format(dates, "%Y_%m"))

time(hot_flag) <- dates
time(dry_flag) <- dates
time(wet_flag) <- dates
time(hd_flag) <- dates
time(hw_flag) <- dates
time(hot_only_flag) <- dates
time(dry_only_flag) <- dates

writeRaster(hot_flag, file.path(output_dir, "hot_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(dry_flag, file.path(output_dir, "dry_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(wet_flag, file.path(output_dir, "wet_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(hd_flag, file.path(output_dir, "hd_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(hw_flag, file.path(output_dir, "hw_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(hot_only_flag, file.path(output_dir, "hot_only_flag_1985_2025.tif"), overwrite = TRUE)
writeRaster(dry_only_flag, file.path(output_dir, "dry_only_flag_1985_2025.tif"), overwrite = TRUE)

make_event_table <- function(flag_rast, dates, out_csv, event_type) {
  if (file.exists(out_csv)) file.remove(out_csv)

  bs <- blocks(flag_rast)
  nc <- ncol(flag_rast)

  readStart(flag_rast)
  on.exit(readStop(flag_rast), add = TRUE)

  first_write <- TRUE

  for (b in seq_len(bs$n)) {
    v <- readValues(flag_rast, row = bs$row[b], nrows = bs$nrows[b], mat = TRUE)
    cell_start <- cellFromRowCol(flag_rast, bs$row[b], 1)
    cell_end <- cellFromRowCol(flag_rast, bs$row[b] + bs$nrows[b] - 1, nc)
    cells <- seq(cell_start, cell_end)

    out_list <- vector("list", nrow(v))
    k <- 0L

    for (i in seq_len(nrow(v))) {
      x <- v[i, ]
      if (all(is.na(x))) next

      x <- ifelse(is.na(x), 0L, as.integer(x > 0))
      if (!any(x == 1L)) next

      r <- rle(x)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1L
      keep <- which(r$values == 1L)

      if (length(keep) == 0) next

      k <- k + 1L
      out_list[[k]] <- data.table(
        cell = cells[i],
        event_no = seq_along(keep),
        start_index = starts[keep],
        end_index = ends[keep],
        duration = r$lengths[keep],
        start_date = as.character(dates[starts[keep]]),
        end_date = as.character(dates[ends[keep]]),
        event_type = event_type
      )
    }

    out_list <- out_list[seq_len(k)]

    if (length(out_list) > 0) {
      out_dt <- rbindlist(out_list)
      fwrite(out_dt, out_csv, append = !first_write)
      first_write <- FALSE
    }
  }
}

make_event_table(hd_flag, dates, file.path(output_dir, "hd_events_1985_2025.csv"), "HD")
make_event_table(hw_flag, dates, file.path(output_dir, "hw_events_1985_2025.csv"), "HW")
make_event_table(hot_only_flag, dates, file.path(output_dir, "hot_only_events_1985_2025.csv"), "HOT_ONLY")
make_event_table(dry_only_flag, dates, file.path(output_dir, "dry_only_events_1985_2025.csv"), "DRY_ONLY")

hd_count <- app(hd_flag, fun = function(x) {
  x <- ifelse(is.na(x), 0L, as.integer(x > 0))
  r <- rle(x)
  sum(r$values == 1L)
})

hw_count <- app(hw_flag, fun = function(x) {
  x <- ifelse(is.na(x), 0L, as.integer(x > 0))
  r <- rle(x)
  sum(r$values == 1L)
})

hd_months <- app(hd_flag, fun = function(x) sum(x > 0, na.rm = TRUE))
hw_months <- app(hw_flag, fun = function(x) sum(x > 0, na.rm = TRUE))

names(hd_count) <- "hd_event_count"
names(hw_count) <- "hw_event_count"
names(hd_months) <- "hd_month_count"
names(hw_months) <- "hw_month_count"

writeRaster(hd_count, file.path(output_dir, "hd_event_count.tif"), overwrite = TRUE)
writeRaster(hw_count, file.path(output_dir, "hw_event_count.tif"), overwrite = TRUE)
writeRaster(hd_months, file.path(output_dir, "hd_month_count.tif"), overwrite = TRUE)
writeRaster(hw_months, file.path(output_dir, "hw_month_count.tif"), overwrite = TRUE)

saveRDS(
  list(
    dates = dates,
    t_hot = t_hot,
    p_dry = p_dry,
    p_wet = p_wet
  ),
  file.path(output_dir, "compound_thresholds_1985_2025.rds")
)
