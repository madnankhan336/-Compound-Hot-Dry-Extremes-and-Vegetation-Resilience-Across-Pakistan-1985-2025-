library(terra)

input_dir <- "PATH/TO/LANDSAT_YEARLY_NDVI_STACKS"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir,
  pattern = "^pakistan_landsat_monthly_ndvi_\\d{4}\\.tif$",
  full.names = TRUE
)

stopifnot(length(files) > 0)

years <- as.integer(sub(".*_(\\d{4})\\.tif$", "\\1", basename(files)))
files <- files[order(years)]

ndvi <- rast(files)

dates <- seq(
  as.Date(sprintf("%d-01-01", min(years))),
  as.Date(sprintf("%d-12-01", max(years))),
  by = "month"
)

stopifnot(nlyr(ndvi) == length(dates))

names(ndvi) <- paste0("NDVI_", format(dates, "%Y_%m"))
time(ndvi) <- dates

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  ndvi <- crop(ndvi, aoi)
  ndvi <- mask(ndvi, aoi)
}

month_id <- as.integer(format(dates, "%m"))

clim_mean <- tapp(ndvi, index = month_id, fun = mean, na.rm = TRUE)
clim_sd <- tapp(ndvi, index = month_id, fun = sd, na.rm = TRUE)

names(clim_mean) <- paste0("mean_", sprintf("%02d", 1:12))
names(clim_sd) <- paste0("sd_", sprintf("%02d", 1:12))

z_layers <- vector("list", length(dates))

for (m in 1:12) {
  idx <- which(month_id == m)
  sd_m <- ifel(clim_sd[[m]] == 0, NA, clim_sd[[m]])
  z_block <- (ndvi[[idx]] - clim_mean[[m]]) / sd_m

  for (j in seq_along(idx)) {
    z_layers[[idx[j]]] <- z_block[[j]]
  }
}

z_ndvi <- rast(z_layers)
names(z_ndvi) <- paste0("NDVI_Z_", format(dates, "%Y_%m"))
time(z_ndvi) <- dates

writeRaster(
  clim_mean,
  file.path(output_dir, "ndvi_monthly_climatology_mean.tif"),
  overwrite = TRUE
)

writeRaster(
  clim_sd,
  file.path(output_dir, "ndvi_monthly_climatology_sd.tif"),
  overwrite = TRUE
)

writeCDF(
  z_ndvi,
  filename = file.path(output_dir, "ndvi_monthly_zscore_1985_2025.nc"),
  varname = "ndvi_z",
  longname = "Monthly Landsat NDVI anomaly z-score",
  unit = "z-score",
  overwrite = TRUE
)

for (y in unique(format(dates, "%Y"))) {
  idx <- which(format(dates, "%Y") == y)

  writeRaster(
    z_ndvi[[idx]],
    file.path(output_dir, paste0("ndvi_zscore_", y, ".tif")),
    overwrite = TRUE
  )
}
