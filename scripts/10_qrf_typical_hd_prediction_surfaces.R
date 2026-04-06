library(terra)
library(data.table)
library(ranger)

model_data_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_model_dataset_complete.csv"
template_raster_path <- "PATH/TO/TEMPLATE_RASTER.tif"
boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"

seed <- 123
ntrees <- 1000
quantiles_out <- c(0.10, 0.50, 0.90)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(model_data_path)
template <- rast(template_raster_path)

if (file.exists(boundary_path)) {
  aoi <- vect(boundary_path)
  template <- crop(template, aoi)
  template <- mask(template, aoi)
}

predictors <- c(
  "ndvi_pre",
  "wb_mean",
  "tmax_mean",
  "elev_m",
  "vpd_mean",
  "pr_mean",
  "soil_min",
  "duration"
)

dt <- dt[
  !is.na(cell) &
    !is.na(resistance) &
    !is.na(pixel_id) &
    complete.cases(dt[, ..predictors])
]

fit_formula <- as.formula(
  paste("resistance ~", paste(predictors, collapse = " + "))
)

qrf <- ranger(
  formula = fit_formula,
  data = dt[, c("resistance", predictors), with = FALSE],
  num.trees = ntrees,
  quantreg = TRUE,
  keep.inbag = TRUE,
  seed = seed,
  respect.unordered.factors = "order",
  write.forest = TRUE
)

saveRDS(
  qrf,
  file.path(output_dir, "qrf_hd_resistance_model.rds")
)

typical_hd <- dt[, .(
  n_events = .N,
  ndvi_pre = mean(ndvi_pre, na.rm = TRUE),
  wb_mean = mean(wb_mean, na.rm = TRUE),
  tmax_mean = mean(tmax_mean, na.rm = TRUE),
  elev_m = mean(elev_m, na.rm = TRUE),
  vpd_mean = mean(vpd_mean, na.rm = TRUE),
  pr_mean = mean(pr_mean, na.rm = TRUE),
  soil_min = mean(soil_min, na.rm = TRUE),
  duration = mean(duration, na.rm = TRUE),
  x = mean(x, na.rm = TRUE),
  y = mean(y, na.rm = TRUE)
), by = .(cell)]

typical_hd <- typical_hd[complete.cases(typical_hd[, ..predictors])]

pred_q <- predict(
  qrf,
  data = typical_hd[, ..predictors],
  type = "quantiles",
  quantiles = quantiles_out
)$predictions

pred_q <- as.data.table(pred_q)
setnames(pred_q, c("p10", "p50", "p90"))

typical_hd <- cbind(typical_hd, pred_q)
typical_hd[, uncertainty := p90 - p10]

typical_hd[, dominant_predictor_proxy := predictors[max.col(
  cbind(
    abs(ndvi_pre),
    abs(wb_mean),
    abs(tmax_mean),
    abs(elev_m),
    abs(vpd_mean),
    abs(pr_mean),
    abs(soil_min),
    abs(duration)
  ),
  ties.method = "first"
)]]

fwrite(
  typical_hd,
  file.path(output_dir, "typical_hd_predictor_summary_by_cell.csv")
)

make_raster_from_cell_values <- function(template, cells, vals, out_name) {
  r <- rast(template)
  values(r) <- NA_real_
  values(r)[cells] <- vals
  names(r) <- out_name
  r
}

make_cat_raster_from_cell_values <- function(template, cells, vals, out_name) {
  r <- rast(template)
  values(r) <- NA_integer_

  levs <- data.table(
    id = seq_along(unique(vals)),
    label = unique(vals)
  )

  id_vals <- levs$id[match(vals, levs$label)]
  values(r)[cells] <- id_vals
  names(r) <- out_name

  levels(r) <- data.frame(value = levs$id, category = levs$label)
  r
}

r_ndvi_pre <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$ndvi_pre, "ndvi_pre_typical_hd")
r_wb <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$wb_mean, "wb_typical_hd")
r_tmax <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$tmax_mean, "tmax_typical_hd")
r_vpd <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$vpd_mean, "vpd_typical_hd")
r_pr <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$pr_mean, "pr_typical_hd")
r_soil <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$soil_min, "soil_min_typical_hd")
r_duration <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$duration, "duration_typical_hd")

r_p10 <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$p10, "hd_resistance_p10")
r_p50 <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$p50, "hd_resistance_p50")
r_p90 <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$p90, "hd_resistance_p90")
r_unc <- make_raster_from_cell_values(template, typical_hd$cell, typical_hd$uncertainty, "hd_resistance_uncertainty")

r_dom <- make_cat_raster_from_cell_values(
  template,
  typical_hd$cell,
  typical_hd$dominant_predictor_proxy,
  "dominant_predictor_proxy"
)

writeRaster(r_ndvi_pre, file.path(output_dir, "ndvi_pre_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_wb, file.path(output_dir, "wb_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_tmax, file.path(output_dir, "tmax_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_vpd, file.path(output_dir, "vpd_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_pr, file.path(output_dir, "pr_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_soil, file.path(output_dir, "soil_min_typical_hd.tif"), overwrite = TRUE)
writeRaster(r_duration, file.path(output_dir, "duration_typical_hd.tif"), overwrite = TRUE)

writeRaster(r_p10, file.path(output_dir, "hd_resistance_p10.tif"), overwrite = TRUE)
writeRaster(r_p50, file.path(output_dir, "hd_resistance_p50.tif"), overwrite = TRUE)
writeRaster(r_p90, file.path(output_dir, "hd_resistance_p90.tif"), overwrite = TRUE)
writeRaster(r_unc, file.path(output_dir, "hd_resistance_uncertainty_p90_p10.tif"), overwrite = TRUE)
writeRaster(r_dom, file.path(output_dir, "dominant_predictor_proxy.tif"), overwrite = TRUE)

pred_stack <- c(r_p10, r_p50, r_p90, r_unc)
writeRaster(
  pred_stack,
  file.path(output_dir, "hd_qrf_prediction_surfaces.tif"),
  overwrite = TRUE
)
