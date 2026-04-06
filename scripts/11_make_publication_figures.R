library(data.table)
library(terra)
library(ggplot2)

cv_path <- "PATH/TO/OUTPUT_DIRECTORY/rf_grouped_cv_predictions.csv"
shap_global_path <- "PATH/TO/OUTPUT_DIRECTORY/Table_S8_global_shap_importance_rf_env.csv"
shap_direction_path <- "PATH/TO/OUTPUT_DIRECTORY/rf_env_shap_direction_distribution.csv"
shap_dependence_path <- "PATH/TO/OUTPUT_DIRECTORY/rf_env_shap_dependence_top_predictors.csv"
pixel_summary_path <- "PATH/TO/OUTPUT_DIRECTORY/all_pixel_summary.csv"
amplification_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_amplification_pixel_summary.csv"

map_p50_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_resistance_p50.tif"
map_unc_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_resistance_uncertainty_p90_p10.tif"
map_wb_path <- "PATH/TO/OUTPUT_DIRECTORY/wb_typical_hd.tif"
map_amp_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_amplification.tif"

boundary_path <- "PATH/TO/STUDY_AREA_BOUNDARY.gpkg"
output_dir <- "PATH/TO/FIGURE_OUTPUT_DIRECTORY"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rmse_fun <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae_fun <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2_fun <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  sse <- sum((obs - pred)^2)
  sst <- sum((obs - mean(obs))^2)
  1 - sse / sst
}

q25 <- function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q75 <- function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)

cv <- fread(cv_path)
shap_global <- fread(shap_global_path)
shap_dir <- fread(shap_direction_path)
shap_dep <- fread(shap_dependence_path)
px <- fread(pixel_summary_path)
amp <- fread(amplification_path)

cv[, model_label := fifelse(
  model == "RF_env", "RF (environment-only)", "RF (environment + x-y)"
)]

cv[, residual := predicted - observed]

px <- px[event_type == "HD"]
px[, elev_band := factor(
  elev_band,
  levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)]

amp[, elev_band := factor(
  fifelse(
    is.na(y), NA_character_, NA_character_
  ),
  levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)]

if ("band" %in% names(amp)) {
  amp[, elev_band := factor(
    band,
    levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
  )]
} else if ("elev_band" %in% names(amp)) {
  amp[, elev_band := factor(
    elev_band,
    levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
  )]
}

cv_scatter <- function(d, label, file_name) {
  met <- d[, .(
    RMSE = rmse_fun(observed, predicted),
    MAE = mae_fun(observed, predicted),
    R2 = r2_fun(observed, predicted)
  )]

  p <- ggplot(d, aes(x = observed, y = predicted)) +
    geom_point(alpha = 0.25, size = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(
      x = "Observed resistance (NDVI anomaly z-score)",
      y = "Predicted resistance (NDVI anomaly z-score)",
      title = label,
      subtitle = paste0(
        "R² = ", sprintf("%.3f", met$R2),
        "; RMSE = ", sprintf("%.3f", met$RMSE),
        "; MAE = ", sprintf("%.3f", met$MAE)
      )
    ) +
    theme_bw(base_size = 12)

  ggsave(
    filename = file.path(output_dir, file_name),
    plot = p,
    width = 7,
    height = 6,
    dpi = 400
  )
}

cv_scatter(
  cv[model == "RF_env"],
  "Grouped cross-validation: RF (environment-only)",
  "Figure_5A_RF_env_scatter.png"
)

cv_scatter(
  cv[model == "RF_env_xy"],
  "Grouped cross-validation: RF (environment + x-y)",
  "Figure_5B_RF_env_xy_scatter.png"
)

p_resid <- ggplot(
  cv,
  aes(x = elev_band, y = residual)
) +
  geom_boxplot(outlier.alpha = 0.15) +
  facet_wrap(~model_label, ncol = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "Elevation band (m a.s.l.)",
    y = "Residual (predicted - observed)",
    title = "Residual distributions across elevation bands"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_5CD_residuals_by_elevation.png"),
  p_resid,
  width = 8,
  height = 8,
  dpi = 400
)

shap_global[, predictor_label := factor(
  predictor_label,
  levels = rev(shap_global[order(mean_abs_shap)]$predictor_label)
)]

p_shap_bar <- ggplot(
  shap_global,
  aes(x = predictor_label, y = mean_abs_shap)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean |SHAP|",
    title = "Global SHAP importance for HD resistance"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_6C_global_shap_importance.png"),
  p_shap_bar,
  width = 8,
  height = 6,
  dpi = 400
)

label_map <- shap_global[, .(variable, predictor_label)]

shap_dir <- merge(shap_dir, label_map, by = "variable", all.x = TRUE)
shap_dir[, predictor_label := factor(
  predictor_label,
  levels = rev(shap_global[order(mean_abs_shap)]$predictor_label)
)]

p_shap_dir <- ggplot(
  shap_dir,
  aes(x = predictor_label, y = shap_value)
) +
  geom_jitter(aes(color = feature_value), width = 0.22, height = 0, alpha = 0.35, size = 0.6) +
  coord_flip() +
  labs(
    x = NULL,
    y = "SHAP value",
    color = "Feature\nvalue",
    title = "Direction and variability of SHAP effects"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_6D_shap_direction_distribution.png"),
  p_shap_dir,
  width = 9,
  height = 6,
  dpi = 400
)

shap_dep <- merge(shap_dep, label_map, by = "variable", all.x = TRUE)
shap_dep[, predictor_label := factor(
  predictor_label,
  levels = unique(shap_dep$predictor_label)
)]

p_dep <- ggplot(
  shap_dep,
  aes(x = feature_value, y = shap_value)
) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
  facet_wrap(~predictor_label, scales = "free_x", ncol = 2) +
  labs(
    x = "Predictor value",
    y = "SHAP value",
    title = "SHAP dependence for leading predictors"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_S3_shap_dependence_top_predictors.png"),
  p_dep,
  width = 10,
  height = 8,
  dpi = 400
)

p_freq <- ggplot(px, aes(x = elev_band, y = n_events)) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(
    x = "Elevation band (m a.s.l.)",
    y = "HD event frequency",
    title = "Hot-dry event frequency by elevation"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_HD_frequency_by_elevation.png"),
  p_freq,
  width = 8,
  height = 5,
  dpi = 400
)

p_res <- ggplot(px, aes(x = elev_band, y = mean_resistance)) +
  geom_boxplot(outlier.alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "Elevation band (m a.s.l.)",
    y = "Resistance (NDVI anomaly z-score)",
    title = "HD resistance by elevation"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_HD_resistance_by_elevation.png"),
  p_res,
  width = 8,
  height = 5,
  dpi = 400
)

p_rec <- ggplot(px, aes(x = elev_band, y = mean_recovery)) +
  geom_boxplot(outlier.alpha = 0.15) +
  labs(
    x = "Elevation band (m a.s.l.)",
    y = "Recovery time (months)",
    title = "HD recovery by elevation"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(output_dir, "Figure_HD_recovery_by_elevation.png"),
  p_rec,
  width = 8,
  height = 5,
  dpi = 400
)

if (nrow(amp) && "amplification" %in% names(amp)) {
  p_amp <- ggplot(amp, aes(x = elev_band, y = amplification)) +
    geom_boxplot(outlier.alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
      x = "Elevation band (m a.s.l.)",
      y = "Amplification",
      title = "HD amplification by elevation"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(output_dir, "Figure_HD_amplification_by_elevation.png"),
    p_amp,
    width = 8,
    height = 5,
    dpi = 400
  )
}

rast_to_df <- function(path, layer_name) {
  r <- rast(path)
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  setDT(d)
  setnames(d, old = names(d)[3], new = "value")
  d[, layer := layer_name]
  d
}

map_df <- rbindlist(list(
  rast_to_df(map_p50_path, "Median predicted resistance (p50)"),
  rast_to_df(map_unc_path, "Predictive dispersion (p90 - p10)"),
  rast_to_df(map_wb_path, "Typical HD climatic water balance"),
  rast_to_df(map_amp_path, "HD amplification")
), use.names = TRUE)

if (file.exists(boundary_path)) {
  bnd <- vect(boundary_path)
  bnd_df <- as.data.frame(crds(bnd, df = TRUE))
  bnd_geom <- bnd
} else {
  bnd_geom <- NULL
}

p_maps <- ggplot() +
  geom_raster(data = map_df, aes(x = x, y = y, fill = value)) +
  facet_wrap(~layer, scales = "free", ncol = 2) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Spatial patterns of HD response and drivers"
  ) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(output_dir, "Figure_7_prediction_surfaces.png"),
  p_maps,
  width = 10,
  height = 8,
  dpi = 400
)
