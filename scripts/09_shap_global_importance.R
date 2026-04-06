library(data.table)
library(ranger)
library(fastshap)

model_path <- "PATH/TO/OUTPUT_DIRECTORY/rf_env_model.rds"
data_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_model_dataset_complete.csv"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"

seed <- 123
total_sample_n <- 6000
nsim <- 200

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rf_env <- readRDS(model_path)
dt <- fread(data_path)

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

label_map <- data.table(
  variable = predictors,
  predictor_label = c(
    "Antecedent vegetation state",
    "Climatic water balance during HD months (pr - PET)",
    "Maximum temperature during HD months",
    "Elevation (m a.s.l.)",
    "Vapour pressure deficit during HD months",
    "Precipitation during HD months",
    "Soil moisture constraint (event minimum)",
    "Event duration (months)"
  )
)

dt <- dt[
  !is.na(resistance) &
    !is.na(elev_band) &
    complete.cases(dt[, ..predictors])
]

dt[, elev_band := factor(
  elev_band,
  levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)]

set.seed(seed)

bands <- levels(dt$elev_band)
per_band <- ceiling(total_sample_n / length(bands))

shap_sample <- rbindlist(lapply(bands, function(b) {
  x <- dt[elev_band == b]
  if (!nrow(x)) return(NULL)
  n_take <- min(per_band, nrow(x))
  x[sample(.N, n_take)]
}), use.names = TRUE, fill = TRUE)

if (nrow(shap_sample) > total_sample_n) {
  shap_sample <- shap_sample[sample(.N, total_sample_n)]
}

X <- as.data.frame(dt[, ..predictors])
X_shap <- as.data.frame(shap_sample[, ..predictors])

pred_fun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}

set.seed(seed)

shap_mat <- fastshap::explain(
  object = rf_env,
  X = X,
  newdata = X_shap,
  pred_wrapper = pred_fun,
  nsim = nsim,
  adjust = TRUE
)

shap_dt <- as.data.table(shap_mat)
shap_dt[, event_id := shap_sample$event_id]
shap_dt[, pixel_id := shap_sample$pixel_id]
shap_dt[, resistance := shap_sample$resistance]
shap_dt[, elev_band := shap_sample$elev_band]
shap_dt[, predicted := pred_fun(rf_env, X_shap)]

setcolorder(
  shap_dt,
  c("event_id", "pixel_id", "elev_band", "resistance", "predicted", predictors)
)

fwrite(
  shap_dt,
  file.path(output_dir, "rf_env_shap_values.csv")
)

global_shap <- data.table(
  variable = predictors,
  mean_abs_shap = vapply(shap_mat[, predictors, drop = FALSE], function(x) mean(abs(x), na.rm = TRUE), numeric(1))
)[order(-mean_abs_shap)]

global_shap <- merge(global_shap, label_map, by = "variable", all.x = TRUE)
global_shap[, rank := seq_len(.N)]

setcolorder(
  global_shap,
  c("rank", "predictor_label", "variable", "mean_abs_shap")
)

fwrite(
  global_shap,
  file.path(output_dir, "Table_S8_global_shap_importance_rf_env.csv")
)

top4 <- global_shap$variable[1:min(4, nrow(global_shap))]

dependence_dt <- rbindlist(lapply(top4, function(v) {
  data.table(
    event_id = shap_sample$event_id,
    pixel_id = shap_sample$pixel_id,
    elev_band = shap_sample$elev_band,
    variable = v,
    feature_value = shap_sample[[v]],
    shap_value = shap_mat[, v]
  )
}), use.names = TRUE)

fwrite(
  dependence_dt,
  file.path(output_dir, "rf_env_shap_dependence_top_predictors.csv")
)

direction_dt <- rbindlist(lapply(predictors, function(v) {
  data.table(
    variable = v,
    feature_value = shap_sample[[v]],
    shap_value = shap_mat[, v]
  )
}), use.names = TRUE)

fwrite(
  direction_dt,
  file.path(output_dir, "rf_env_shap_direction_distribution.csv")
)
