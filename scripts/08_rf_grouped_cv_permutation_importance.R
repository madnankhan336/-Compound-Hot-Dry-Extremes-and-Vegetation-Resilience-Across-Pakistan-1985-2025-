library(data.table)
library(ranger)

input_path <- "PATH/TO/OUTPUT_DIRECTORY/hd_model_dataset_complete.csv"
output_dir <- "PATH/TO/OUTPUT_DIRECTORY"

seed <- 123
nfolds <- 5
ntrees <- 1000

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_path)

dt <- dt[
  !is.na(resistance) &
    !is.na(pixel_id) &
    !is.na(ndvi_pre) &
    !is.na(wb_mean) &
    !is.na(tmax_mean) &
    !is.na(elev_m) &
    !is.na(vpd_mean) &
    !is.na(pr_mean) &
    !is.na(soil_min) &
    !is.na(duration) &
    !is.na(x) &
    !is.na(y)
]

dt[, elev_band := factor(
  elev_band,
  levels = c("<500", "500–1500", "1500–2500", "2500–3500", "≥3500")
)]

predictors_env <- c(
  "ndvi_pre",
  "wb_mean",
  "tmax_mean",
  "elev_m",
  "vpd_mean",
  "pr_mean",
  "soil_min",
  "duration"
)

predictors_xy <- c(predictors_env, "x", "y")

make_formula <- function(y, x) {
  as.formula(paste(y, "~", paste(x, collapse = " + ")))
}

rmse_fun <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae_fun <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
r2_fun <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (!length(obs)) return(NA_real_)
  sse <- sum((obs - pred)^2)
  sst <- sum((obs - mean(obs))^2)
  if (sst == 0) return(NA_real_)
  1 - sse / sst
}

metric_table <- function(obs, pred) {
  data.table(
    RMSE = rmse_fun(obs, pred),
    MAE = mae_fun(obs, pred),
    R2 = r2_fun(obs, pred)
  )
}

set.seed(seed)

u_pixels <- unique(dt$pixel_id)
u_pixels <- sample(u_pixels)

fold_id <- rep(seq_len(nfolds), length.out = length(u_pixels))
fold_map <- data.table(pixel_id = u_pixels, fold = fold_id)

setkey(dt, pixel_id)
setkey(fold_map, pixel_id)

dt <- fold_map[dt]

fit_fold_model <- function(train_dt, test_dt, predictors) {
  form <- make_formula("resistance", predictors)

  mod <- ranger(
    formula = form,
    data = train_dt[, c("resistance", predictors), with = FALSE],
    num.trees = ntrees,
    importance = "none",
    seed = seed,
    respect.unordered.factors = "order",
    write.forest = TRUE
  )

  pred <- predict(mod, data = test_dt[, predictors, with = FALSE])$predictions

  data.table(
    event_id = test_dt$event_id,
    pixel_id = test_dt$pixel_id,
    fold = test_dt$fold,
    elev_band = test_dt$elev_band,
    observed = test_dt$resistance,
    predicted = pred
  )
}

cv_env <- vector("list", nfolds)
cv_xy <- vector("list", nfolds)

for (k in seq_len(nfolds)) {
  train_dt <- dt[fold != k]
  test_dt <- dt[fold == k]

  cv_env[[k]] <- fit_fold_model(train_dt, test_dt, predictors_env)
  cv_env[[k]][, model := "RF_env"]

  cv_xy[[k]] <- fit_fold_model(train_dt, test_dt, predictors_xy)
  cv_xy[[k]][, model := "RF_env_xy"]
}

cv_env <- rbindlist(cv_env)
cv_xy <- rbindlist(cv_xy)

cv_all <- rbindlist(list(cv_env, cv_xy), use.names = TRUE)

fwrite(
  cv_all,
  file.path(output_dir, "rf_grouped_cv_predictions.csv")
)

overall_metrics <- cv_all[, metric_table(observed, predicted), by = model]

fwrite(
  overall_metrics,
  file.path(output_dir, "rf_grouped_cv_overall_metrics.csv")
)

elev_metrics <- cv_all[, metric_table(observed, predicted), by = .(model, elev_band)]

fwrite(
  elev_metrics,
  file.path(output_dir, "rf_grouped_cv_metrics_by_elevation.csv")
)

fit_full_model <- function(data, predictors) {
  ranger(
    formula = make_formula("resistance", predictors),
    data = data[, c("resistance", predictors), with = FALSE],
    num.trees = ntrees,
    importance = "permutation",
    seed = seed,
    respect.unordered.factors = "order",
    write.forest = TRUE
  )
}

rf_env <- fit_full_model(dt, predictors_env)
rf_xy <- fit_full_model(dt, predictors_xy)

perm_env <- data.table(
  variable = names(rf_env$variable.importance),
  permutation_importance = as.numeric(rf_env$variable.importance)
)[order(-permutation_importance)]

perm_xy <- data.table(
  variable = names(rf_xy$variable.importance),
  permutation_importance = as.numeric(rf_xy$variable.importance)
)[order(-permutation_importance)]

label_map <- data.table(
  variable = c(
    "ndvi_pre", "wb_mean", "tmax_mean", "elev_m",
    "vpd_mean", "pr_mean", "soil_min", "duration",
    "x", "y"
  ),
  predictor_label = c(
    "Antecedent vegetation state",
    "Climatic water balance during HD months (pr - PET)",
    "Maximum temperature during HD months",
    "Elevation (m a.s.l.)",
    "Vapour pressure deficit during HD months",
    "Precipitation during HD months",
    "Soil moisture constraint (event minimum)",
    "Event duration (months)",
    "Easting / longitude coordinate",
    "Northing / latitude coordinate"
  )
)

perm_env <- merge(perm_env, label_map, by = "variable", all.x = TRUE)
perm_xy <- merge(perm_xy, label_map, by = "variable", all.x = TRUE)

perm_env[, rank := seq_len(.N)]
perm_xy[, rank := seq_len(.N)]

setcolorder(perm_env, c("rank", "predictor_label", "variable", "permutation_importance"))
setcolorder(perm_xy, c("rank", "predictor_label", "variable", "permutation_importance"))

fwrite(
  perm_env,
  file.path(output_dir, "rf_env_permutation_importance.csv")
)

fwrite(
  perm_xy,
  file.path(output_dir, "rf_env_xy_permutation_importance.csv")
)

saveRDS(
  rf_env,
  file.path(output_dir, "rf_env_model.rds")
)

saveRDS(
  rf_xy,
  file.path(output_dir, "rf_env_xy_model.rds")
)

summary_table <- copy(overall_metrics)
summary_table[model == "RF_env", model := "RF (environment-only)"]
summary_table[model == "RF_env_xy", model := "RF (environment + x-y)"]

fwrite(
  summary_table,
  file.path(output_dir, "Table_S4_grouped_cv_overall_metrics.csv")
)

summary_elev <- copy(elev_metrics)
summary_elev[model == "RF_env", model := "RF (environment-only)"]
summary_elev[model == "RF_env_xy", model := "RF (environment + x-y)"]

fwrite(
  summary_elev,
  file.path(output_dir, "Table_S5_grouped_cv_metrics_by_elevation.csv")
)

perm_env_out <- copy(perm_env)
perm_xy_out <- copy(perm_xy)

fwrite(
  perm_env_out,
  file.path(output_dir, "Table_S6_permutation_importance_rf_env.csv")
)

fwrite(
  perm_xy_out,
  file.path(output_dir, "Table_S7_permutation_importance_rf_env_xy.csv")
)
