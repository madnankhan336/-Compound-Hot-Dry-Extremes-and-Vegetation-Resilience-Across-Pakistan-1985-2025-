scripts_dir <- "PATH/TO/YOUR/SCRIPTS_DIRECTORY"

scripts <- c(
  "03_ndvi_monthly_zscore.R",
  "04_detect_compound_events.R",
  "05_event_metrics_resistance_recovery.R",
  "06_elevation_band_summaries.R",
  "07_prepare_hd_model_dataset.R",
  "08_rf_grouped_cv_permutation_importance.R",
  "09_shap_global_importance.R",
  "10_qrf_typical_hd_prediction_surfaces.R",
  "11_make_publication_figures.R",
  "12_robustness_sensitivity_analysis.R"
)

script_paths <- file.path(scripts_dir, scripts)

missing_scripts <- script_paths[!file.exists(script_paths)]
if (length(missing_scripts) > 0) {
  stop(
    paste(
      "These script files are missing:",
      paste(missing_scripts, collapse = "\n")
    )
  )
}

run_log <- data.frame(
  script = scripts,
  start_time = as.POSIXct(NA),
  end_time = as.POSIXct(NA),
  status = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(script_paths)) {
  run_log$start_time[i] <- Sys.time()

  status_i <- tryCatch({
    source(script_paths[i], echo = FALSE, chdir = TRUE)
    "completed"
  }, error = function(e) {
    message("\nError in: ", scripts[i])
    message(e$message)
    "failed"
  })

  run_log$end_time[i] <- Sys.time()
  run_log$status[i] <- status_i

  if (status_i == "failed") {
    break
  }
}

print(run_log)

write.csv(
  run_log,
  file = file.path(scripts_dir, "pipeline_run_log.csv"),
  row.names = FALSE
)
