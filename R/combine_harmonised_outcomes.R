combine_harmonised_outcomes <- function(df_harmonised, saving_path, save_format = "RData"){
  # Combine all harmonized outcomes into a single dataframe
  harmonised_files_raw <- list.files(path = saving_path, pattern = "harmonised_.*_raw.RData", full.names = TRUE)
  cat("Found harmonized files in 0-1 scale: ", harmonised_files_raw, "\n")
  harmonised_files_original <- list.files(path = saving_path, pattern = "harmonised_.*_original.RData", full.names = TRUE)
  cat("Found harmonized files in original outcome scale: ", harmonised_files_original, "\n")

  predicted_files_raw <- list.files(path = saving_path, pattern = "predicted_.*_raw.RData", full.names = TRUE)
  cat("Found predicted files in 0-1 scale: ", predicted_files_raw, "\n")
  predicted_files_original <- list.files(path = saving_path, pattern = "predicted_.*_original.RData", full.names = TRUE)
  cat("Found predicted files in original outcome scale: ", predicted_files_original, "\n")
  
all_files <- unique(c(predicted_files_raw, predicted_files_original, harmonised_files_raw, harmonised_files_original))

for (file in all_files) {
    fname <- basename(file)
    env <- new.env()
    load(file, envir = env)

    if (grepl("^predicted_.*_raw\\.RData$", fname)) {
        feature_name <- sub("predicted_(.*)_raw\\.RData", "\\1", fname)
        if (exists("y_pred", envir = env)) df_harmonised[[feature_name]] <- env$y_pred
        else warning(sprintf("y_pred not found in %s", file))
    } else if (grepl("^predicted_.*_original\\.RData$", fname)) {
        feature_name <- sub("predicted_(.*)_original\\.RData", "\\1", fname)
        if (exists("y_pred_original", envir = env)) df_harmonised[[feature_name]] <- env$y_pred_original
        else warning(sprintf("y_pred_original not found in %s", file))
    } else if (grepl("^harmonised_.*_raw\\.RData$", fname)) {
        feature_name <- sub("harmonised_(.*)_raw\\.RData", "\\1", fname)
        if (exists("y_harmonised", envir = env)) df_harmonised[[feature_name]] <- env$y_harmonised
        else warning(sprintf("y_harmonised not found in %s", file))
    } else if (grepl("^harmonised_.*_original\\.RData$", fname)) {
        feature_name <- sub("harmonised_(.*)_original\\.RData", "\\1", fname)
        if (exists("y_harmonised_original", envir = env)) df_harmonised[[feature_name]] <- env$y_harmonised_original
        else warning(sprintf("y_harmonised_original not found in %s", file))
    } else {
        warning(sprintf("Unrecognized file pattern: %s", fname))
    }
}
  
  df_harmonised <- as.data.frame(df_harmonised)
  
  # Save the combined harmonized dataframe
  saving_data(df_harmonised, "df_combined_harmonised_realdata", saving_path, save_format = save_format)
  
  return(df_harmonised)
}