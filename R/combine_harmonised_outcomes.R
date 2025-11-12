#' Combine Harmonized Outcomes
#' This function combines all harmonized outcomes saved in the specified directory into a single dataframe.
#' It looks for files matching the patterns "harmonised_.*_raw.RData", '"harmonised_.*_original.RData", 
#' "predicted_.*_raw.RData", and "predicted_.*_original.RData".
#' The combined dataframe is then saved to disk.
#'

combine_harmonised_outcomes <- function(saving_path, save_format = "RData"){
  
  env <- new.env()
  load(file = paste0(saving_path, "/filtered_realdata_df.", save_format), envir = env)
  df_harmonised <- env[[ls(env)[1]]]
  cat("Loaded original dataframe with dimensions: ", dim(df_harmonised), "\n")
  print(head(df_harmonised))
  
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
    y <- env[[ls(env)[1]]]
    
    if (grepl("^predicted_.*_raw\\.RData$", fname)) {
      feature_name <- sub("(predicted_.*_raw)\\.RData", "\\1", fname)
      cat("Adding feature: ", feature_name, "\n")
      df_harmonised[[feature_name]] <- y
    } else if (grepl("^predicted_.*_original\\.RData$", fname)) {
      feature_name <- sub("(predicted_.*_original)\\.RData", "\\1", fname)
      cat("Adding feature: ", feature_name, "\n")
      df_harmonised[[feature_name]] <- y
    } else if (grepl("^harmonised_.*_raw\\.RData$", fname)) {
      feature_name <- sub("(harmonised_.*_raw)\\.RData", "\\1", fname)
      cat("Adding feature: ", feature_name, "\n")
      df_harmonised[[feature_name]] <- y
    } else if (grepl("^harmonised_.*_original\\.RData$", fname)) {
      feature_name <- sub("(harmonised_.*_original)\\.RData", "\\1", fname)
      cat("Adding feature: ", feature_name, "\n")
      df_harmonised[[feature_name]] <- y
    } else {
      warning(sprintf("Unrecognized file pattern: %s", fname))
    }
  }
  
  df_harmonised <- as.data.frame(df_harmonised)
  
  # Save the combined harmonized dataframe
  saving_data(df_harmonised, "df_combined_harmonised_realdata", saving_path, save_format = save_format)
  
  return(df_harmonised)
}
