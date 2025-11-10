
# This function normalizes the data and extracts the features to be harmonized.
# It applies quantile normalization to biological and IQM covariates while keeping the outcome variable.
#
# Arguments:
# - data_bio: Data frame of biological covariates
# - data_iqm: Data frame of IQM covariates
# - outcomes_col: Vector of column names corresponding to outcome variables
# - id_col: name of subject ID column

normalise_data <- function(data_bio, data_iqm, outcomes_col, id_col, site_col, var_scaling){
  
  # Apply quantile normalization to biological covariates, excluding num_ID and outcome columns
  norm_data_bio <- as.data.frame(quantile_normalize_bart(data_bio[, -which(names(data_bio) %in% c(id_col, outcomes_col))]))
  
  # Apply quantile normalization to IQM covariates, excluding num_ID
  norm_data_iqm <- as.data.frame(quantile_normalize_bart(data_iqm[, -which(names(data_iqm) %in% c(id_col, site_col))]))
  
  if(length(site_col) > 0){
    print("Adding site information to IQM data")
    cat("Number of sites:", length(unique(data_iqm[[site_col]])), "\n")
    norm_data_iqm[[site_col]] <- as.numeric(data_iqm[[site_col]])  # Re-add site column
  }

  # Keep only the outcome variables and num_ID
  Y <- data_bio %>%
    dplyr::select(outcomes_col, id_col)

  Y_outcomes <- as.matrix(Y[, 1:(ncol(Y)-1)])

  # Get means and sds for each outcome column
  original_means <- colMeans(Y_outcomes)
  original_sds <- colSds(Y_outcomes)

  Y_norm <- sweep(Y_outcomes, 2, original_means, FUN = "-")
  Y_norm <- sweep(Y_norm, 2, original_sds, FUN = "/")
  
  # Return normalized data
  return(list("norm_data_bio" = norm_data_bio, "norm_data_iqm" = norm_data_iqm, "Y" = Y, "Y_norm" = Y_norm,
              "original_means" = original_means, "original_sds" = original_sds))
}