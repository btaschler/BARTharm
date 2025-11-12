# Function to save data in various formats
# - data: The data object to be saved (dataframe, matrix, vector, etc.)
# - file_name: The base name of the file (without extension)
# - saving_path: The directory path where the file will be saved
# - save_format: The format to save the data in ("RData", "rds", "csv", "tsv")

saving_data <- function(data, file_name, saving_path, save_format = "RData") {
  
  # Ensure saving_path ends with a slash
  if (!grepl("/$", saving_path)) {
    saving_path <- paste0(saving_path, "/")
  }
  
  if (save_format == "RData") {
    save(data, file = paste0(saving_path, file_name, ".RData"))
  } else if (save_format == "rds") {
    saveRDS(data, file = paste0(saving_path, file_name, ".rds"))
  } else if (save_format == "csv") {
    write.csv(data, file = paste0(saving_path, file_name, ".csv"), row.names = FALSE)
  } else if (save_format == "tsv") {
    write.table(data, file = paste0(saving_path, file_name, ".tsv"), sep = "\t", row.names = FALSE)
  } else {
    stop("Unsupported save format. Please use 'RData', 'rds', csv', or 'tsv'.")
  }
  
  message("Data saved in ", save_format, " format to: ", saving_path)
}
