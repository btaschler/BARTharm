# Loading libraries 
library(caret)
library(SoftBart)
library(dplyr)
library(matrixStats)
library(MCMCpack)

# Loading functions
source("R/bartharm.R")
source("R/bartharm_inference.R")
source("R/get_data.R")
source("R/simulate_data.R")
source("R/normalise_data.R")
source("R/load_data.R")
source("R/saving_data.R")
source("R/combine_harmonised_outcomes.R")

# Define path to real dataset (must be an .RData file containing a data frame)
file_path <- 'data/real_data.RData' # example Real Data file (can also use .csv or .tsv)

# Directory where results and harmonized outputs will be saved
saving_path <- "results/"

# Specify column names for biological covariates (e.g., age and sex)
bio_col <- c("Age", "Sex")

# Specify column names for image quality metrics (IQMs)
iqm_col <- c("snv", "cnr", "qi_1", "qi_2") 

# Specify column names for the outcomes to be harmonized
outcomes_col <- c("NBV1", "NBV2")

# Specify the name of the subject ID column
id_col <- c("num_ID")

# Specify the name of the scanner ID column (if available)
site_col <- c() # leave empty if scanner IDs are not available

# Specify the saving format
save_format = "RData" # can also be rds, csv or tsv

# Run BARTharm harmonization on real data
# - Loads and normalizes data
# - Runs BART-based Gibbs sampling to separate nuisance (mu) and signal (tau)
# - Saves harmonized outcomes and posterior draws

df_harmonised <- bartharm(
  simulate_data = FALSE,             # Use real data (not simulated)
  saving_path = saving_path,         # Output directory
  save_format = save_format,         # Saving format
  file_path = file_path,             # Input dataset path
  bio_col = bio_col,                 # Biological covariates
  iqm_col = iqm_col,                 # IQM covariates
  outcomes_col = outcomes_col,       # Outcome variables to harmonize
  id_col = id_col,                   # Subject ID column
  site_col = site_col,               # Site ID column
  num_iter = 500,                    # Number of MCMC iterations
  burn_in = 50,                      # Burn-in iterations
  thinning_interval = 2,             # Thinning
  num_tree_mu = 200,                 # Trees in mu forest (IQMs)
  num_tree_tau = 50,                 # Trees in tau forest (bio features)
  beta_mu = 2, beta_tau = 2,         # BART prior parameters
  gamma_mu = 0.95, gamma_tau = 0.95, # BART prior parameters
  var_scaling = FALSE                # Whether to harmonize variance (is set to TRUE, site_col must be provided)
)

# If running harmonisation sequentially, i.e., length(outcomes_col) > 1, then df_harmonised already contains all harmonized outcomes

# If running each outcome in parallelel, i.e., length(outcomes_col) = 1, you can combine the results into a single data frame with the following function:
# This has to be ran once all pararlel jobs are finished and all harmonized outcome files are saved in 'saving_path'
df_harmonised <- combine_harmonised_outcomes(saving_path, save_format)

# The final harmonized data frame is saved at: saving_path/df_combined_harmonised_realdata.<save_format>

  
