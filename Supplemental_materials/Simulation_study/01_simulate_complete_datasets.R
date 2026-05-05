library(dplyr)
library(rblimp)
# (Adjust the file path if your Blimp installation is located elsewhere)
rblimp::set_blimp("C:/Program Files/Blimp/blimp.exe")
blimp_path <- "C:/Program Files/Blimp/blimp.exe"
# Adjust working directory
setwd("C:/Monte Carlo Simulation")


# -------------------------------------------------------------------------
# 1. SETUP
# -------------------------------------------------------------------------
n_sims <- 1000
base_seed <- 12345

cat("Starting data generation...\n")

# -------------------------------------------------------------------------
# 2. DATA GENERATION LOOP
# -------------------------------------------------------------------------
for (i in 1:n_sims) {
  
  current_seed <- base_seed + i
  file_name <- paste0("complete_data_", i, ".csv")
  
  blimp_syntax <- paste0("
SIMULATE:
    id2(2) = 100;      
    id1(1) = 50 * 100; 

    define:
        b_w = 0.4;     
        b_b = 0.6;     
        var_u0 = 1.0;  
        var_e = 1.0;   
        aux_b_weight = 0.6; 
        aux_w_weight = 1.2; 

    id2:
        x_b = normal(0, 1);
        u0  = normal(0, var_u0);
        y_b = x_b * b_b + u0;
        aux_b = normal(aux_b_weight * y_b + aux_b_weight * x_b, 1.0);

    id1:
        x_w = normal(0, 1);
        y_w = normal(x_w * b_w, var_e);
        aux_w = normal(aux_w_weight * y_w - aux_w_weight * x_w, 1.0);
        x = x_w + x_b;
        y = y_w + y_b; 
        aux = aux_w + aux_b;

VARIABLES: id2 y x aux;
CLUSTERID: id2;
SEED: ", current_seed, ";
SAVE: dataset = ", file_name, ";
  ")
  
  # Write the syntax to a temporary .imp file
  writeLines(blimp_syntax, "temp_generation.imp")
  
  # Send the file directly to the Blimp executable
  system(paste(shQuote(blimp_path), "temp_generation.imp"))

  if (i %% 50 == 0) {
    cat("Successfully generated dataset:", i, "\n")
  }
}

file.remove("temp_generation.imp")
cat("Data generation complete!\n")



############################################################
############missing data rate#################################

n_sims <- 1000
total_expected_missing <- 0

cat("Calculating average missing rate...\n")

for (i in 1:n_sims) {
  # Read the complete dataset
  file_name <- paste0("complete_data_", i, ".csv")
  dat <- read.csv(file_name)
  
  # Calculate the exact probability of missingness used in the simulation
  prob_missing <- plogis(-0.5 + 2.5 * dat$aux)
  
  # Add the mean probability of this dataset to the running total
  total_expected_missing <- total_expected_missing + mean(prob_missing)
}

# Calculate the grand average across all 1000 datasets
grand_average_missing <- total_expected_missing / n_sims

cat("--------------------------------------------------\n")
cat("The overall average missingness rate is:", round(grand_average_missing * 100, 1), "%\n")
cat("--------------------------------------------------\n")