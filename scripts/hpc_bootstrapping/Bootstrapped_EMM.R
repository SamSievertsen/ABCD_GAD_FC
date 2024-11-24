library(dplyr)
library(stats)
library(lme4)
library(lmerTest)
library(sampling)
library(rsample)
library(modelbased)
library(multcomp)
library(emmeans)
library(boot)

## Load Necessary Data
# Set-up: Merge in analysis data
analysis_five_grouped_imaging_data <- read.csv("./data/ABCD/Sam_ABCD_Project/Data/analysis_five_grouped_imaging_data.csv")

# Set-up: make numerical values numeric
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data %>%
  mutate_at(vars(8:15), as.numeric)

# Set-up: make factor values factor type
analysis_five_grouped_imaging_data$eventname <- as.factor(analysis_five_grouped_imaging_data$eventname)
analysis_five_grouped_imaging_data$group <- as.factor(analysis_five_grouped_imaging_data$group)
analysis_five_grouped_imaging_data$sex <- as.factor(analysis_five_grouped_imaging_data$sex)
analysis_five_grouped_imaging_data$site_name <- as.factor(analysis_five_grouped_imaging_data$site_name)
analysis_five_grouped_imaging_data$subjectkey <- as.factor(analysis_five_grouped_imaging_data$subjectkey)

# Set-up: Remove any rows with NA or empty values and retain subjects who have both baseline and followup scans only (to perform subtraction)
analysis_five_grouped_imaging_data <- na.omit(analysis_five_grouped_imaging_data)
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data[rowSums(analysis_five_grouped_imaging_data == "") == 0, ]
analysis_five_grouped_imaging_data <-  analysis_five_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") &
           any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()

# Set-up: Set the reference level of the diagnostic group to "controls" and the eventname to "baseline" 
analysis_five_grouped_imaging_data$group <- relevel(analysis_five_grouped_imaging_data$group, ref = "Control")
analysis_five_grouped_imaging_data$eventname <- relevel(analysis_five_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")

## Function Objective:
# The bootstrap_analysis function is designed to perform bootstrap resampling on a specified statistical model fitted to a dataset. Its purpose is to estimate model parameters, particularly focusing on computing estimated marginal means (EMMs) and their associated contrasts across different groups and time points.
## Input:
# model_formula: The formula representing the statistical model to be fitted. This includes the response variable and its predictors.
# data: The dataset containing the variables used in the model, including the response variable, predictors, and grouping/time variables.
# n_bootstraps: Number of bootstrap samples to draw from the data for estimating model parameters (10,000).
## Output:
# The output of the function is a list of bootstrap results, which contains estimates of interest (such as EMMs) computed from each bootstrap sample.
# The bootstrap_analysis function internally calls a custom bootstrap function (bootstrap_func) that performs the following steps:
# Randomly samples subjects with replacement from the dataset until at least one subject is sampled for each group defined in the dataset.
# Constructs a new dataset (sampled_data) by selecting rows corresponding to the sampled subjects and specific time points (e.g., baseline and follow-up).
# Fits a specified linear mixed-effects model (lmer) using the lmerTest package to the sampled data.
# Computes estimated marginal means (EMMs) using the emmeans package for each group at different time points, adjusting for false discovery rate (FDR).
# Stores the computed EMMs and their contrasts in a data frame (emm_df) and extracts these estimates into a vector (emm_vector) for further analysis.
# Finally, the bootstrap_analysis function conducts the bootstrap resampling procedure using the boot function, applying the bootstrap_func to each bootstrap sample (R = n_bootstraps). The resulting list of bootstrapped estimates (boot_results) is returned as the output of the function

## Emmeans Function
suppressMessages({bootstrap_analysis <- function(model_formula, data, n_bootstraps = 10000) {
  # Function to fit the model, compute EMMs, and post-hoc contrasts
  bootstrap_func <- function(data, indices) {
    # Initialize a flag to check if at least one subject has been sampled for each group
    all_groups_sampled <- FALSE
    # Repeat sampling until at least one subject is sampled for each group
    while (!all_groups_sampled) {
      # Sample subjects with replacement
      sampled_subjects <- sample(unique(data$subjectkey), size = 2074, replace = TRUE)
      # Create an empty dataframe to store the sampled data
      sampled_data <- data.frame(matrix(ncol = ncol(data), nrow = 0))
      # For each sampled subject, select one row for each time point
      for (subj in sampled_subjects) {
        # Select one row for baseline
        baseline_row <- data[data$subjectkey == subj & data$eventname == "baseline_year_1_arm_1", ]
        # Select one row for year 1 follow-up
        year2_followup_row <- data[data$subjectkey == subj & data$eventname == "2_year_follow_up_y_arm_1", ]
        # Combine the selected rows
        combined_row <- rbind(baseline_row, year2_followup_row)
        # Add the combined row to the sampled data
        sampled_data <- rbind(sampled_data, combined_row)
      }
      # Check if at least one subject has been sampled for each group
      all_groups_sampled <- all(c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter") %in% sampled_data$group)
      # If at least one subject is sampled for each group, proceed with model fitting
      if (all_groups_sampled) {
        # Print the number of unique values in the "group" variable
        cat("Number of unique values in 'group' variable:", length(unique(sampled_data$group)), "\n")
        group_levels <- levels(sampled_data$group)
        cat("Unique values in 'group' variable:", group_levels, "\n")
        # Print the dimensions of sampled_data
        print(dim(sampled_data))
        print(str(sampled_data))
        # Fit the model
        model <- lmerTest::lmer(model_formula, data = sampled_data)
        # Compute estimated marginal means
        emm <- emmeans(model, pairwise ~ eventname | group, adjust = "FDR")
        emm_df <- data.frame(emm$contrasts)
        emm_vector <- c("Control" = NA, "Continuous GAD" = NA, "GAD Converter" = NA, "GAD Remitter" = NA)
      } else {
        cat("Resampling required: Not all groups have at least one sampled subject.\n")
      }
    }
    # Loop through the groups in the specified order
    for (group_name in names(emm_vector)) {
      # Check if there is an estimate for the current group
      if (group_name %in% emm_df$group) {
        # Get the estimate for the current group
        estimate <- emm_df$estimate[emm_df$group == group_name]
        # Store the estimate in the appropriate position in the emm_vector
        emm_vector[group_name] <- estimate
      } else {
        # Print a message if the group is not found in the df
        message(paste("No estimate found for", group_name))
      }
    }
    emm_vector <- emm_df$estimate
    return(emm_vector)
  }
  # Perform bootstrapping
  boot_results <- boot(data, bootstrap_func, R = n_bootstraps)
  return(boot_results)
 }
})
# Models and their formulas
models <- list(Within_VAN_repeated_measures_MEM = "rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey)")
# Perform bootstrap analysis for each model
suppressMessages({results <- list()
for (model_name in names(models)) {
  model_formula <- models[[model_name]]
  boot_results <- bootstrap_analysis(model_formula, analysis_five_grouped_imaging_data)
  results[[model_name]] <- boot_results
}})
results[[model_name]] <- boot_results
boot_ci <- lapply(1:4, function(i) {
  boot.ci(boot_results, index = i)
})

# Save both the bootstrapped model results and bootstrapped confidence intervals as R data objects
save(boot_results, file = paste0("./Data/ABCD_Analysis_Five_Bootstrapped_EMM_050124.RData"))
save(boot_ci, file = paste0("./Data/ABCD_Analysis_Five_Bootstrapped_EMM_CI_050124.RData"))
