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

# Set-up: make factor values factor type
analysis_five_grouped_imaging_data$eventname <- as.factor(analysis_five_grouped_imaging_data$eventname)
analysis_five_grouped_imaging_data$group <- as.factor(analysis_five_grouped_imaging_data$group)
analysis_five_grouped_imaging_data$sex <- as.factor(analysis_five_grouped_imaging_data$sex)
analysis_five_grouped_imaging_data$site_name <- as.factor(analysis_five_grouped_imaging_data$site_name)
analysis_five_grouped_imaging_data$subjectkey <- as.factor(analysis_five_grouped_imaging_data$subjectkey)

# Set-up: Set the reference level of the diagnostic group to "controls" and the eventname to "baseline" 
analysis_five_grouped_imaging_data$group <- relevel(analysis_five_grouped_imaging_data$group, ref = "Control")
analysis_five_grouped_imaging_data$eventname <- relevel(analysis_five_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")

### Emmeans linear slopes contrast function 
# Define a function to compute estimated marginal means and post-hoc contrasts using bootstrapping
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
        # Compute post-hoc contrasts
        contrasts <- contrast(emm[[1]], interaction = c("poly", "pairwise"), by = NULL)
        # Store the contrast slope_contrast_results in a df
        contrasts_df <- data.frame(contrasts)
        # Create a named vector to store the pairwise estimates derived from the contrast
        contrasts_vector <- c("Control - Continuous GAD" = NA, "Control - GAD Converter" = NA, "Control - GAD Remitter" = NA, "Continuous GAD - GAD Converter" = NA, "Continuous GAD - GAD Remitter" = NA, "GAD Converter - GAD Remitter" = NA)
      } else {
        cat("Resampling required: Not all groups have at least one sampled subject.\n")
      }
    }
    # Loop through the groups in the specified order
    for (group_name in names(contrasts_vector)) {
      # Check if there is an estimate for the current group
      if (group_name %in% contrasts_df$group_pairwise) {
        # Get the estimate for the current group
        estimate <- contrasts_df$estimate[contrasts_df$group_pairwise == group_name]
        # Store the estimate in the appropriate position in the contrasts_vector
        contrasts_vector[group_name] <- estimate
      } else {
        # Print a message if the group is not found in the df
        message(paste("No estimate found for", group_name))
      }
    }
    contrasts_vector <- contrasts_df$estimate
    return(contrasts_vector)
  }
  # Perform bootstrapping
  boot_slope_contrast_results <- boot(data, bootstrap_func, R = n_bootstraps)
  return(boot_slope_contrast_results)
}
})
# Models and their formulas
models <- list(Within_VAN_repeated_measures_MEM = "rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey)")
# Perform bootstrap analysis for each model
suppressMessages({slope_contrast_results <- list()
for (model_name in names(models)) {
  model_formula <- models[[model_name]]
  boot_slope_contrast_results <- bootstrap_analysis(model_formula, analysis_five_grouped_imaging_data)
  slope_contrast_results[[model_name]] <- boot_slope_contrast_results
}})
# Store all bootstrapped model slope_contrast_results in the boot list object
slope_contrast_results[[model_name]] <- boot_slope_contrast_results
# Bootstrap confidence intervals based on the distribution of estimates
slope_contrast_boot_ci <- lapply(1:6, function(i) {
  boot.ci(boot_slope_contrast_results, index = i)
})


# Save both the bootstrapped model results and bootstrapped confidence intervals as R data objects
save(slope_contrast_results, file = paste0("./Data/ABCD_Analysis_Five_Bootstrapped_EMM_Slope_Contrast_050124.RData"))
save(slope_contrast_boot_ci, file = paste0("./Data/ABCD_Analysis_Five_Bootstrapped_EMM_Slope_Contrast_CI_050124.RData"))
