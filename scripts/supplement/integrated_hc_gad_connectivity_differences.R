## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(tidyr)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(car)
library(ggseg)
library(ggseg3d)
library(ggsegGordon)
library(gridGraphics) 
library(patchwork)
options(digits = 8, scipen = 999) 

# Read in required data 
# Merged GAD + HC Sample
GAD_HC_sample <- read.csv("./data_processed/supplement/GAD_HC_resampled_merged_groups.csv")

# Family ID Data
family_ID_data <-  read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(src_subject_id, rel_family_id))

# Resting state fMRI Data
cleaned_qcd_rsfMRI_data <- read.csv("./data_processed/rsfMRI_data_qcd_unsubset.csv")


## Data Wrangling ##

#1. Merge the family ID data into the GAD + HC sample
#1.1 Join the family ID data with the GAD + HC sample
GAD_HC_sample <- left_join(GAD_HC_sample, family_ID_data)

#1.2 Retain columns of interest in the desired order
GAD_HC_sample_for_merging <- GAD_HC_sample %>% 
  dplyr::select(c(src_subject_id, rel_family_id, interview_age, sex, site_name, eventname, group, GAD_timepoint, analysis_group, subgroup))

#2. Join the rsfMRI data with the HC + GAD sample
GAD_HC_group_connectivity_analysis_data <- left_join(GAD_HC_sample_for_merging, cleaned_qcd_rsfMRI_data)

#3. Prep data for analyses
#3.1 Convert relevant columns to numeric type
GAD_HC_group_connectivity_analysis_data <- mutate_at(GAD_HC_group_connectivity_analysis_data, vars(11:105), as.numeric)

#3.21 Convert relevant columns to factor type
GAD_HC_group_connectivity_analysis_data <- mutate_at(GAD_HC_group_connectivity_analysis_data, vars(c(rel_family_id, sex, site_name, eventname, group, analysis_group)), as.factor)

#3.22 Set the reference level of the group variable for analysis purposes
GAD_HC_group_connectivity_analysis_data$group <- relevel(GAD_HC_group_connectivity_analysis_data$group, ref = "HC")


## Analysis ##

#1. Analyze connectivity differences in metrics of interest between the current GAD and HC groups 
#1.1 Establish the range of the dependent variables 
GAD_HC_connectivity_analysis_dp_col_range <- 12:105

#1.2 Initialize an empty dataframe to store analysis values
GAD_HC_connectivity_analysis_raw_results <- data.frame(
  column_name = character(), 
  IV = character(),
  estimate = numeric(), 
  std_error = numeric(), 
  t_value = numeric(), 
  f_value = numeric(), 
  df = numeric(), 
  residual_df = numeric(), 
  p_value = numeric(), 
  stringsAsFactors = FALSE)

#1.3 Run the linear model for the site-visit group 
for (column_number in GAD_HC_connectivity_analysis_dp_col_range) {
  
  #1.31 Get the column name
  dependent_variable <- colnames(GAD_HC_group_connectivity_analysis_data)[column_number]
  print(dependent_variable)
  
  #1.32 Fit the linear regression model
  GAD_HC_connectivity_analysis_lm <- lmerTest::lmer(GAD_HC_group_connectivity_analysis_data[, dependent_variable] ~ group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|rel_family_id), na.action = na.omit, data = GAD_HC_group_connectivity_analysis_data)
  
  #1.33 Run the ANCOVA Model on the linear regression to extract omnibus group effect(s)
  GAD_HC_connectivity_analysis_ANCOVA <- car::Anova(GAD_HC_connectivity_analysis_lm, type = "II", test.statistic = "F")
  
  #1.34 Loop through fixed effects in the model to extract relevant results
  for (independent_variable in row.names(GAD_HC_connectivity_analysis_ANCOVA)) {
    if (independent_variable == "rsfmri_c_ngd_meanmotion") {
      
      #1.341 Extract + store parameters of interest for continuous variables
      summary_lm <- summary(GAD_HC_connectivity_analysis_lm)
      estimate <- summary_lm$coefficients[independent_variable, "Estimate"]
      std_error <- summary_lm$coefficients[independent_variable, "Std. Error"]
      t_value <- summary_lm$coefficients[independent_variable, "t value"]
      p_value <- summary_lm$coefficients[independent_variable, "Pr(>|t|)"]
      GAD_HC_connectivity_analysis_raw_results <- rbind(GAD_HC_connectivity_analysis_raw_results, data.frame(
        column_name = dependent_variable,
        IV = independent_variable,
        estimate = estimate,
        std_error = std_error,
        t_value = t_value,
        f_value = NA,
        df = NA,
        residual_df = NA,
        p_value = p_value))
    } else {
      
      #1.342 Extract + store parameters of interest for categorical variables
      f_value <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "F"]
      df <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Df"]
      residual_df <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Df.res"]
      p_value <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Pr(>F)"] 
      GAD_HC_connectivity_analysis_raw_results <- rbind(GAD_HC_connectivity_analysis_raw_results, data.frame(
        column_name = dependent_variable,
        IV = independent_variable,
        estimate = NA,
        std_error = NA,
        t_value = NA,
        f_value = f_value,
        df = df,
        residual_df = residual_df,
        p_value = p_value
      ))
    }
  }
}


## Clean Analysis Results ##

#1. Isolate the GAD vs HC group analysis results, FDR correct, determine/store significant results, and re-join with the rest of the results
#1.1 Subset the data based on the relevant DX vs CN IV variable strings
GAD_HC_connectivity_analysis_p_adjust <- subset(GAD_HC_connectivity_analysis_raw_results, grepl("group", GAD_HC_connectivity_analysis_raw_results$IV))

#1.21 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
GAD_HC_connectivity_analysis_p_adjust$p_adjusted <- p.adjust(GAD_HC_connectivity_analysis_p_adjust$p_value, method = "fdr")

#1.22 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
GAD_HC_connectivity_analysis_significant_results <- subset(GAD_HC_connectivity_analysis_p_adjust, GAD_HC_connectivity_analysis_p_adjust$p_adjusted <= 0.05)

#1.23 Specify the significant FC metrics
GAD_HC_connectivity_analysis_significant_fc_metrics <- GAD_HC_connectivity_analysis_significant_results$column_name
  
#1.31 Subset the original results dataframe to include only the rows where 'column_name' matches significant FC metrics
GAD_HC_connectivity_analysis_significant_results_full_models <- GAD_HC_connectivity_analysis_raw_results[GAD_HC_connectivity_analysis_raw_results$column_name %in% GAD_HC_connectivity_analysis_significant_fc_metrics, ]

#1.32 Merge the FDR corrected p values back into the full significant model results where applicable
GAD_HC_connectivity_analysis_significant_results_full_models <- left_join(GAD_HC_connectivity_analysis_significant_results_full_models, GAD_HC_connectivity_analysis_significant_results)

#1.33 Merge the FDR corrected p values back into the full all model results where applicable
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values <- left_join(GAD_HC_connectivity_analysis_raw_results, GAD_HC_connectivity_analysis_p_adjust)

#1.4 Pivot the full model significant results data to be wider (for copying into results tables)
#1.41 Convert necessary columns to numeric to avoid type issues
GAD_HC_connectivity_analysis_significant_results_full_models <- GAD_HC_connectivity_analysis_significant_results_full_models %>%
  mutate(across(c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted), as.numeric))

#1.42 Pivot the full model all results data to be wider (for copying into results tables)
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <- GAD_HC_connectivity_analysis_results_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = column_name, 
    names_from = IV, 
    values_from = c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted),
    names_glue = "{IV}_{.value}")

#1.5 Subset columns of interest and store them in the desired order
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <-
  GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted %>%
  dplyr::select(
    c(column_name,
      group_f_value,
      group_df,
      group_residual_df,
      group_p_value,
      group_p_adjusted,
      rsfmri_c_ngd_meanmotion_estimate,
      rsfmri_c_ngd_meanmotion_t_value,
      rsfmri_c_ngd_meanmotion_std_error,
      rsfmri_c_ngd_meanmotion_p_value,
      sex_f_value,
      sex_df,
      sex_residual_df,
      sex_p_value,
      eventname_f_value,
      eventname_df,
      eventname_residual_df,
      eventname_p_value))


## Output ##

#1. Save the full model results as a csv file
write.csv(GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted, "./results/GAD_HC_connectivity_analysis_results.csv", row.names = FALSE)
