## Setup ##

# Load packages for loading, wrangling, mutating, and visualizing data and set environmental variables
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(lme4)
library(lmerTest)
library(polycor)
library(Cairo)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(modelbased)
library(cowplot)
suppressWarnings({library(gridExtra)})
suppressWarnings({library(sjPlot)})
suppressWarnings({library(sjmisc)})
suppressWarnings({library(sjlabelled)})
options(scipen = 999, digits = 6)

# Set-up: Create a dataframe containing GAD subjectkey's and their diagnosis status (baseline dx, followup dx, or both)
analysis_four_data <- full_join(baseline_GAD, followup_GAD)
both_GAD <- both_GAD_random_assignment %>%
dplyr::select(-group)
both_GAD$group <- rep("GAD_Both")
analysis_four_data <- full_join(analysis_four_data, both_GAD)
analysis_four_data <- analysis_four_data %>%
dplyr::select(-c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"))
write.csv(analysis_four_data, "C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_GAD_subjectkeys_groups.csv", row.names = FALSE)
analysis_four_controls <- control_sample %>%
dplyr::select(c(subjectkey, group))
write.csv(analysis_four_controls, "C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_CN_subjectkeys_groups.csv", row.names = FALSE)
# Set-up: Read in analysis four (analysis 2 group connectivity difference) subjectkey's and merge
analysis_four_GAD_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_GAD_subjectkeys_groups.csv")
analysis_four_CN_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_CN_subjectkeys_groups.csv")
analysis_four_sample_groups <- full_join(analysis_four_CN_sub_groups, analysis_four_GAD_sub_groups)
# Set-up: Join the grouped subjectkeys with the rsfMRI imaging data. Keep only the vars of interest 
analysis_four_grouped_imaging_data <- merge(analysis_four_sample_groups, ABCD_rsfMRI_Data_no_desc)
analysis_four_grouped_imaging_data$age_in_years <- floor(as.numeric(as.character(analysis_four_grouped_imaging_data$interview_age)) / 12)
analysis_four_grouped_imaging_data <- analysis_four_grouped_imaging_data %>% 
  dplyr::select(c(subjectkey, eventname, group, interview_age, age_in_years, sex, site_name, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))
analysis_four_grouped_imaging_data<- analysis_four_grouped_imaging_data %>%
  mutate_at(vars(8:13), as.numeric)
# Set-up: Remove any rows with NA or empty values and retain subjects who have both baseline and followup scans only (to perform subtraction)
analysis_four_grouped_imaging_data <- na.omit(analysis_four_grouped_imaging_data)
analysis_four_grouped_imaging_data <- analysis_four_grouped_imaging_data[rowSums(analysis_four_grouped_imaging_data == "") == 0, ]
analysis_four_grouped_imaging_data <-  analysis_four_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") &
           any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()



# Set-up: bring in CBCL data and remove first row (contains variable descriptions)
analysis_four_CBCL_data <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/abcd_cbcls01.csv")
analysis_four_CBCL_data <- analysis_four_CBCL_data[-1,]
analysis_four_CBCL_data$interview_age <- as.numeric(analysis_four_CBCL_data$interview_age)
analysis_four_CBCL_data <- analysis_four_CBCL_data %>%
  mutate_at(vars(10:90), as.numeric)
# Set-up: create a copy of the GAD / Control sample subjects/data from analysis 2
analysis_four_repeat_sample <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/site_visit_analysis_data.csv")
analysis_four_repeat_sample$interview_age <- as.numeric(analysis_four_repeat_sample$interview_age)
# Set-up: Load in and join the ABCD family_id data to the analysis data
abcd_family_id_data <- read.delim("C:/Users/Sam Sievertsen/Desktop/ABCD/ABCC_Package_1203705_Tabulated-Data/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, rel_family_id))
names(abcd_family_id_data)[names(abcd_family_id_data) == "rel_family_id"] <- "family_id"
analysis_four_repeat_sample <- left_join(analysis_four_repeat_sample, abcd_family_id_data)
# Set-up: combine subjects case & control status data with the resting state statistics yielded to be significant and CBCL data
analysis_four_imaging_CBCL_data <- left_join(analysis_four_repeat_sample, analysis_four_CBCL_data)
analysis_four_imaging_CBCL_data <- analysis_four_imaging_CBCL_data %>%
  dplyr::select(c(subjectkey, eventname, family_id, analysis_group, group, interview_age, age_in_years, sex, site_name, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_syn_anxdep_t, cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t))
analysis_four_imaging_CBCL_data$subjectkey <- as.factor(analysis_four_imaging_CBCL_data$subjectkey)
analysis_four_imaging_CBCL_data$family_id <- as.factor(analysis_four_imaging_CBCL_data$family_id)
analysis_four_imaging_CBCL_data$analysis_group <- as.factor(analysis_four_imaging_CBCL_data$analysis_group)
analysis_four_imaging_CBCL_data$eventname <- as.factor(analysis_four_imaging_CBCL_data$eventname)
analysis_four_imaging_CBCL_data$sex <- as.factor(analysis_four_imaging_CBCL_data$sex)
analysis_four_imaging_CBCL_data$site_name <- as.factor(analysis_four_imaging_CBCL_data$site_name)
analysis_four_imaging_CBCL_data$cbcl_scr_dsm5_anxdisord_t <- as.numeric(analysis_four_imaging_CBCL_data$cbcl_scr_dsm5_anxdisord_t)
analysis_four_imaging_CBCL_data$cbcl_scr_syn_anxdep_t <- as.numeric(analysis_four_imaging_CBCL_data$cbcl_scr_syn_anxdep_t)
analysis_four_imaging_CBCL_data$cbcl_scr_syn_internal_t <- as.numeric(analysis_four_imaging_CBCL_data$cbcl_scr_syn_internal_t)
analysis_four_imaging_CBCL_data$cbcl_scr_syn_external_t <- as.numeric(analysis_four_imaging_CBCL_data$cbcl_scr_syn_external_t)

# Set-up: Exclude any subjects from the dataframe who have NA or empty values
CBCL_columns_to_check_for_empties <- c("cbcl_scr_dsm5_anxdisord_t", "cbcl_scr_syn_anxdep_t", "cbcl_scr_syn_internal_t", "cbcl_scr_syn_external_t")
for (col in CBCL_columns_to_check_for_empties) {
  analysis_four_imaging_CBCL_data <- analysis_four_imaging_CBCL_data[!(is.na(analysis_four_imaging_CBCL_data[, col]) | analysis_four_imaging_CBCL_data[, col] == "" | analysis_four_imaging_CBCL_data[, col] == " "), ]
}

# # For each network from analysis 2 (significance), is there an association between connectivity and symptom level? First, check whether group membership predicts symptoms
# #3. Create a columns range for the CBCL columns 
# analysis_four_sanity_dp_col_range <- 15:17
# #3.1 Create an empty dataframe to store analysis values
# analysis_four_sanity_results_df <- data.frame(column_name = character(), p_value = numeric(), beta = numeric())
# #3.11 Run the lm for the site-visit group 
# for (col_num in analysis_four_sanity_dp_col_range) {
#   #3.12 Get the column name
#   analysis_col_name <- colnames(analysis_four_imaging_CBCL_data)[col_num]
#   print(analysis_col_name)
#   #3.13 Run the linear regression model
#   analysis_four_sanity_lm_analysis <- lm(analysis_four_imaging_CBCL_data[,analysis_col_name] ~ relevel(analysis_group, ref = "control") + site_name + sex + eventname, na.action = na.omit, data = analysis_four_imaging_CBCL_data)
#   #3.14 Get the p-value and beta for the column
#   print(summary(analysis_four_sanity_lm_analysis))
#   analysis_four_sanity_analysis_p_value <- summary(analysis_four_sanity_lm_analysis)$coefficients[,4]
#   analysis_four_sanity_analysis_beta <- summary(analysis_four_sanity_lm_analysis)$coefficients[,1]
#   #3.15 Store the results of the lm
#   analysis_four_sanity_results_df <- rbind(analysis_four_sanity_results_df, data.frame(column_name = analysis_col_name, p_value = analysis_four_sanity_analysis_p_value, beta = analysis_four_sanity_analysis_beta))
# }
# #3.21 Add a new column in the results df containing the independent variables
# analysis_four_sanity_results_df$IV <- rownames(analysis_four_sanity_results_df)
# #3.22 Subset the data based on the relevant DX vs CN IV variable strings
# analysis_four_sanity_p_adjust_subset <- subset(analysis_four_sanity_results_df, grepl("group", analysis_four_sanity_results_df$IV))
# #3.23 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
# analysis_four_sanity_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_sanity_p_adjust_subset$p_value, method = "fdr")
# #3.24 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
# analysis_four_sanity_significant_results <- subset(analysis_four_sanity_p_adjust_subset, analysis_four_sanity_p_adjust_subset$p_adjusted <= 0.05)
# 
# # Test whether symptom predicts connectivity
# #4.11 Create a columns range for the CBCL columns
# cbcl_col_range <- 17:19
# connectivity_col_range <- 10:16
# #4.12 Create an empty dataframe to store analysis values
# analysis_four_symptom_results_df <- data.frame(column_name = character(), p_value = numeric(), beta = numeric())
# #4.13 Run the lm for each combination of CBCL metric and connectivity metric
# for (cbcl_col_num in cbcl_col_range) {
#   #4.131 Get the CBCL metric column name
#   cbcl_col_name <- colnames(analysis_four_imaging_CBCL_data)[cbcl_col_num]
#   for (connectivity_col_num in connectivity_col_range) {
#     #4.132 Get the connectivity metric column name
#     connectivity_col_name <- colnames(analysis_four_imaging_CBCL_data)[connectivity_col_num]
#     #4.133 Run the linear regression model
#     analysis_four_symptom_lm_analysis <- lm(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, " + site_name + sex + eventname"), na.action = na.omit, data = analysis_four_imaging_CBCL_data)
#     #4.144 Get the p-value and beta for the column
#     analysis_four_symptom_analysis_p_value <- summary(analysis_four_symptom_lm_analysis)$coefficients[, 4]
#     analysis_four_symptom_analysis_beta <- summary(analysis_four_symptom_lm_analysis)$coefficients[, 1]
#     #4.145 Store the results of the lm
#     analysis_four_symptom_results_df <- rbind(
#       analysis_four_symptom_results_df,
#       data.frame(column_name = paste0(connectivity_col_name),
#                  p_value = analysis_four_symptom_analysis_p_value,
#                  beta = analysis_four_symptom_analysis_beta))
#   }
# }
# #4.21 Add a new column in the results df containing the independent variables
# analysis_four_symptom_results_df$IV <- rownames(analysis_four_symptom_results_df)
# #4.22 Subset the data based on the relevant DX vs CN IV variable strings
# analysis_four_symptom_p_adjust_subset <- subset(analysis_four_symptom_results_df, grepl("cbcl", analysis_four_symptom_results_df$IV))
# #4.23 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
# analysis_four_symptom_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_symptom_p_adjust_subset$p_value, method = "fdr")
# #4.24 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
# analysis_four_symptom_significant_results <- subset(analysis_four_symptom_p_adjust_subset, analysis_four_symptom_p_adjust_subset$p_adjusted <= 0.05)

#5.Test whether symptom level between groups predicts connectivity 
#5.11 Create a columns range for the CBCL columns
#5.11 Create a columns range for the CBCL columns
cbcl_col_range <- 16:18
connectivity_col_range <- 11:15
#5.12 Create an empty dataframe to store analysis values
analysis_four_symptom_group_int_results_df <- data.frame(
  cbcl_column_name = character(),
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
#5.13 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  #5.131 Get the CBCL metric column name
  cbcl_col_name <- colnames(analysis_four_imaging_CBCL_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #5.132 Get the connectivity metric column name
    connectivity_col_name <- colnames(analysis_four_imaging_CBCL_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    #5.133 Run the linear regression model
    analysis_four_symptom_group_int_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, "*analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = analysis_four_imaging_CBCL_data)
    #5.134 Run the ANCOVA Model
    analysis_four_symptom_group_int_ANCOVA_analysis <- car::Anova(analysis_four_symptom_group_int_mlm_analysis, type = "III", test.statistic = "F")
    #5.135 Loop through fixed effects in the model
    for (iv_name in row.names(analysis_four_symptom_group_int_ANCOVA_analysis)) {
      summary_lm <- summary(analysis_four_symptom_group_int_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #5.1351 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          analysis_four_symptom_group_int_results_df <- rbind(
            analysis_four_symptom_group_int_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
              column_name = connectivity_col_name,
              IV = iv_name,
              estimate = estimate,
              std_error = std_error,
              t_value = t_value,
              f_value = NA,
              df = NA,
              residual_df = NA,
              p_value = p_value))
        }
      } else {
        #5.1352 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- analysis_four_symptom_group_int_ANCOVA_analysis[iv_name, "F"]
          df <- analysis_four_symptom_group_int_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- analysis_four_symptom_group_int_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- analysis_four_symptom_group_int_ANCOVA_analysis[iv_name, "Pr(>F)"]
          analysis_four_symptom_group_int_results_df <- rbind(
            analysis_four_symptom_group_int_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
              column_name = connectivity_col_name,
              IV = iv_name,
              estimate = NA,
              std_error = NA,
              t_value = NA,
              f_value = f_value,
              df = df,
              residual_df = residual_df,
              p_value = p_value))
        }
      }
    }
  }
}
#5.21 Subset the data based on the relevant IV variable strings for FDR correction of p values
analysis_four_symptom_group_int_p_adjust_subset <- subset(analysis_four_symptom_group_int_results_df, grepl(":analysis_group", IV))
analysis_four_symptom_group_int_cbcl_p_adjust_subset <- subset(analysis_four_symptom_group_int_results_df, grepl("cbcl_", IV))
analysis_four_symptom_group_int_cbcl_p_adjust_subset <- analysis_four_symptom_group_int_cbcl_p_adjust_subset[is.na(analysis_four_symptom_group_int_cbcl_p_adjust_subset$f_value), ]
#5.221 Create a new column conducting an FDR (p adjustment) on the derived p values
analysis_four_symptom_group_int_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_symptom_group_int_p_adjust_subset$p_value, method = "fdr")
analysis_four_symptom_group_int_cbcl_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_symptom_group_int_cbcl_p_adjust_subset$p_value, method = "fdr")
#5.222 Merge the FDR corrected p-values together
analysis_four_symptom_group_int_p_adjust_subset <- full_join(analysis_four_symptom_group_int_p_adjust_subset, analysis_four_symptom_group_int_cbcl_p_adjust_subset)
#5.23 Create a new dataframe where only significant results are stored
analysis_four_symptom_group_int_significant_results <- subset(analysis_four_symptom_group_int_p_adjust_subset, p_adjusted <= 0.05)
#5.24 Join the FDR corrected p-values with the rest of the model results 
analysis_four_symptom_group_int_merged_adjusted_p_values <- left_join(analysis_four_symptom_group_int_results_df, analysis_four_symptom_group_int_p_adjust_subset)
#5.25 Remove all numbers from the strings in the column IV (for later merging purposes)
analysis_four_symptom_group_int_merged_adjusted_p_values$IV <- gsub("\\d+", "", analysis_four_symptom_group_int_merged_adjusted_p_values$IV)
#5.261 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
analysis_four_symptom_group_int_merged_adjusted_p_values <- analysis_four_symptom_group_int_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "analysis_group" ~ "Diagnosis Group",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint",
      grepl("^cbcl_.*:analysis_group", IV) ~ "CBCL Score*Diagnosis Group Interaction",
      TRUE ~ IV))
#5.262 Remove rows with Intercept values
analysis_four_symptom_group_int_merged_adjusted_p_values <- analysis_four_symptom_group_int_merged_adjusted_p_values[analysis_four_symptom_group_int_merged_adjusted_p_values$IV != "(Intercept)", ]
#5.27 Pivot the full model results to be wider
analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted <- analysis_four_symptom_group_int_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
<<<<<<< HEAD
#5.28 Reorder columns to match the order in your paper
analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted <- analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score*Diagnosis Group Interaction_f_value`, `CBCL Score*Diagnosis Group Interaction_df`, `CBCL Score*Diagnosis Group Interaction_residual_df`, `CBCL Score*Diagnosis Group Interaction_p_value`, `CBCL Score*Diagnosis Group Interaction_p_adjusted`, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `Diagnosis Group_f_value`, `Diagnosis Group_df`, `Diagnosis Group_residual_df`, `Diagnosis Group_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#5.29 Write the full model results as a csv file
write.csv(analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted, "C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_merged_Group_Full_Model_Results.csv")
=======
#5.28 Reorder columns to match the order in the paper
analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted <- analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score*Diagnosis Group Interaction_f_value`, `CBCL Score*Diagnosis Group Interaction_df`, `CBCL Score*Diagnosis Group Interaction_residual_df`, `CBCL Score*Diagnosis Group Interaction_p_value`, `CBCL Score*Diagnosis Group Interaction_p_adjusted`, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `CBCL Score_p_adjusted`, `Diagnosis Group_f_value`, `Diagnosis Group_df`, `Diagnosis Group_residual_df`, `Diagnosis Group_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#5.29 Write the full model results as a csv file
write.csv(analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/CBCL_merged_Group_Full_Model_Results.csv", row.names = FALSE)
>>>>>>> b8644c7 (interim code updates)
#5.3 Plot the FC ~ CBCL*DX Group Data
#5.301 Create a copy of the imaging CBCL data for plotting 
analysis_four_imaging_CBCL_data_for_plotting <- analysis_four_imaging_CBCL_data
#5.302 Alter the analysis group variable to match paper formatting
analysis_four_imaging_CBCL_data_for_plotting$analysis_group <- as.character(analysis_four_imaging_CBCL_data_for_plotting$analysis_group)
analysis_four_imaging_CBCL_data_for_plotting$analysis_group[analysis_four_imaging_CBCL_data_for_plotting$analysis_group == "control"] <- "HC"
analysis_four_imaging_CBCL_data_for_plotting$analysis_group <- as.factor(analysis_four_imaging_CBCL_data_for_plotting$analysis_group)
analysis_four_imaging_CBCL_data_for_plotting$analysis_group <- relevel(analysis_four_imaging_CBCL_data_for_plotting$analysis_group, ref = "HC")
#5.31 Define the Merged group CBCL scores and connectivity metrics
Merged_cbcl_scores <- c("cbcl_scr_syn_anxdep_t", "cbcl_scr_syn_external_t", "cbcl_scr_syn_internal_t")
Merged_cbcl_labels <- c("Anxious Depressed", "Externalizing", "Internalizing")
Merged_connectivity_metrics <- c("rsfmri_cor_ngd_cerc_scs_cdelh", "rsfmri_cor_ngd_cerc_scs_aglh", "rsfmri_c_ngd_vta_ngd_vta", "rsfmri_cor_ngd_df_scs_ptlh", "rsfmri_cor_ngd_sa_scs_ptlh")
<<<<<<< HEAD
Merged_connectivity_labels <- c("CON - Left Caudate", "CON - Left Amygdala", "Within-VAN", "DMN - Left Putamen", "SA - Left Putamen")
=======
Merged_connectivity_labels <- c("CON - Left Caudate", "CON - Left Amygdala", "Within-VAN", "DMN - Left Putamen", "SN - Left Putamen")
>>>>>>> b8644c7 (interim code updates)
#5.32 Create an empty list to store the plots
Merged_plots_list <- list()
#5.33 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(Merged_cbcl_scores)) {
  for (j in seq_along(Merged_connectivity_metrics)) {
    cbcl_score <- Merged_cbcl_scores[i]
    cbcl_label <- Merged_cbcl_labels[i]
    connectivity_metric <- Merged_connectivity_metrics[j]
    connectivity_label <- Merged_connectivity_labels[j]
<<<<<<< HEAD
    #5.331 Create each merged scatterplot
    Merged_scatterplot <- ggplot(analysis_four_imaging_CBCL_data_for_plotting, aes_string(x = cbcl_score, y = connectivity_metric, color = "analysis_group")) +
=======
    #5.331 Extract the beta estimate from the dataframe
    beta_value <- analysis_four_symptom_group_int_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric,
             cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)
    #5.332 Create each merged scatterplot
    Merged_scatterplot <- ggplot(analysis_four_imaging_CBCL_data_for_plotting, 
                                 aes_string(x = cbcl_score, y = connectivity_metric, color = "analysis_group")) +
>>>>>>> b8644c7 (interim code updates)
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      scale_color_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8, margin = margin(0.01, unit = "cm"), angle = 0),
            axis.text.y = element_text(size = 8, margin = margin(0.01, unit = "cm"), angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))
<<<<<<< HEAD
    #5.332 Store each plot in the list of plots to facet
=======
    #5.333 Annotate each plot with beta value
    Merged_scatterplot <- Merged_scatterplot +
      annotate("text",
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,     # Adjust horizontal alignment
               vjust = 1.1)   # Adjust vertical alignment
    #5.334 Store each plot in the list of plots to facet
>>>>>>> b8644c7 (interim code updates)
    Merged_plots_list[[length(Merged_plots_list) + 1]] <- Merged_scatterplot
  }
}
#5.34 Combine all plots into a facet plot
Merged_facet_plot <- do.call(grid.arrange, c(Merged_plots_list, nrow = 5))
#5.35 Create a single plot with a legend from one of the plots
<<<<<<< HEAD
Merged_facet_legend_plot <- ggplot(analysis_four_imaging_CBCL_data_for_plotting, aes_string(x = Merged_cbcl_scores[1], y = Merged_connectivity_metrics[1], color = "analysis_group")) +
=======
Merged_facet_legend_plot <- ggplot(analysis_four_imaging_CBCL_data_for_plotting, 
                                   aes_string(x = Merged_cbcl_scores[1], 
                                              y = Merged_connectivity_metrics[1], 
                                              color = "analysis_group")) +
>>>>>>> b8644c7 (interim code updates)
  geom_point() +
  scale_color_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) +
  theme(legend.position = "bottom") +
  labs(color = "Diagnostic Group")
#5.36 Extract the legend from the single plot
Merged_facet_legend <- cowplot::get_legend(Merged_facet_legend_plot)
<<<<<<< HEAD
#5.37 Combine the facet plot and the extracted legend
Merged_facet_plot_with_legend <- plot_grid(Merged_facet_plot, Merged_facet_legend, ncol = 1, rel_heights = c(1, 0.1))
print(Merged_facet_plot_with_legend)
#5.38 Save the final facet plot including the legend
ggsave("Merged_cbcl_fc_facet_plot_with_legend.png", Merged_facet_plot_with_legend, bg = "white", width = 8.5, height = 11, units = "in")
=======
#5.371 Combine the facet plot and the extracted legend
Merged_facet_plot_with_legend <- plot_grid(Merged_facet_plot, Merged_facet_legend, ncol = 1, rel_heights = c(1, 0.1))
#5.372 Print the final plot
print(Merged_facet_plot_with_legend)
#5.38 Save the final facet plot including the legend
ggsave("Merged_cbcl_fc_facet_plot_with_legend.pdf", Merged_facet_plot_with_legend, bg = "white", width = 8.5, height = 11, units = "in", device = "pdf")
>>>>>>> b8644c7 (interim code updates)

#6. Determine if there are any within GAD group differences in connectivity based on CBCL symptom dimension expression 
#6.1 First, subset just the GAD data
analysis_four_GAD_CBCL_imaging_data <- subset(analysis_four_imaging_CBCL_data, analysis_group == "GAD")
analysis_four_GAD_CBCL_imaging_data$site_name <- as.factor(analysis_four_GAD_CBCL_imaging_data$site_name)
analysis_four_GAD_CBCL_imaging_data$eventname <- as.factor(analysis_four_GAD_CBCL_imaging_data$eventname)
analysis_four_GAD_CBCL_imaging_data$sex <- as.factor(analysis_four_GAD_CBCL_imaging_data$sex)
# 6.21 Create a columns range for the CBCL columns
cbcl_col_range <- 16:18
connectivity_col_range <- 11:15
#6.22 Create an empty dataframe to store analysis values
analysis_four_GAD_results_df <- data.frame(cbcl_column_name = character(),
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
#6.23 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  #6.231 Get the CBCL metric column name
  cbcl_col_name <- colnames(analysis_four_GAD_CBCL_imaging_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #6.232 Get the connectivity metric column name
    connectivity_col_name <- colnames(analysis_four_GAD_CBCL_imaging_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    #6.233 Run the linear regression model
    analysis_four_GAD_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, " + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = analysis_four_GAD_CBCL_imaging_data)
    #6.234 Run the ANCOVA Model
    analysis_four_GAD_ANCOVA_analysis <- car::Anova(analysis_four_GAD_mlm_analysis, type = "II", test.statistic = "F")
    #6.235 Loop through fixed effects in the model
    for (iv_name in row.names(analysis_four_GAD_ANCOVA_analysis)) {
      summary_lm <- summary(analysis_four_GAD_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #6.2351 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          analysis_four_GAD_results_df <- rbind(
            analysis_four_GAD_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        #6.2352 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- analysis_four_GAD_ANCOVA_analysis[iv_name, "F"]
          df <- analysis_four_GAD_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- analysis_four_GAD_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- analysis_four_GAD_ANCOVA_analysis[iv_name, "Pr(>F)"]
          analysis_four_GAD_results_df <- rbind(
            analysis_four_GAD_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}
#6.31 Subset the data based on the relevant IV variable strings for FDR correction of p values
analysis_four_GAD_cbcl_p_adjust_subset <- subset(analysis_four_GAD_results_df, grepl("cbcl_", IV))
#6.32 Create a new column conducting an FDR (p adjustment) on the derived p values
analysis_four_GAD_cbcl_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_GAD_cbcl_p_adjust_subset$p_value, method = "fdr")
#6.33 Create a new dataframe where only significant results are stored
analysis_four_GAD_significant_results <- subset(analysis_four_GAD_cbcl_p_adjust_subset, p_adjusted <= 0.05)
#6.34 Join the FDR corrected p-values with the rest of the model results 
analysis_four_GAD_merged_adjusted_p_values <- left_join(analysis_four_GAD_results_df, analysis_four_GAD_cbcl_p_adjust_subset)
#6.35 Remove all numbers from the strings in the column IV (for later merging purposes)
analysis_four_GAD_merged_adjusted_p_values$IV <- gsub("\\d+", "", analysis_four_GAD_merged_adjusted_p_values$IV)
#6.361 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
analysis_four_GAD_merged_adjusted_p_values <- analysis_four_GAD_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))
#6.362 Remove rows with Intercept values
analysis_four_GAD_merged_adjusted_p_values <- analysis_four_GAD_merged_adjusted_p_values[analysis_four_GAD_merged_adjusted_p_values$IV != "(Intercept)", ]
#6.37 Pivot the full model results to be wider
analysis_four_GAD_merged_adjusted_p_values_pivoted <- analysis_four_GAD_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
#6.38 Reorder columns to match the order in your paper
analysis_four_GAD_merged_adjusted_p_values_pivoted <- analysis_four_GAD_merged_adjusted_p_values_pivoted %>%
<<<<<<< HEAD
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#6.354 Write the full model results as a csv file
write.csv(analysis_four_GAD_merged_adjusted_p_values_pivoted, "C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_GAD_Group_Full_Model_Results.csv")
=======
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `CBCL Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#6.354 Write the full model results as a csv file
write.csv(analysis_four_GAD_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/CBCL_GAD_Group_Full_Model_Results.csv", row.names = FALSE)
>>>>>>> b8644c7 (interim code updates)
# #6.4 Correlate and plot the results of the most significantly associated CBCL - Connectivity metrics (none survived correction)
# GAD_cbcl_scr_syn_external_t <- analysis_four_GAD_CBCL_imaging_data$cbcl_scr_syn_external_t
# GAD_rsfmri_cor_ngd_cerc_scs_ptlh <- analysis_four_GAD_CBCL_imaging_data$rsfmri_cor_ngd_cerc_scs_ptlh
# GAD_complete_cases <- complete.cases(GAD_cbcl_scr_syn_external_t, GAD_rsfmri_cor_ngd_cerc_scs_ptlh)
# CBCL_External_SA_PTLH_Correlation <- cor(GAD_cbcl_scr_syn_external_t[GAD_complete_cases], GAD_rsfmri_cor_ngd_cerc_scs_ptlh[GAD_complete_cases])
# ggplot(analysis_four_GAD_CBCL_imaging_data, aes(x = cbcl_scr_syn_external_t, y = rsfmri_cor_ngd_cerc_scs_ptlh)) +
#   geom_point() +  # scatterplot
#   geom_smooth(method = "lm", se = TRUE) +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(size = 12, angle = 0, margin = margin(0.3, unit = "cm")),
#         axis.text.y = element_text(size = 12, angle = 0),
#         axis.title.x = element_text(size = 12, margin = margin(t = 10)),
#         axis.title.y = element_text(size = 12, margin = margin(r = 10)),
#         legend.position = "top") +
#   labs(x = "CBCL Externalizing Symptoms T-Score", y = "Cingulo Opercular Network - LH Putamen Connectivity") +
#   annotate("text", x = 60,
#            y = 0.4,
#            label = paste("r =", round(CBCL_External_SA_PTLH_Correlation, 4)),
#            hjust = 0, vjust = 0)
#6.5 Plot the GAD group CBCL analysis results
#6.501 Read in the beta estimate data for models of interest 
<<<<<<< HEAD
analysis_four_GAD_merged_adjusted_p_values_pivoted <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_GAD_Group_Full_Model_Results.csv")
#6.502 Convert the beta values of interest to numeric for in plot formatting
analysis_four_GAD_merged_adjusted_p_values_pivoted$beta_cbcl_score <- as.numeric(analysis_four_GAD_merged_adjusted_p_values_pivoted$beta_cbcl_score)
=======
#analysis_four_GAD_merged_adjusted_p_values_pivoted <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_GAD_Group_Full_Model_Results.csv")
#6.502 Convert the beta values of interest to numeric for in plot formatting
analysis_four_GAD_merged_adjusted_p_values_pivoted$`CBCL Score_estimate` <- as.numeric(analysis_four_GAD_merged_adjusted_p_values_pivoted$`CBCL Score_estimate`)
>>>>>>> b8644c7 (interim code updates)
#6.51 Define the GAD group CBCL scores and connectivity metrics
GAD_cbcl_scores <-
  c("cbcl_scr_syn_anxdep_t",
    "cbcl_scr_syn_external_t",
    "cbcl_scr_syn_internal_t")
GAD_cbcl_labels <-
  c("Anxious Depressed",
    "Externalizing",
    "Internalizing")
GAD_connectivity_metrics <-
  c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh")
GAD_connectivity_labels <-
  c("CON - Left Caudate",
    "CON - Left Amygdala",
    "Within-VAN",
    "DMN - Left Putamen",
<<<<<<< HEAD
    "SA - Left Putamen")
=======
    "SN - Left Putamen")
>>>>>>> b8644c7 (interim code updates)
#6.52 Create an empty list to store the plots
GAD_plots_list <- list()
#6.53 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(GAD_cbcl_scores)) {
  for (j in seq_along(GAD_connectivity_metrics)) {
    cbcl_score <- GAD_cbcl_scores[i]
    cbcl_label <- GAD_cbcl_labels[i]
    connectivity_metric <- GAD_connectivity_metrics[j]
    connectivity_label <- GAD_connectivity_labels[j]
    #6.531 Extract the beta estimate from the dataframe
    beta_value <-
      analysis_four_GAD_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric,
<<<<<<< HEAD
             cbcl_col_name == cbcl_score) %>%
      pull(beta_cbcl_score)
=======
             cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)
>>>>>>> b8644c7 (interim code updates)
    #6.532 Create GAD_scatterplot
    GAD_scatterplot <-
      ggplot(analysis_four_GAD_CBCL_imaging_data, aes_string(x = cbcl_score, y = connectivity_metric)) +
      geom_point(color = "#D81159") +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black",
                                     size = 1,
                                     linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8,
                                       margin = margin(0.01, unit = "cm"),
                                       angle = 0),
            axis.text.y = element_text(
              size = 8,
              margin = margin(0.01, unit = "cm"),
              angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))
    #6.533 Annotate with beta value
    GAD_scatterplot <- GAD_scatterplot +
      annotate("text",
<<<<<<< HEAD
               x = 65,
               y = 0.55,
               label = paste("β =", round(as.numeric(beta_value), 4)),
               size = 2.5,
               hjust = 0,
               vjust = 0.5)
=======
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,   # Adjust horizontal alignment
               vjust = 1.1)
>>>>>>> b8644c7 (interim code updates)
    #6.534 Store the plot in the list
    GAD_plots_list[[length(GAD_plots_list) + 1]] <-
      GAD_scatterplot
  }
}
#6.54 Combine all plots into a facet plot
GAD_facet_plot <-
  do.call(grid.arrange, c(GAD_plots_list, nrow = 5))
#6.55 Print the facet plot
print(GAD_facet_plot)
#6.56 Save the plot
<<<<<<< HEAD
ggsave("GAD_group_cbcl_fc_facet_plot.png", GAD_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300)
=======
ggsave("GAD_group_cbcl_fc_facet_plot.pdf", GAD_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300, device = "pdf")
>>>>>>> b8644c7 (interim code updates)


#7. Determine if there are any within Control group differences in connectivity based on CBCL symptom dimension expression 
#7.1 First, subset just the Control data
analysis_four_control_CBCL_imaging_data <- subset(analysis_four_imaging_CBCL_data, analysis_group == "control")
analysis_four_control_CBCL_imaging_data$site_name <- as.factor(analysis_four_control_CBCL_imaging_data$site_name)
analysis_four_control_CBCL_imaging_data$eventname <- as.factor(analysis_four_control_CBCL_imaging_data$eventname)
analysis_four_control_CBCL_imaging_data$sex <- as.factor(analysis_four_control_CBCL_imaging_data$sex)
#7.21 Create a columns range for the CBCL columns
#7.22 Create an empty dataframe to store analysis values
analysis_four_control_results_df <-
  data.frame(cbcl_column_name = character(),
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
#7.23 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  #7.231 Get the CBCL metric column name
  cbcl_col_name <- colnames(analysis_four_control_CBCL_imaging_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #7.232 Get the connectivity metric column name
    connectivity_col_name <- colnames(analysis_four_control_CBCL_imaging_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    #7.233 Run the linear regression model
    analysis_four_control_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, " + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = analysis_four_control_CBCL_imaging_data)
    #7.234 Run the ANCOVA Model
    analysis_four_control_ANCOVA_analysis <- car::Anova(analysis_four_control_mlm_analysis, type = "II", test.statistic = "F")
    #7.235 Loop through fixed effects in the model
    for (iv_name in row.names(analysis_four_control_ANCOVA_analysis)) {
      summary_lm <- summary(analysis_four_control_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #7.2351 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          analysis_four_control_results_df <- rbind(
            analysis_four_control_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        #7.2352 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- analysis_four_control_ANCOVA_analysis[iv_name, "F"]
          df <- analysis_four_control_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- analysis_four_control_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- analysis_four_control_ANCOVA_analysis[iv_name, "Pr(>F)"]
          analysis_four_control_results_df <- rbind(
            analysis_four_control_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}
#7.31 Subset the data based on the relevant IV variable strings for FDR correction of p values
analysis_four_control_cbcl_p_adjust_subset <- subset(analysis_four_control_results_df, grepl("cbcl_", IV))
#7.32 Create a new column conducting an FDR (p adjustment) on the derived p values
analysis_four_control_cbcl_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_control_cbcl_p_adjust_subset$p_value, method = "fdr")
#7.33 Create a new dataframe where only significant results are stored
analysis_four_control_significant_results <- subset(analysis_four_control_cbcl_p_adjust_subset, p_adjusted <= 0.05)
#7.34 Join the FDR corrected p-values with the rest of the model results 
analysis_four_control_merged_adjusted_p_values <- left_join(analysis_four_control_results_df, analysis_four_control_cbcl_p_adjust_subset)
#7.35 Remove all numbers from the strings in the column IV (for later merging purposes)
analysis_four_control_merged_adjusted_p_values$IV <- gsub("\\d+", "", analysis_four_control_merged_adjusted_p_values$IV)
#7.361 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
analysis_four_control_merged_adjusted_p_values <- analysis_four_control_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))
#7.362 Remove rows with Intercept values
analysis_four_control_merged_adjusted_p_values <- analysis_four_control_merged_adjusted_p_values[analysis_four_control_merged_adjusted_p_values$IV != "(Intercept)", ]
#7.37 Pivot the full model results to be wider
analysis_four_control_merged_adjusted_p_values_pivoted <- analysis_four_control_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
#7.38 Reorder columns to match the order in your paper
analysis_four_control_merged_adjusted_p_values_pivoted <- analysis_four_control_merged_adjusted_p_values_pivoted %>%
<<<<<<< HEAD
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#7.354 Write the full model results as a csv file
write.csv(analysis_four_control_merged_adjusted_p_values_pivoted, "C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_control_Group_Full_Model_Results.csv")
=======
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `CBCL Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))
#7.354 Write the full model results as a csv file
write.csv(analysis_four_control_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/CBCL_control_Group_Full_Model_Results.csv", row.names = FALSE)
>>>>>>> b8644c7 (interim code updates)
# #7.4 Plot nominally significant results
# #7.41 CO-Amygdala and Externalizing Symptoms
# control_cbcl_scr_syn_external_t <- analysis_four_control_CBCL_imaging_data$cbcl_scr_syn_external_t
# control_rsfmri_cor_ngd_cerc_scs_aglh <- analysis_four_control_CBCL_imaging_data$rsfmri_cor_ngd_cerc_scs_aglh
# control_complete_cases <- complete.cases(control_cbcl_scr_syn_external_t, control_rsfmri_cor_ngd_cerc_scs_aglh)
# CBCL_External_SA_PTLH_Correlation <- cor(control_cbcl_scr_syn_external_t[control_complete_cases], control_rsfmri_cor_ngd_cerc_scs_aglh[control_complete_cases])
# ggplot(analysis_four_control_CBCL_imaging_data, aes(x = cbcl_scr_syn_external_t, y = rsfmri_cor_ngd_cerc_scs_aglh)) +
#   geom_point() +  # scatterplot
#   geom_smooth(method = "lm", se = TRUE) +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(size = 12, angle = 0, margin = margin(0.3, unit = "cm")),
#         axis.text.y = element_text(size = 12, angle = 0),
#         axis.title.x = element_text(size = 12, margin = margin(t = 10)),
#         axis.title.y = element_text(size = 12, margin = margin(r = 10)),
#         legend.position = "top") +
#   labs(x = "CBCL Externalizing Symptoms T-Score", y = "CO Network - LH Amygdala Connectivity") +
#   annotate("text", x = 55,
#            y = 0.4,
#            label = paste("r =", round(CBCL_External_SA_PTLH_Correlation, 4)),
#            hjust = 0, vjust = 0)
# #7.42 CO Network - RH Hippocampus
# control_cbcl_scr_syn_external_t <- analysis_four_control_CBCL_imaging_data$cbcl_scr_syn_external_t
# control_rsfmri_cor_ngd_cerc_scs_hprh <- analysis_four_control_CBCL_imaging_data$rsfmri_cor_ngd_cerc_scs_hprh
# control_complete_cases <- complete.cases(control_cbcl_scr_syn_external_t, control_rsfmri_cor_ngd_cerc_scs_hprh)
# CBCL_External_SA_PTLH_Correlation <- cor(control_cbcl_scr_syn_external_t[control_complete_cases], control_rsfmri_cor_ngd_cerc_scs_hprh[control_complete_cases])
# ggplot(analysis_four_control_CBCL_imaging_data, aes(x = cbcl_scr_syn_external_t, y = rsfmri_cor_ngd_cerc_scs_hprh)) +
#   geom_point() +  # scatterplot
#   geom_smooth(method = "lm", se = TRUE) +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(size = 12, angle = 0, margin = margin(0.3, unit = "cm")),
#         axis.text.y = element_text(size = 12, angle = 0),
#         axis.title.x = element_text(size = 12, margin = margin(t = 10)),
#         axis.title.y = element_text(size = 12, margin = margin(r = 10)),
#         legend.position = "top") +
#   labs(x = "CBCL Externalizing Symptoms T-Score", y = "CO Network - RH Hippocampus Connectivity") +
#   annotate("text", x = 55,
#            y = 0.4,
#            label = paste("r =", round(CBCL_External_SA_PTLH_Correlation, 4)),
#            hjust = 0, vjust = 0)
# #7.43 Within VTA Network Anx Dep
# control_cbcl_scr_syn_anxdep_t <- analysis_four_control_CBCL_imaging_data$cbcl_scr_syn_anxdep_t
# control_rsfmri_c_ngd_vta_ngd_vta <- analysis_four_control_CBCL_imaging_data$rsfmri_c_ngd_vta_ngd_vta
# control_complete_cases <- complete.cases(control_cbcl_scr_syn_anxdep_t, control_rsfmri_c_ngd_vta_ngd_vta)
# CBCL_External_SA_PTLH_Correlation <- cor(control_cbcl_scr_syn_anxdep_t[control_complete_cases], control_rsfmri_c_ngd_vta_ngd_vta[control_complete_cases])
# ggplot(analysis_four_control_CBCL_imaging_data, aes(x = cbcl_scr_syn_anxdep_t, y = rsfmri_c_ngd_vta_ngd_vta)) +
#   geom_point() +  # scatterplot
#   geom_smooth(method = "lm", se = TRUE) +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(size = 12, angle = 0, margin = margin(0.3, unit = "cm")),
#         axis.text.y = element_text(size = 12, angle = 0),
#         axis.title.x = element_text(size = 12, margin = margin(t = 10)),
#         axis.title.y = element_text(size = 12, margin = margin(r = 10)),
#         legend.position = "top") +
#   labs(x = "CBCL Anx-Dep Symptoms T-Score", y = "Within VTA Network Connectivity") +
#   annotate("text", x = 55,
#            y = 0.4,
#            label = paste("r =", round(CBCL_External_SA_PTLH_Correlation, 4)),
#            hjust = 0, vjust = 0)
# #7.44 CO LH Putamen - Anx Dep
# control_cbcl_scr_syn_anxdep_t <- analysis_four_control_CBCL_imaging_data$cbcl_scr_syn_anxdep_t
# control_rsfmri_cor_ngd_cerc_scs_ptlh <- analysis_four_control_CBCL_imaging_data$rsfmri_cor_ngd_cerc_scs_ptlh
# control_complete_cases <- complete.cases(control_cbcl_scr_syn_anxdep_t, control_rsfmri_cor_ngd_cerc_scs_ptlh)
# CBCL_External_SA_PTLH_Correlation <- cor(control_cbcl_scr_syn_anxdep_t[control_complete_cases], control_rsfmri_cor_ngd_cerc_scs_ptlh[control_complete_cases])
# ggplot(analysis_four_control_CBCL_imaging_data, aes(x = cbcl_scr_syn_anxdep_t, y = rsfmri_cor_ngd_cerc_scs_ptlh)) +
#   geom_point() +  # scatterplot
#   geom_smooth(method = "lm", se = TRUE) +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(size = 12, angle = 0, margin = margin(0.3, unit = "cm")),
#         axis.text.y = element_text(size = 12, angle = 0),
#         axis.title.x = element_text(size = 12, margin = margin(t = 10)),
#         axis.title.y = element_text(size = 12, margin = margin(r = 10)),
#         legend.position = "top") +
#   labs(x = "CBCL Anx-Dep Symptoms T-Score", y = "CO Network - LH Putamen Connectivity") +
#   annotate("text", x = 55,
#            y = 0.4,
#            label = paste("r =", round(CBCL_External_SA_PTLH_Correlation, 4)),
#            hjust = 0, vjust = 0)
#7.5 Plot the control group CBCL analysis results
#7.501 Read in the beta estimate data for models of interest 
<<<<<<< HEAD
analysis_four_control_merged_adjusted_p_values_pivoted <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_control_Group_Full_Model_Results.csv")
#7.502 Convert the beta values of interest to numeric for in plot formatting
analysis_four_control_merged_adjusted_p_values_pivoted$beta_cbcl_score <- as.numeric(analysis_four_control_merged_adjusted_p_values_pivoted$beta_cbcl_score)
=======
#analysis_four_control_merged_adjusted_p_values_pivoted <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_control_Group_Full_Model_Results.csv")
#7.502 Convert the beta values of interest to numeric for in plot formatting
analysis_four_control_merged_adjusted_p_values_pivoted$`CBCL Score_estimate` <- as.numeric(analysis_four_control_merged_adjusted_p_values_pivoted$`CBCL Score_estimate`)
>>>>>>> b8644c7 (interim code updates)
#7.51 Define the control group CBCL scores and connectivity metrics
control_cbcl_scores <-
  c("cbcl_scr_syn_anxdep_t",
    "cbcl_scr_syn_external_t",
    "cbcl_scr_syn_internal_t")
control_cbcl_labels <-
  c("Anxious Depressed",
    "Externalizing",
    "Internalizing")
control_connectivity_metrics <-
  c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh")
control_connectivity_labels <-
  c("CON - Left Caudate",
    "CON - Left Amygdala",
    "Within-VAN",
    "DMN - Left Putamen",
<<<<<<< HEAD
    "SA - Left Putamen")
=======
    "SN - Left Putamen")
>>>>>>> b8644c7 (interim code updates)
#7.52 Create an empty list to store the plots
control_plots_list <- list()
#7.53 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(control_cbcl_scores)) {
  for (j in seq_along(control_connectivity_metrics)) {
    cbcl_score <- control_cbcl_scores[i]
    cbcl_label <- control_cbcl_labels[i]
    connectivity_metric <- control_connectivity_metrics[j]
    connectivity_label <- control_connectivity_labels[j]
    #7.531 Extract the beta estimate from the dataframe
    beta_value <-
      analysis_four_control_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric,
<<<<<<< HEAD
             cbcl_col_name == cbcl_score) %>%
      pull(beta_cbcl_score)
=======
             cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)
>>>>>>> b8644c7 (interim code updates)
    #7.532 Create control_scatterplot
    control_scatterplot <-
      ggplot(analysis_four_control_CBCL_imaging_data, aes_string(x = cbcl_score, y = connectivity_metric)) +
      geom_point(color = "#218380") +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black",
                                     size = 1,
                                     linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8,
                                       margin = margin(0.01, unit = "cm"),
                                       angle = 0),
            axis.text.y = element_text(
              size = 8,
              margin = margin(0.01, unit = "cm"),
              angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))
    #7.533 Annotate with beta value
    control_scatterplot <- control_scatterplot +
      annotate("text",
<<<<<<< HEAD
               x = 65,
               y = 0.5,
               label = paste("β =", round(as.numeric(beta_value), 4)),
               size = 2.5,
               hjust = 0,
               vjust = 0.5)
=======
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,   # Adjust horizontal alignment
               vjust = 1.1)
>>>>>>> b8644c7 (interim code updates)
    #7.534 Store the plot in the list
    control_plots_list[[length(control_plots_list) + 1]] <-
      control_scatterplot
  }
}
#7.54 Combine all plots into a facet plot
control_facet_plot <-
  do.call(grid.arrange, c(control_plots_list, nrow = 5))
#7.55 Print the facet plot
print(control_facet_plot)
#7.56 Save the plot
<<<<<<< HEAD
ggsave("control_group_cbcl_fc_facet_plot.png", control_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300)
=======
ggsave("control_group_cbcl_fc_facet_plot.pdf", control_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300, device = "pdf")
>>>>>>> b8644c7 (interim code updates)

#8.1 Create a copy of the GAD-Control merged imaging and CBCL dataframe
analysis_four_merged_CBCL_imaging_data <- analysis_four_imaging_CBCL_data
#8.21 Create a columns range for the CBCL columns
cbcl_col_range <- 16:18
connectivity_col_range <- 11:15
#8.22 Create an empty dataframe to store analysis values
analysis_four_merged_results_df <- data.frame(IV= character(), column_name = character(), cbcl_col_name = character(), f_value = numeric(), df = numeric(), residual_df = numeric(), p_value = numeric())
#8.23 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  #8.231 Get the CBCL metric column name
  cbcl_col_name <- colnames(analysis_four_merged_CBCL_imaging_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #8.232 Get the connectivity metric column name
    connectivity_col_name <- colnames(analysis_four_merged_CBCL_imaging_data)[connectivity_col_num]
    print(connectivity_col_name)
    #8.233 Run the linear regression model
    analysis_four_merged_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, "* analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = analysis_four_merged_CBCL_imaging_data)
    #8.234 Run the ANCOVA Model
    analysis_four_merged_ANCOVA_analysis <- car::Anova(analysis_four_merged_mlm_analysis, type = "III", test.statistic="F")
    #8.235 Get the f-value for the column
    analysis_four_merged_ANCOVA_analysis_f_value <- analysis_four_merged_ANCOVA_analysis$`F`
    #8.236 Get the DF and residual DF for the column 
    analysis_four_merged_ANCOVA_analysis_df <- analysis_four_merged_ANCOVA_analysis$Df
    analysis_four_merged_ANCOVA_analysis_residual_df <- analysis_four_merged_ANCOVA_analysis$Df.res
    #8.237 Get the p-value for the column
    analysis_four_merged_ANCOVA_analysis_p_value <- analysis_four_merged_ANCOVA_analysis$`Pr(>F)`
    #8.238 Store the results of the ANCOVA
    analysis_four_merged_results_df <- rbind(
      analysis_four_merged_results_df,
      data.frame(IV = rownames(analysis_four_merged_ANCOVA_analysis),
                 column_name = paste0(connectivity_col_name),
                 cbcl_col_name = cbcl_col_name,
                 f_value = analysis_four_merged_ANCOVA_analysis_f_value,
                 df = analysis_four_merged_ANCOVA_analysis_df,
                 residual_df = analysis_four_merged_ANCOVA_analysis_residual_df,
                 p_value = analysis_four_merged_ANCOVA_analysis_p_value))
  }
}
#8.31 Subset the data based on the relevant DX vs CN IV variable strings
analysis_four_merged_p_adjust_subset <- subset(analysis_four_merged_results_df, grepl(":", analysis_four_merged_results_df$IV))
#8.32 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
analysis_four_merged_p_adjust_subset$p_adjusted <- p.adjust(analysis_four_merged_p_adjust_subset$p_value, method = "fdr")
#8.33 Create a new dataframe where only significant results from the analysis variable (DX vs CN symptom interaction) are stored
analysis_four_merged_significant_results <- subset(analysis_four_merged_p_adjust_subset, analysis_four_merged_p_adjust_subset$p_adjusted <= 0.05)
#8.34 Print the model output to a csv file 
#8.341 Join the FDR corrected p-values with the rest of the model results 
analysis_four_merged_merged_adjusted_p_values <- left_join(analysis_four_merged_results_df, analysis_four_merged_p_adjust_subset)
#8.342 Remove all numbers from the strings in the column IV (for later merging purposes)
analysis_four_merged_merged_adjusted_p_values$IV <- gsub("\\d+", "", analysis_four_merged_merged_adjusted_p_values$IV)
#8.343 Pivot the full model all results data to be wider (for copying into results tables)
analysis_four_merged_merged_adjusted_p_values_pivoted <- analysis_four_merged_merged_adjusted_p_values %>%
  pivot_wider(id_cols = c(column_name, cbcl_col_name), names_from = IV, values_from = c(f_value, df, residual_df, p_value, p_adjusted))
#8.35 Write the full model results as a csv file
write.csv(analysis_four_merged_merged_adjusted_p_values_pivoted, "C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/CBCL_merged_Group_Full_Model_Results.csv", row.names = FALSE)
