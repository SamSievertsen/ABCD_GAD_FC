#Setup: Load packages for loading, wrangling, mutating, and visualizing data
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(sampling)
library(rsample)
library(splitstackshape)
library(polycor)
library(Cairo)
library(circlize)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(chorddiag)
library(modelbased)
options(digits = 4)

#Setup: Load BPM Data 
ABCD_youth_BPM <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/abcd_yssbpm01.csv") 
ABCD_youth_BPM <- ABCD_youth_BPM[-1, ]
Imaging_CBCL_Data <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_imaging_CBCL_data.csv")
Control_Imaging_CBCL_Data <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_control_CBCL_imaging_data.csv")
GAD_Imaging_CBCL_Data <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_GAD_CBCL_imaging_data.csv")
abcd_family_id_data <- read.delim("/Users/samsievertsen/Desktop/ABCC_Package_1203705_Tabulated-Data/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, rel_family_id))
names(abcd_family_id_data)[names(abcd_family_id_data) == "rel_family_id"] <- "family_id"
Imaging_CBCL_Data <- left_join(Imaging_CBCL_Data, abcd_family_id_data)
Control_Imaging_CBCL_Data <- left_join(Control_Imaging_CBCL_Data, abcd_family_id_data)
GAD_Imaging_CBCL_Data <- left_join(GAD_Imaging_CBCL_Data, abcd_family_id_data)

#1. Data cleaning
#1.11 Create a cleaned (complete) version of the BPM and PoA data containing non NA/empty scores of interest that have been derived without missing questions
ABCD_youth_BPM_filtered <- ABCD_youth_BPM %>% 
  dplyr::select(c(subjectkey, eventname, bpm_y_scr_internal_t, bpm_y_scr_internal_nm, bpm_y_scr_external_t, bpm_y_scr_external_nm, poa_y_ss_sum, poa_y_ss_sum_nm))
ABCD_youth_PoA_filtered <- ABCD_youth_BPM %>% 
  dplyr::select(c(subjectkey, eventname, poa_y_ss_sum, poa_y_ss_sum_nm))
#1.12 Keep rows where the BPM and PoA scores are not NA or empty, and where the number of questions missing that went into the calculation of each BPM score is 0
ABCD_youth_BPM_complete_data <- subset(ABCD_youth_BPM_filtered, !is.na(bpm_y_scr_internal_t) & bpm_y_scr_internal_t != "" & !is.na(bpm_y_scr_external_t) & bpm_y_scr_external_t != "" & bpm_y_scr_internal_nm == 0 & bpm_y_scr_external_nm == 0)
ABCD_youth_PoA_complete_data <- subset(ABCD_youth_PoA_filtered, !is.na(poa_y_ss_sum) & poa_y_ss_sum != "" & poa_y_ss_sum_nm == 0)
#1.13 Keep only the columns of interest
ABCD_youth_BPM_clean <- ABCD_youth_BPM_complete_data %>% 
  dplyr::select(c(subjectkey, eventname, bpm_y_scr_internal_t, bpm_y_scr_external_t))
ABCD_youth_PoA_clean <- ABCD_youth_PoA_complete_data %>% 
  dplyr::select(c(subjectkey, eventname, poa_y_ss_sum))
#1.2 Merge the whole group, Control only, and GAD only CBCL/Imaging dataframes with the youth BPM data
#1.21 Whole sample 
Imaging_CBCL_BPM_Data <- left_join(Imaging_CBCL_Data, ABCD_youth_BPM_clean)
Imaging_CBCL_PoA_Data <- left_join(Imaging_CBCL_Data, ABCD_youth_PoA_clean)
#1.22 Control Sample (PoA excluded from this point on as no subjects have scores @ relevant timepoint)
Control_Imaging_CBCL_BPM_Data <- left_join(Control_Imaging_CBCL_Data, ABCD_youth_BPM_clean)
#1.23 GAD Sample (PoA excluded from this point on as no subjects have scores @ relevant timepoint)
GAD_Imaging_CBCL_BPM_Data <- left_join(GAD_Imaging_CBCL_Data, ABCD_youth_BPM_clean)
#1.3 Change the relevant columns to the correct data type for modeling
#1.31 Whole sample
Imaging_CBCL_BPM_Data$analysis_group <- as.factor(Imaging_CBCL_BPM_Data$analysis_group)
Imaging_CBCL_BPM_Data$group <- as.factor(Imaging_CBCL_BPM_Data$group)
Imaging_CBCL_BPM_Data$site_name <- as.factor(Imaging_CBCL_BPM_Data$site_name)
Imaging_CBCL_BPM_Data$family_id <- as.factor(Imaging_CBCL_BPM_Data$family_id)
Imaging_CBCL_BPM_Data$sex <- factor(Imaging_CBCL_BPM_Data$sex)
Imaging_CBCL_BPM_Data$eventname <- factor(Imaging_CBCL_BPM_Data$eventname)
Imaging_CBCL_BPM_Data$scanner_model <- as.factor(Imaging_CBCL_BPM_Data$scanner_model)
Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t <- as.numeric(Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t)
Imaging_CBCL_BPM_Data$bpm_y_scr_external_t <- as.numeric(Imaging_CBCL_BPM_Data$bpm_y_scr_external_t)
#1.32 Control Sample
Control_Imaging_CBCL_BPM_Data$analysis_group <- as.factor(Control_Imaging_CBCL_BPM_Data$analysis_group)
Control_Imaging_CBCL_BPM_Data$group <- as.factor(Control_Imaging_CBCL_BPM_Data$group)
Control_Imaging_CBCL_BPM_Data$site_name <- as.factor(Control_Imaging_CBCL_BPM_Data$site_name)
Control_Imaging_CBCL_BPM_Data$family_id <- as.factor(Control_Imaging_CBCL_BPM_Data$family_id)
Control_Imaging_CBCL_BPM_Data$sex <- factor(Control_Imaging_CBCL_BPM_Data$sex)
Control_Imaging_CBCL_BPM_Data$eventname <- factor(Control_Imaging_CBCL_BPM_Data$eventname)
Control_Imaging_CBCL_BPM_Data$scanner_model <- as.factor(Control_Imaging_CBCL_BPM_Data$scanner_model)
Control_Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t <- as.numeric(Control_Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t)
Control_Imaging_CBCL_BPM_Data$bpm_y_scr_external_t <- as.numeric(Control_Imaging_CBCL_BPM_Data$bpm_y_scr_external_t)
#1.33 GAD Sample
GAD_Imaging_CBCL_BPM_Data$analysis_group <- as.factor(GAD_Imaging_CBCL_BPM_Data$analysis_group)
GAD_Imaging_CBCL_BPM_Data$group <- as.factor(GAD_Imaging_CBCL_BPM_Data$group)
GAD_Imaging_CBCL_BPM_Data$site_name <- as.factor(GAD_Imaging_CBCL_BPM_Data$site_name)
GAD_Imaging_CBCL_BPM_Data$family_id <- as.factor(GAD_Imaging_CBCL_BPM_Data$family_id)
GAD_Imaging_CBCL_BPM_Data$sex <- factor(GAD_Imaging_CBCL_BPM_Data$sex)
GAD_Imaging_CBCL_BPM_Data$eventname <- factor(GAD_Imaging_CBCL_BPM_Data$eventname)
GAD_Imaging_CBCL_BPM_Data$scanner_model <- as.factor(GAD_Imaging_CBCL_BPM_Data$scanner_model)
GAD_Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t <- as.numeric(GAD_Imaging_CBCL_BPM_Data$bpm_y_scr_internal_t)
GAD_Imaging_CBCL_BPM_Data$bpm_y_scr_external_t <- as.numeric(GAD_Imaging_CBCL_BPM_Data$bpm_y_scr_external_t)

#2. Modeling 
#2.1 Test whether symptom predicts connectivity in the whole group
#2.11 Create a columns range for the BPM columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14
#2.12 Create an empty dataframe to store analysis values
<<<<<<< HEAD
whole_sample_BPM_con_results_df <- data.frame(IV = character(), column_name = character(), bpm_col_name = character(), f_value = numeric(), df = numeric(), residual_df = numeric(), p_value = numeric())
#2.13 Run the lm for each combination of BPM metric and connectivity metric
=======
whole_sample_BPM_con_results_df <- data.frame(
  bpm_column_name = character(),
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
#2.13 Run the lm for each combination of bpm metric and connectivity metric
>>>>>>> b8644c7 (interim code updates)
for (bpm_col_num in bpm_col_range) {
  #2.131 Get the bpm metric column name
  bpm_col_name <- colnames(Imaging_CBCL_BPM_Data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #2.132 Get the connectivity metric column name
    connectivity_col_name <- colnames(Imaging_CBCL_BPM_Data)[connectivity_col_num]
<<<<<<< HEAD
    print(connectivity_col_name)
    #2.133 Run the linear regression model
    whole_sample_BPM_con_lm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, "*analysis_group + sex + age_in_years + rsfmri_c_ngd_meanmotion + (1|site_name) + (1|family_id)"), na.action = na.omit, data = Imaging_CBCL_BPM_Data)
    #2.134 Run the ANCOVA
    whole_sample_BPM_con_ANCOVA_analysis <- car::Anova(whole_sample_BPM_con_lm_analysis, type = "III", test.statistic="F")
    #2.135 Get the f-value for the column
    whole_sample_BPM_con_ANCOVA_analysis_f_value <- whole_sample_BPM_con_ANCOVA_analysis$`F`
    #2.136 Get the DF and residual DF for the column 
    whole_sample_BPM_con_ANCOVA_analysis_df <- whole_sample_BPM_con_ANCOVA_analysis$Df
    whole_sample_BPM_con_ANCOVA_analysis_residual_df <- whole_sample_BPM_con_ANCOVA_analysis$Df.res
    #2.137 Get the p-value for the column
    whole_sample_BPM_con_ANCOVA_analysis_p_value <- whole_sample_BPM_con_ANCOVA_analysis$`Pr(>F)`
    #2.138 Store the results of the ANCOVA
    whole_sample_BPM_con_results_df <- rbind(
      whole_sample_BPM_con_results_df,
      data.frame(IV = rownames(whole_sample_BPM_con_ANCOVA_analysis),
                 column_name = paste0(connectivity_col_name),
                 bpm_col_name = bpm_col_name,
                 f_value = whole_sample_BPM_con_ANCOVA_analysis_f_value,
                 df = whole_sample_BPM_con_ANCOVA_analysis_df,
                 residual_df = whole_sample_BPM_con_ANCOVA_analysis_residual_df,
                 p_value = whole_sample_BPM_con_ANCOVA_analysis_p_value))
  }
}
#2.14 Remove all numbers from the strings in the column IV (for later merging purposes)
whole_sample_BPM_con_results_df$IV <- gsub("\\d+", "", whole_sample_BPM_con_results_df$IV)
#2.15 Subset the data based on the relevant DX vs CN IV variable strings
whole_sample_BPM_con_p_adjust_subset <- subset(whole_sample_BPM_con_results_df, grepl("t:analysis_", whole_sample_BPM_con_results_df$IV))
#2.15 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
whole_sample_BPM_con_p_adjust_subset$p_adjusted <- p.adjust(whole_sample_BPM_con_p_adjust_subset$p_value, method = "fdr")
#2.16 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
whole_sample_BPM_con_significant_results <- subset(whole_sample_BPM_con_p_adjust_subset, whole_sample_BPM_con_p_adjust_subset$p_adjusted <= 0.05)
#2.17 Merge the FDR corrected P values back into the full all model results
whole_sample_BPM_con_results_df_merged_adjusted_p_values <- left_join(whole_sample_BPM_con_results_df, whole_sample_BPM_con_p_adjust_subset)
#2.18 Pivot the full model results data to be wider (for copying into results tables)
whole_sample_BPM_con_results_df_merged_adjusted_p_values_pivoted <- whole_sample_BPM_con_results_df_merged_adjusted_p_values %>%
  pivot_wider(id_cols = c(column_name, bpm_col_name), names_from = IV, values_from = c(f_value, df, residual_df, p_value, p_adjusted))
#2.19 Write the output full model results to a csv file
write.csv(whole_sample_BPM_con_results_df_merged_adjusted_p_values_pivoted, "C:/Users/Sam Sievertsen/Desktop/SamResearch/Results/BPM_Merged_Group_Full_Model_Results.csv", row.names = FALSE)
=======
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    #2.133 Run the linear regression model
    whole_sample_BPM_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, "*analysis_group + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id)"), na.action = na.omit, data = Imaging_CBCL_BPM_Data)
    #2.134 Run the ANCOVA Model
    whole_sample_BPM_con_ANCOVA_analysis <- car::Anova(whole_sample_BPM_con_mlm_analysis, type = "III", test.statistic = "F")
    #2.135 Loop through fixed effects in the model
    for (iv_name in row.names(whole_sample_BPM_con_ANCOVA_analysis)) {
      summary_lm <- summary(whole_sample_BPM_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #2.1351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          whole_sample_BPM_con_results_df <- rbind(
            whole_sample_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
        #2.1352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- whole_sample_BPM_con_ANCOVA_analysis[iv_name, "F"]
          df <- whole_sample_BPM_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- whole_sample_BPM_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- whole_sample_BPM_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          whole_sample_BPM_con_results_df <- rbind(
            whole_sample_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
#2.14 Subset the data based on the relevant IV variable strings for FDR correction of p values
whole_sample_BPM_con_p_adjust_subset <- subset(whole_sample_BPM_con_results_df, grepl(":analysis_group", IV))
whole_sample_BPM_con_bpm_p_adjust_subset <- subset(whole_sample_BPM_con_results_df, grepl("bpm_", IV))
whole_sample_BPM_con_bpm_p_adjust_subset <- whole_sample_BPM_con_bpm_p_adjust_subset[is.na(whole_sample_BPM_con_bpm_p_adjust_subset$f_value), ]
#2.151 Create a new column conducting an FDR (p adjustment) on the derived p values
whole_sample_BPM_con_p_adjust_subset$p_adjusted <- p.adjust(whole_sample_BPM_con_p_adjust_subset$p_value, method = "fdr")
whole_sample_BPM_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(whole_sample_BPM_con_bpm_p_adjust_subset$p_value, method = "fdr")
#2.152 Merge the FDR corrected p-values together
whole_sample_BPM_con_p_adjust_subset <- full_join(whole_sample_BPM_con_p_adjust_subset, whole_sample_BPM_con_bpm_p_adjust_subset)
#2.16 Create a new dataframe where only significant results are stored
whole_sample_BPM_con_significant_results <- subset(whole_sample_BPM_con_p_adjust_subset, p_adjusted <= 0.05)
#2.17 Join the FDR corrected p-values with the rest of the model results 
whole_sample_BPM_con_merged_adjusted_p_values <- left_join(whole_sample_BPM_con_results_df, whole_sample_BPM_con_p_adjust_subset)
#2.181 Remove all numbers from the strings in the column IV (for later merging purposes)
whole_sample_BPM_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", whole_sample_BPM_con_merged_adjusted_p_values$IV)
#2.182 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
whole_sample_BPM_con_merged_adjusted_p_values <- whole_sample_BPM_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "analysis_group" ~ "Diagnosis Group",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint",
      grepl("^bpm_.*:analysis_group", IV) ~ "bpm Score*Diagnosis Group Interaction",
      TRUE ~ IV))
#2.183 Remove rows with Intercept values
whole_sample_BPM_con_merged_adjusted_p_values <- whole_sample_BPM_con_merged_adjusted_p_values[whole_sample_BPM_con_merged_adjusted_p_values$IV != "(Intercept)", ]
#2.191 Pivot the full model results to be wider
whole_sample_BPM_con_merged_adjusted_p_values_pivoted <- whole_sample_BPM_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
#2.192 Reorder columns to match the order in the paper
whole_sample_BPM_con_merged_adjusted_p_values_pivoted <- whole_sample_BPM_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score*Diagnosis Group Interaction_f_value`, `bpm Score*Diagnosis Group Interaction_df`, `bpm Score*Diagnosis Group Interaction_residual_df`, `bpm Score*Diagnosis Group Interaction_p_value`, `bpm Score*Diagnosis Group Interaction_p_adjusted`, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `Diagnosis Group_f_value`, `Diagnosis Group_df`, `Diagnosis Group_residual_df`, `Diagnosis Group_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))
#2.193 Write the full model results as a csv file
write.csv(whole_sample_BPM_con_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/bpm_merged_Group_Full_Model_Results.csv", row.names = FALSE)
>>>>>>> b8644c7 (interim code updates)


#2.2 Test for differences in predicting connectivity based on BPM symptoms in the control group
#2.21 Create a columns range for the BPM columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14
#2.22 Create an empty dataframe to store analysis values
control_group_BPM_con_results_df <- data.frame(bpm_column_name = character(),
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
#2.23 Run the lm for each combination of bpm metric and connectivity metric
for (bpm_col_num in bpm_col_range) {
  #2.231 Get the bpm metric column name
  bpm_col_name <- colnames(Control_Imaging_CBCL_BPM_Data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #2.232 Get the connectivity metric column name
    connectivity_col_name <- colnames(Control_Imaging_CBCL_BPM_Data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    #2.233 Run the linear regression model
    control_group_BPM_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, " + rsfmri_c_ngd_meanmotion + sex + (1|site_name)"), na.action = na.omit, data = Control_Imaging_CBCL_BPM_Data)
    #2.234 Run the ANCOVA Model
    control_group_BPM_con_ANCOVA_analysis <- car::Anova(control_group_BPM_con_mlm_analysis, type = "II", test.statistic = "F")
    #2.235 Loop through fixed effects in the model
    for (iv_name in row.names(control_group_BPM_con_ANCOVA_analysis)) {
      summary_lm <- summary(control_group_BPM_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #2.2351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          control_group_BPM_con_results_df <- rbind(
            control_group_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
        #2.2352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- control_group_BPM_con_ANCOVA_analysis[iv_name, "F"]
          df <- control_group_BPM_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- control_group_BPM_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- control_group_BPM_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          control_group_BPM_con_results_df <- rbind(
            control_group_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
#2.24 Subset the data based on the relevant IV variable strings for FDR correction of p values
control_group_BPM_con_bpm_p_adjust_subset <- subset(control_group_BPM_con_results_df, grepl("bpm_", IV))
#2.25 Create a new column conducting an FDR (p adjustment) on the derived p values
control_group_BPM_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(control_group_BPM_con_bpm_p_adjust_subset$p_value, method = "fdr")
#2.26 Create a new dataframe where only significant results are stored
control_group_BPM_con_significant_results <- subset(control_group_BPM_con_bpm_p_adjust_subset, p_adjusted <= 0.05)
#2.27 Join the FDR corrected p-values with the rest of the model results 
control_group_BPM_con_merged_adjusted_p_values <- left_join(control_group_BPM_con_results_df, control_group_BPM_con_bpm_p_adjust_subset)
#2.28 Remove all numbers from the strings in the column IV (for later merging purposes)
control_group_BPM_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", control_group_BPM_con_merged_adjusted_p_values$IV)
#2.291 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
control_group_BPM_con_merged_adjusted_p_values <- control_group_BPM_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))
#2.292 Remove rows with Intercept values
control_group_BPM_con_merged_adjusted_p_values <- control_group_BPM_con_merged_adjusted_p_values[control_group_BPM_con_merged_adjusted_p_values$IV != "(Intercept)", ]
#2.2101 Pivot the full model results to be wider
control_group_BPM_con_merged_adjusted_p_values_pivoted <- control_group_BPM_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
#2.2102 Reorder columns to match the order in your paper
control_group_BPM_con_merged_adjusted_p_values_pivoted <- control_group_BPM_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))
#2.211 Write the full model results as a csv file
write.csv(control_group_BPM_con_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/BPM_control_Group_Full_Model_Results.csv", row.names = FALSE)




#2.3 Test for differences in predicting connectivity based on BPM symptoms in the GAD group
#2.31 Create a columns range for the BPM columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14
#2.32 Create an empty dataframe to store analysis values
GAD_group_BPM_con_results_df <- data.frame(bpm_column_name = character(),
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
#2.33 Run the lm for each combination of bpm metric and connectivity metric
for (bpm_col_num in bpm_col_range) {
  #2.331 Get the bpm metric column name
  bpm_col_name <- colnames(GAD_Imaging_CBCL_BPM_Data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    #2.332 Get the connectivity metric column name
    connectivity_col_name <- colnames(GAD_Imaging_CBCL_BPM_Data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    #2.333 Run the linear regression model
    GAD_group_BPM_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, " + rsfmri_c_ngd_meanmotion + sex + (1|site_name)"), na.action = na.omit, data = GAD_Imaging_CBCL_BPM_Data)
    #2.334 Run the ANCOVA Model
    GAD_group_BPM_con_ANCOVA_analysis <- car::Anova(GAD_group_BPM_con_mlm_analysis, type = "II", test.statistic = "F")
    #2.335 Loop through fixed effects in the model
    for (iv_name in row.names(GAD_group_BPM_con_ANCOVA_analysis)) {
      summary_lm <- summary(GAD_group_BPM_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        #2.3351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          GAD_group_BPM_con_results_df <- rbind(
            GAD_group_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
        #2.3352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- GAD_group_BPM_con_ANCOVA_analysis[iv_name, "F"]
          df <- GAD_group_BPM_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- GAD_group_BPM_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- GAD_group_BPM_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          GAD_group_BPM_con_results_df <- rbind(
            GAD_group_BPM_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
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
#2.34 Subset the data based on the relevant IV variable strings for FDR correction of p values
GAD_group_BPM_con_bpm_p_adjust_subset <- subset(GAD_group_BPM_con_results_df, grepl("bpm_", IV))
#2.35 Create a new column conducting an FDR (p adjustment) on the derived p values
GAD_group_BPM_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(GAD_group_BPM_con_bpm_p_adjust_subset$p_value, method = "fdr")
#2.36 Create a new dataframe where only significant results are stored
GAD_group_BPM_con_significant_results <- subset(GAD_group_BPM_con_bpm_p_adjust_subset, p_adjusted <= 0.05)
#2.37 Join the FDR corrected p-values with the rest of the model results 
GAD_group_BPM_con_merged_adjusted_p_values <- left_join(GAD_group_BPM_con_results_df, GAD_group_BPM_con_bpm_p_adjust_subset)
#2.38 Remove all numbers from the strings in the column IV (for later merging purposes)
GAD_group_BPM_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", GAD_group_BPM_con_merged_adjusted_p_values$IV)
#2.391 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
GAD_group_BPM_con_merged_adjusted_p_values <- GAD_group_BPM_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))
#2.392 Remove rows with Intercept values
GAD_group_BPM_con_merged_adjusted_p_values <- GAD_group_BPM_con_merged_adjusted_p_values[GAD_group_BPM_con_merged_adjusted_p_values$IV != "(Intercept)", ]
#2.3101 Pivot the full model results to be wider
GAD_group_BPM_con_merged_adjusted_p_values_pivoted <- GAD_group_BPM_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")
#2.3102 Reorder columns to match the order in your paper
GAD_group_BPM_con_merged_adjusted_p_values_pivoted <- GAD_group_BPM_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))
#2.311 Write the full model results as a csv file
write.csv(GAD_group_BPM_con_merged_adjusted_p_values_pivoted, "/Users/samsievertsen/Desktop/SamResearch/Results/BPM_GAD_Group_Full_Model_Results.csv", row.names = FALSE)


