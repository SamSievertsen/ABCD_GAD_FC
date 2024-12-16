## Set Up ##

# Load in necessary packages and configure environmental variables
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(rcompanion)
library(forcats)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggseg)
library(ggseg3d)
library(ggsegGordon)
library(gridGraphics) 
library(patchwork)
options(digits = 8, scipen = 999)  

# Read in the mental health/medication treatment data for analysis
treatment_data <- read.csv("./data_processed/supplement/treatment_data.csv")

# Read in the grouped HC vs GAD connectivity difference data
site_visit_analysis_data <- read.csv("./data_processed/main_analysis/site_visit_analysis_data.csv")

# Read in the grouped repeated measures connectivity change/difference data
repeated_measures_analysis_data <- read.csv("./data_processed/main_analysis/repeated_measures_grouped_imaging_data.csv")


## Data Wrangling ## 

#1. Merge the treatment data with the HC + current GAD group connectivity difference analysis data
group_con_difference_treatment_data <- left_join(site_visit_analysis_data, treatment_data, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))


#2. Merge the treatment data with the repeated measures data
repeated_measures_treatment_data <- left_join(repeated_measures_analysis_data, treatment_data, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))

#3. Aggregate the repeated measures treatment data to check if treatment was received across timepoints
repeated_measures_aggregated_treatment_data <- repeated_measures_treatment_data %>%
  group_by(subjectkey, group) %>%
  summarize(ever_treated = any(received_any_treatment == 1), .groups = "drop") %>% 
  filter(group != "Control")

#4. Create a binary variable for newly initiated services and dropped services between baseline and 2 year follow up in the GAD group
repeated_measures_treatment_initiation_data <- repeated_measures_treatment_data %>% 
  dplyr::select(c(subjectkey, group, treatment_status)) %>% 
  group_by(subjectkey) %>% 
  mutate(new_services = ifelse(treatment_status == "tx_after_baseline", 1, 0),
         dropped_services = ifelse(treatment_status == "tx_before_baseline", 1, 0)) %>% 
  filter(group != "Control") %>% 
  slice(1)

#5. Create a version of the data for assessing whether change in treatment is associated with change in connectivity from baseline to 2Y followup
#5.1 Subset the columns and subjects of interest
repeated_measures_treatment_connectivity_change_data <- repeated_measures_treatment_data %>% 
  filter(group != "Control") %>% 
  dplyr::select(c(subjectkey, family_id, eventname, treatment_status, sex, site_name, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))

#5.2 Create a variable indicating whether subjects received treatment between baseline and follow up or not
repeated_measures_treatment_connectivity_change_data$treatment_bw_baseline_2y <- if_else(repeated_measures_treatment_connectivity_change_data$treatment_status == "baseline_and_followup_tx" | repeated_measures_treatment_connectivity_change_data$treatment_status == "tx_after_baseline", "tx_bw_baseline_2y", "no_tx_bw_baseline_2y")

#5.3 Ensure columns of interest are the correct data type for the analysis
#5.31 Factor type 
repeated_measures_treatment_connectivity_change_data <- repeated_measures_treatment_connectivity_change_data %>%
  mutate(
    subjectkey = as.factor(subjectkey),
    family_id = as.factor(family_id),
    eventname = as.factor(eventname),
    treatment_status = fct_relevel(as.factor(treatment_status), "no_tx_history"),
    treatment_bw_baseline_2y = fct_relevel(as.factor(treatment_bw_baseline_2y), "no_tx_bw_baseline_2y"),
    sex = as.factor(sex),
    site_name = as.factor(site_name))

#5.32 Numeric type
repeated_measures_treatment_connectivity_change_data <- mutate_at(repeated_measures_treatment_connectivity_change_data, vars(7:12), as.numeric)

#6. Create a whole sample dataframe with the 2 and 4 group treatment variables to assess whether GAD subgroup x Treatment subgroup interaction is associated with within-VAN connectivity
#6.1 Create a copy of the repeated measures treatment data
repeated_measures_within_VAN_treatment_analysis_data <- repeated_measures_treatment_data

#6.2 Create the 2 group treatment variable
repeated_measures_within_VAN_treatment_analysis_data <-
  repeated_measures_within_VAN_treatment_analysis_data %>%
  mutate(treatment_bw_baseline_2y = if_else(
      treatment_status == "baseline_and_followup_tx" |
        treatment_status == "tx_after_baseline", "tx_bw_baseline_2y", "no_tx_bw_baseline_2y"))

#6.3 Retain only columns + subjects of interest
repeated_measures_within_VAN_treatment_analysis_data <- repeated_measures_within_VAN_treatment_analysis_data %>%
  dplyr::select(c(subjectkey, family_id, eventname, sex, site_name, group, treatment_status, treatment_bw_baseline_2y, rsfmri_c_ngd_meanmotion, rsfmri_c_ngd_vta_ngd_vta)) %>% 
  filter(group != "Control")

#6.4 Ensure columns are the correct data type
repeated_measures_within_VAN_treatment_analysis_data <- repeated_measures_within_VAN_treatment_analysis_data %>%
  mutate(
    subjectkey = as.factor(subjectkey),
    family_id = as.factor(family_id),
    eventname = as.factor(eventname),
    group = as.factor(group),
    treatment_status = fct_relevel(as.factor(treatment_status), "no_tx_history"),
    treatment_bw_baseline_2y = fct_relevel(as.factor(treatment_bw_baseline_2y), "no_tx_bw_baseline_2y"),
    sex = as.factor(sex),
    site_name = as.factor(site_name),
    rsfmri_c_ngd_meanmotion = as.numeric(rsfmri_c_ngd_meanmotion),
    rsfmri_c_ngd_vta_ngd_vta = as.numeric(rsfmri_c_ngd_vta_ngd_vta))


## Data Analysis ## 

#1. Calculate the N and rate of each type of treatment for the GAD + HC group to be added to table 1
#1.11 Create a grouped descriptive stats summary of treatment
grouped_con_difference_treatment_summary_stats <- group_con_difference_treatment_data %>% 
  group_by(analysis_group) %>% 
  summarize(group_n = n(),
            n_outpatient_tx = sum(outpatient == 1, na.rm = TRUE),
            percent_outpatient_tx = (n_outpatient_tx / group_n) * 100,
            n_partial_hospitalization_tx = sum(partial_hospitalization == 1, na.rm = TRUE),
            percent_partial_hospitalization_tx = (n_partial_hospitalization_tx / group_n) * 100,
            n_inpatient_tx = sum(inpatient == 1, na.rm = TRUE),
            percent_inpatient_tx = (n_inpatient_tx / group_n) * 100,
            n_psychotherapy_tx = sum(psychotherapy == 1, na.rm = TRUE),
            percent_psychotherapy_tx = (n_psychotherapy_tx / group_n) * 100,
            n_medication_management_tx = sum(medication_management == 1, na.rm = TRUE),
            percent_medication_management_tx = (n_medication_management_tx / group_n) * 100,
            n_any_mental_health_tx = sum(any_mental_health_tx == 1, na.rm = TRUE),
            percent_any_mental_health_tx = (n_any_mental_health_tx / group_n) * 100,
            n_antidepressants_med = sum(taking_antidepressant == 1, na.rm = TRUE),
            percent_antidepressants_med = (n_antidepressants_med / group_n) * 100,
            n_anxiolytics_med = sum(taking_anxiolytic == 1, na.rm = TRUE),
            percent_anxiolytics_med = (n_anxiolytics_med / group_n) * 100,
            n_antidep_or_anxio_med = sum(taking_medication == 1, na.rm = TRUE),
            percent_antidep_or_anxio_med = (n_antidep_or_anxio_med / group_n) * 100,
            n_receiving_any_tx = sum(received_any_treatment == 1, na.rm = TRUE),
            percent_receiving_any_tx = (n_receiving_any_tx / group_n) * 100)

#1.12 Pivot the treatment summary stats into a long version of the table
grouped_con_difference_treatment_summary_stats_long <- grouped_con_difference_treatment_summary_stats %>%
  pivot_longer(
    cols = -analysis_group,
    names_to = "variable",      
    values_to = "value")

#1.13 Pivot the long version of the table back to a wider format, such that each diagnostic group gets their own row
grouped_con_difference_treatment_summary_stats_cleaned <- grouped_con_difference_treatment_summary_stats_long %>%
  pivot_wider(
    names_from = analysis_group,
    values_from = value)


#1.21 Create a whole sample descriptive stats summary of treatment
whole_sample_con_difference_treatment_summary_stats <- group_con_difference_treatment_data %>% 
  summarize(group_n = n(),
            n_outpatient_tx = sum(outpatient == 1, na.rm = TRUE),
            percent_outpatient_tx = (n_outpatient_tx / group_n) * 100,
            n_partial_hospitalization_tx = sum(partial_hospitalization == 1, na.rm = TRUE),
            percent_partial_hospitalization_tx = (n_partial_hospitalization_tx / group_n) * 100,
            n_inpatient_tx = sum(inpatient == 1, na.rm = TRUE),
            percent_inpatient_tx = (n_inpatient_tx / group_n) * 100,
            n_psychotherapy_tx = sum(psychotherapy == 1, na.rm = TRUE),
            percent_psychotherapy_tx = (n_psychotherapy_tx / group_n) * 100,
            n_medication_management_tx = sum(medication_management == 1, na.rm = TRUE),
            percent_medication_management_tx = (n_medication_management_tx / group_n) * 100,
            n_any_mental_health_tx = sum(any_mental_health_tx == 1, na.rm = TRUE),
            percent_any_mental_health_tx = (n_any_mental_health_tx / group_n) * 100,
            n_antidepressants_med = sum(taking_antidepressant == 1, na.rm = TRUE),
            percent_antidepressants_med = (n_antidepressants_med / group_n) * 100,
            n_anxiolytics_med = sum(taking_anxiolytic == 1, na.rm = TRUE),
            percent_anxiolytics_med = (n_anxiolytics_med / group_n) * 100,
            n_antidep_or_anxio_med = sum(taking_medication == 1, na.rm = TRUE),
            percent_antidep_or_anxio_med = (n_antidep_or_anxio_med / group_n) * 100,
            n_receiving_any_tx = sum(received_any_treatment == 1, na.rm = TRUE),
            percent_receiving_any_tx = (n_receiving_any_tx / group_n) * 100)

#1.22 Pivot the whole sample summary stats into a long version of the table
whole_sample_con_difference_treatment_summary_stats_long <- whole_sample_con_difference_treatment_summary_stats %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "whole_sample")


#1.3 Merge the grouped and whole sample treatment summary stats dataframes
group_con_difference_treatment_summary_stats <- left_join(grouped_con_difference_treatment_summary_stats_cleaned, whole_sample_con_difference_treatment_summary_stats_long)


#2. Calculate the N and rate of each type of treatment for the repeated measures subgroups to be added to the sample characteristics table for that analysis
#2.11 Create a grouped descriptive stats summary of treatment
grouped_repeated_measures_treatment_summary_stats <- repeated_measures_treatment_data %>% 
  group_by(group, eventname) %>% 
  summarize(group_n = n(),
            n_outpatient_tx = sum(outpatient == 1, na.rm = TRUE),
            percent_outpatient_tx = (n_outpatient_tx / group_n) * 100,
            n_partial_hospitalization_tx = sum(partial_hospitalization == 1, na.rm = TRUE),
            percent_partial_hospitalization_tx = (n_partial_hospitalization_tx / group_n) * 100,
            n_inpatient_tx = sum(inpatient == 1, na.rm = TRUE),
            percent_inpatient_tx = (n_inpatient_tx / group_n) * 100,
            n_psychotherapy_tx = sum(psychotherapy == 1, na.rm = TRUE),
            percent_psychotherapy_tx = (n_psychotherapy_tx / group_n) * 100,
            n_medication_management_tx = sum(medication_management == 1, na.rm = TRUE),
            percent_medication_management_tx = (n_medication_management_tx / group_n) * 100,
            n_any_mental_health_tx = sum(any_mental_health_tx == 1, na.rm = TRUE),
            percent_any_mental_health_tx = (n_any_mental_health_tx / group_n) * 100,
            n_antidepressants_med = sum(taking_antidepressant == 1, na.rm = TRUE),
            percent_antidepressants_med = (n_antidepressants_med / group_n) * 100,
            n_anxiolytics_med = sum(taking_anxiolytic == 1, na.rm = TRUE),
            percent_anxiolytics_med = (n_anxiolytics_med / group_n) * 100,
            n_antidep_or_anxio_med = sum(taking_medication == 1, na.rm = TRUE),
            percent_antidep_or_anxio_med = (n_antidep_or_anxio_med / group_n) * 100,
            n_receiving_any_tx = sum(received_any_treatment == 1, na.rm = TRUE),
            percent_receiving_any_tx = (n_receiving_any_tx / group_n) * 100,
            n_no_tx_history = sum(treatment_status == "no_tx_history", na.rm = TRUE),
            percent_no_tx_history = (n_no_tx_history / group_n) * 100,
            n_baseline_and_followup_tx = sum(treatment_status == "baseline_and_followup_tx", na.rm = TRUE),
            percent_baseline_and_followup_tx = (n_baseline_and_followup_tx / group_n) * 100,
            n_tx_before_baseline = sum(treatment_status == "tx_before_baseline", na.rm = TRUE),
            percent_tx_before_baseline = (n_tx_before_baseline / group_n) * 100,
            n_tx_after_baseline = sum(treatment_status == "tx_after_baseline", na.rm = TRUE),
            percent_tx_after_baseline = (n_tx_after_baseline / group_n) * 100,
            n_insufficient_data = sum(is.na(treatment_status)),
            percent_insufficient_data = (n_insufficient_data / group_n) * 100)

#2.12 Pivot the treatment summary stats into a long version of the table
grouped_repeated_measures_treatment_summary_stats_long <- grouped_repeated_measures_treatment_summary_stats %>%
  pivot_longer(
    cols = -c(group, eventname),
    names_to = "variable",      
    values_to = "value")

#2.13 Pivot the long version of the table back to a wider format, such that each diagnostic group gets their own row
grouped_repeated_measures_treatment_summary_stats_cleaned <- grouped_repeated_measures_treatment_summary_stats_long %>%
  pivot_wider(
    names_from = c(group, eventname),
    values_from = value)


#2.21 Create a whole sample descriptive stats summary of treatment
whole_sample_repeated_measures_treatment_summary_stats <- repeated_measures_treatment_data %>% 
  group_by(eventname) %>% 
  summarize(group_n = n(),
            n_outpatient_tx = sum(outpatient == 1, na.rm = TRUE),
            percent_outpatient_tx = (n_outpatient_tx / group_n) * 100,
            n_partial_hospitalization_tx = sum(partial_hospitalization == 1, na.rm = TRUE),
            percent_partial_hospitalization_tx = (n_partial_hospitalization_tx / group_n) * 100,
            n_inpatient_tx = sum(inpatient == 1, na.rm = TRUE),
            percent_inpatient_tx = (n_inpatient_tx / group_n) * 100,
            n_psychotherapy_tx = sum(psychotherapy == 1, na.rm = TRUE),
            percent_psychotherapy_tx = (n_psychotherapy_tx / group_n) * 100,
            n_medication_management_tx = sum(medication_management == 1, na.rm = TRUE),
            percent_medication_management_tx = (n_medication_management_tx / group_n) * 100,
            n_any_mental_health_tx = sum(any_mental_health_tx == 1, na.rm = TRUE),
            percent_any_mental_health_tx = (n_any_mental_health_tx / group_n) * 100,
            n_antidepressants_med = sum(taking_antidepressant == 1, na.rm = TRUE),
            percent_antidepressants_med = (n_antidepressants_med / group_n) * 100,
            n_anxiolytics_med = sum(taking_anxiolytic == 1, na.rm = TRUE),
            percent_anxiolytics_med = (n_anxiolytics_med / group_n) * 100,
            n_antidep_or_anxio_med = sum(taking_medication == 1, na.rm = TRUE),
            percent_antidep_or_anxio_med = (n_antidep_or_anxio_med / group_n) * 100,
            n_receiving_any_tx = sum(received_any_treatment == 1, na.rm = TRUE),
            percent_receiving_any_tx = (n_receiving_any_tx / group_n) * 100,
            n_no_tx_history = sum(treatment_status == "no_tx_history", na.rm = TRUE),
            percent_no_tx_history = (n_no_tx_history / group_n) * 100,
            n_baseline_and_followup_tx = sum(treatment_status == "baseline_and_followup_tx", na.rm = TRUE),
            percent_baseline_and_followup_tx = (n_baseline_and_followup_tx / group_n) * 100,
            n_tx_before_baseline = sum(treatment_status == "tx_before_baseline", na.rm = TRUE),
            percent_tx_before_baseline = (n_tx_before_baseline / group_n) * 100,
            n_tx_after_baseline = sum(treatment_status == "tx_after_baseline", na.rm = TRUE),
            percent_tx_after_baseline = (n_tx_after_baseline / group_n) * 100,
            n_insufficient_data = sum(is.na(treatment_status)),
            percent_insufficient_data = (n_insufficient_data / group_n) * 100)

#2.22 Pivot the whole sample summary stats into a long version of the table
whole_sample_repeated_measures_treatment_summary_stats_long <- whole_sample_repeated_measures_treatment_summary_stats %>%
  pivot_longer(
    cols = -eventname,
    names_to = "variable",
    values_to = "whole_sample")

#2.23 Pivot the long version of the table back to a wider format
whole_sample_repeated_measures_treatment_summary_stats_cleaned <- whole_sample_repeated_measures_treatment_summary_stats_long %>%
  pivot_wider(
    names_from = eventname,
    values_from = whole_sample)

#2.24 Change the column names to reflect that the table corresponds to the whole repeated measures sample
whole_sample_repeated_measures_treatment_summary_stats_cleaned <- whole_sample_repeated_measures_treatment_summary_stats_cleaned %>% 
  rename("whole_sample_2_year_follow_up_y_arm_1" = "2_year_follow_up_y_arm_1",
         "whole_sample_baseline_year_1_arm_1" = "baseline_year_1_arm_1")

#2.3 Merge the grouped and whole sample treatment summary stats dataframes
repeated_measures_treatment_summary_stats <- left_join(grouped_repeated_measures_treatment_summary_stats_cleaned, whole_sample_repeated_measures_treatment_summary_stats_cleaned)


#3. Analysis of differences in services received between the 3 GAD subgroups
#3.1 Create a contingency table for whether each subject had ever received treatment by 2 year follow up or not
repeated_measures_contingency_table <- table(repeated_measures_aggregated_treatment_data$group, repeated_measures_aggregated_treatment_data$ever_treated)

#3.2 Run a fishers exact test (given that some of the cells have N's < 5) to determine if the frequency of treatment by GAD subgroup is significantly different
repeated_measures_contingency_test <- fisher.test(repeated_measures_contingency_table)


#4. Determine whether the GAD subgroups differ in new service initiation or service between baseline and 2-year follow-up
#4.1 New treatment initiation
#4.11 Create a contingency table for the newly initiated services
repeated_measures_treatment_initiation_table <- table(repeated_measures_treatment_initiation_data$group, repeated_measures_treatment_initiation_data$new_services)

#4.12 Perform a Fisher's exact test to determine if the frequency of treatment initiation between baseline and 2 year follow up by GAD subgroup is significantly different
repeated_measures_treatment_initiation_test <- fisher.test(repeated_measures_treatment_initiation_table)

#4.13 Given that the results of the Fisher's exact test were significant, run a pairwise Fisher's on all individual levels of the treatment x GAD subgroup variables
pairwise_result <- pairwiseNominalIndependence(repeated_measures_treatment_initiation_table, method = "fdr")

#4.2 Dropping services after baseline
#4.21 Create a contingency table for the dropped treatment services
repeated_measures_dropped_treatment_table <- table(repeated_measures_treatment_initiation_data$group, repeated_measures_treatment_initiation_data$dropped_services)

#4.22 Perform a Fisher's exact test to determine if the frequency of dropped treatment between baseline and 2 year follow up by GAD subgroup is significantly different
repeated_measures_dropped_treatment_test <- fisher.test(repeated_measures_dropped_treatment_table)


#5. Assess whether connectivity change is associated with a change in treatment status in the 5 significant connectivity variables
#5.1 Establish the range of the dependent variables 
treatment_connectivity_analysis_dp_col_range <- 8:12

#5.2 Initialize an empty dataframe to store analysis values
treatment_connectivity_analysis_raw_results <- data.frame(
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

#5.3 Run the linear model for the site-visit group 
for (column_number in treatment_connectivity_analysis_dp_col_range) {
  
  #5.31 Get the column name
  dependent_variable <- colnames(repeated_measures_treatment_connectivity_change_data)[column_number]
  print(dependent_variable)
  
  #5.32 Fit the linear regression model
  treatment_connectivity_analysis_lm <- lmerTest::lmer(repeated_measures_treatment_connectivity_change_data[, dependent_variable] ~ treatment_bw_baseline_2y*eventname + rsfmri_c_ngd_meanmotion + sex+ (1|subjectkey) + (1|site_name) + (1|family_id), na.action = na.omit, data = repeated_measures_treatment_connectivity_change_data)
  
  #5.33 Run the ANCOVA Model on the linear regression to extract omnibus group effect(s)
  treatment_connectivity_analysis_ANCOVA <- car::Anova(treatment_connectivity_analysis_lm, type = "III", test.statistic = "F")
  
  #5.34 Loop through fixed effects in the model to extract relevant results
  for (independent_variable in row.names(treatment_connectivity_analysis_ANCOVA)) {
    if (independent_variable == "rsfmri_c_ngd_meanmotion") {
      
      #5.341 Extract + store parameters of interest for continuous variables
      summary_lm <- summary(treatment_connectivity_analysis_lm)
      estimate <- summary_lm$coefficients[independent_variable, "Estimate"]
      std_error <- summary_lm$coefficients[independent_variable, "Std. Error"]
      t_value <- summary_lm$coefficients[independent_variable, "t value"]
      p_value <- summary_lm$coefficients[independent_variable, "Pr(>|t|)"]
      treatment_connectivity_analysis_raw_results <- rbind(treatment_connectivity_analysis_raw_results, data.frame(
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
      
      #5.342 Extract + store parameters of interest for categorical variables
      f_value <- treatment_connectivity_analysis_ANCOVA[independent_variable, "F"]
      df <- treatment_connectivity_analysis_ANCOVA[independent_variable, "Df"]
      residual_df <- treatment_connectivity_analysis_ANCOVA[independent_variable, "Df.res"]
      p_value <- treatment_connectivity_analysis_ANCOVA[independent_variable, "Pr(>F)"] 
      treatment_connectivity_analysis_raw_results <- rbind(treatment_connectivity_analysis_raw_results, data.frame(
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

#5.41 Subset the data based on the relevant IV variable strings for FDR correction of p values
#5.411 Subset the data for interaction terms containing ":eventname" in the IV variable
treatment_connectivity_analysis_intx_p_adjust_subset <- subset(treatment_connectivity_analysis_raw_results, grepl(":eventname", IV))

#5.412 Subset the data for main effects containing only "treatment_bw_baseline_2y" in the IV variable
treatment_connectivity_analysis_p_adjust_subset <- subset(
  treatment_connectivity_analysis_raw_results,  grepl("^treatment_bw_baseline_2y$", IV))

#5.42 Create a new column conducting an FDR (p adjustment) on the derived p values
#5.421 Interaction term p values
treatment_connectivity_analysis_intx_p_adjust_subset$p_adjusted <- p.adjust(treatment_connectivity_analysis_intx_p_adjust_subset$p_value, method = "fdr")

#5.422 Main effect p values
treatment_connectivity_analysis_p_adjust_subset$p_adjusted <- p.adjust(treatment_connectivity_analysis_p_adjust_subset$p_value, method = "fdr")

#5.43 Merge the FDR corrected p-values together
treatment_connectivity_analysis_intx_p_adjust_subset <- full_join(treatment_connectivity_analysis_intx_p_adjust_subset, treatment_connectivity_analysis_p_adjust_subset)

#5.44 Join the FDR corrected p-values with the rest of the model results 
treatment_connectivity_analysis_merged_adjusted_p_values <- left_join(treatment_connectivity_analysis_raw_results, treatment_connectivity_analysis_intx_p_adjust_subset)

#5.45 Remove all numbers from the strings in the column IV (for later merging purposes)
treatment_connectivity_analysis_merged_adjusted_p_values$IV <- gsub("\\d+", "", treatment_connectivity_analysis_merged_adjusted_p_values$IV)

#5.461 Rename IV values to match the expected paper format
treatment_connectivity_analysis_merged_adjusted_p_values <- treatment_connectivity_analysis_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      IV == "treatment_bw_baseline_y" ~ "Treatment Between Baseline and 2Y",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint",
      IV == "treatment_bw_baseline_y:eventname" ~ "Treatment Between Baseline and 2Y*Timepoint Interaction",
      TRUE ~ IV))

#5.462 Remove rows with Intercept values
treatment_connectivity_analysis_merged_adjusted_p_values <- treatment_connectivity_analysis_merged_adjusted_p_values[treatment_connectivity_analysis_merged_adjusted_p_values$IV != "(Intercept)", ]

#5.47 Pivot the full model results to be in wide format
treatment_connectivity_analysis_merged_adjusted_p_values_pivoted <- treatment_connectivity_analysis_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = column_name,
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#5.48 Reorder columns to match the order in the paper
treatment_connectivity_analysis_merged_adjusted_p_values_pivoted <-
  treatment_connectivity_analysis_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(
    column_name,
    `Treatment Between Baseline and 2Y*Timepoint Interaction_f_value`,
    `Treatment Between Baseline and 2Y*Timepoint Interaction_df`,
    `Treatment Between Baseline and 2Y*Timepoint Interaction_residual_df`,
    `Treatment Between Baseline and 2Y*Timepoint Interaction_p_value`,
    `Treatment Between Baseline and 2Y*Timepoint Interaction_p_adjusted`,
    `Treatment Between Baseline and 2Y_f_value`,
    `Treatment Between Baseline and 2Y_df`,
    `Treatment Between Baseline and 2Y_residual_df`,
    `Treatment Between Baseline and 2Y_p_value`,
    `Treatment Between Baseline and 2Y_p_adjusted`, 
    `FD Motion_estimate`,
    `FD Motion_t_value`,
    `FD Motion_std_error`,
    `FD Motion_p_value`,
    `Sex_f_value`,
    `Sex_df`,
    `Sex_residual_df`,
    `Sex_p_value`,
    `Timepoint_f_value`,
    `Timepoint_df`,
    `Timepoint_residual_df`,
    `Timepoint_p_value`))



#6. Assess whether treatment subgroup interacting with GAD subgroup is significantly associated with a change in within-VAN connectivity
#6.1 Run the within-VAN GAD subgroup model with the 2 group treatment between baseline and 2Y variable as an interaction term
#6.11 Fit the repeated measures model for the 2 group tx variable
within_VAN_repeated_measures_2_group_model <- lmerTest::lmer(rsfmri_c_ngd_vta_ngd_vta ~ group*eventname*treatment_bw_baseline_2y + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_within_VAN_treatment_analysis_data)

#6.12 Apply a type III ANCOVA to the repeated measures model to assess for omnibus significance
within_VAN_repeated_measures_2_group_ANCOVA <- car::Anova(within_VAN_repeated_measures_2_group_model, type = "III", test.statistic="F")

#6.13 2 group model post hoc 
#6.131 Estimated marginal means within groups over time points
within_VAN_repeated_measures_2_group_EMM <- emmeans(within_VAN_repeated_measures_2_group_model, pairwise ~ eventname | group | treatment_bw_baseline_2y)

#6.132 Generate the emmeans contrast for comparing slopes between groups over time (AKA relationship between estimated means between timepoints)
within_VAN_repeated_measures_2_group_posthoc_contrasts <- contrast(within_VAN_repeated_measures_2_group_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL, adjust = "none")

#6.133 Generate standard error derived confidence intervals for the model results 
within_VAN_repeated_measures_2_group_posthoc_contrasts_CI <- confint(within_VAN_repeated_measures_2_group_posthoc_contrasts)


#6.2 Run the within-VAN GAD subgroup model with the 4 group treatment between baseline and 2Y variable as an interaction term
#6.21 Fit the repeated measures model for the 4 group tx variable
within_VAN_repeated_measures_4_group_model <- lmerTest::lmer(rsfmri_c_ngd_vta_ngd_vta ~ group*eventname*treatment_status + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_within_VAN_treatment_analysis_data)

#6.22 Apply a type III ANCOVA to the repeated measures model to assess for omnibus significance
within_VAN_repeated_measures_4_group_ANCOVA <- car::Anova(within_VAN_repeated_measures_4_group_model, type = "III", test.statistic="F")

#6.23 4 group model post hoc 
#6.231 Estimated marginal means within groups over time points
within_VAN_repeated_measures_4_group_EMM <- emmeans(within_VAN_repeated_measures_4_group_model, pairwise ~ eventname | group | treatment_status)

#6.232 Generate the emmeans contrast for comparing slopes between groups over time (AKA relationship between estimated means between timepoints)
within_VAN_repeated_measures_4_group_posthoc_contrasts <- contrast(within_VAN_repeated_measures_4_group_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL, adjust = "none")

#6.233 Generate standard error derived confidence intervals for the model results 
within_VAN_repeated_measures_4_group_posthoc_contrasts_CI <- confint(within_VAN_repeated_measures_4_group_posthoc_contrasts)


## Output ## 

#1. Write the group connectivity difference treatment summary stats as a csv file
write.csv(group_con_difference_treatment_summary_stats, "./results/hc_gad_group_con_difference_treatment_summary_stats.csv", row.names = FALSE)

#2. Write the repeated measures treatment summary stats as a csv file
write.csv(repeated_measures_treatment_summary_stats, "./results/repeated_measures_treatment_summary_stats.csv", row.names = FALSE)

#3. Write the connectivity ~ treatment subgroup repeated measures analysis results as a csv file
write.csv(treatment_connectivity_analysis_merged_adjusted_p_values_pivoted, "./results/treatment_group_connectivity_analysis_results.csv", row.names = FALSE)