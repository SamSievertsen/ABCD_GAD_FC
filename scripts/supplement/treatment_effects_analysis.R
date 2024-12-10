## Set Up ##

# Load in necessary packages and configure environmental variables
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
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
repeated_measures_treatment_summary_stats <- left_join(grouped_con_difference_treatment_summary_stats_cleaned, whole_sample_con_difference_treatment_summary_stats_long)


## Output ## 

#1. Write the group connectivity difference treatment summary stats as a csv file
write.csv(group_con_difference_treatment_summary_stats, "./results/hc_gad_group_con_difference_treatment_summary_stats.csv", row.names = FALSE)

#2. Write the repeated measures treatment summary stats as a csv file
write.csv(repeated_measures_treatment_summary_stats, "./results/repeated_measures_treatment_summary_stats.csv", row.names = FALSE)
