## Setup ## 

# Load packages for loading, wrangling, mutating, and visualizing data
library(dplyr)
library(writexl)
library(readxl)
library(stringr)
options(digits = 9)

# Read in the control sample data
control_sample <- read.csv("./data_processed/main_analysis/resampled_hc_sample.csv")

# Read in the GAD group data
GAD_sample <- read.csv("./data_processed/main_analysis/gad_group.csv")

# Read in the cleaned, QC'd + QA'd imaging data
clinical_imaging_data <- read.csv("./data_processed/main_analysis/subset_qcd_imaging_data.csv")

# Read in the Family ID data
abcd_family_id_data <- read.delim("./data_raw/acspsw03.txt") %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::select(c(subjectkey, rel_family_id))


## Data Wrangling ##

#1. Merge the dataframe containing GAD cases with matched controls (based on site-visit distributions)
#1.1 Remove irrelevant columns from the control sample dataframe
control_sample <- control_sample %>% dplyr::select(-X)

#1.2 Merge the control and imaging dataframes
control_sample <- merge(control_sample, clinical_imaging_data, 
                        by = c("subjectkey", "site_name", "eventname"), 
                        all.x = TRUE) %>%
  dplyr::select(c(subjectkey, eventname, site_name, group, sex, interview_age))

#1.3 Calculate the age of participants in years
control_sample$age_in_years <- floor((as.numeric(control_sample$interview_age)) / 12)

#1.4 Merge the control and GAD samples into a single dataframe
case_control_sample <- merge(control_sample, GAD_sample, all = TRUE)

#1.5 Create a column to classify subjects into analysis groups. Baseline GAD and follow-up GAD subjects are classified as "GAD", and control subjects remain classified as "control"
case_control_sample$analysis_group <- if_else(case_control_sample$group == "control", "control", "GAD")

#2. Create a cleaned dataframe containing relevant DVs and IVs for site-visit GAD & Control subjects
#2.1 Merge case-control sample with clinical imaging data, excluding the column `ksads_10_869_p`
site_visit_analysis_data <- case_control_sample %>%
  dplyr::select(subjectkey, eventname, group, sex, interview_age, site_name, age_in_years, analysis_group) %>%
  dplyr::left_join(
    dplyr::select(clinical_imaging_data, -ksads_10_869_p),
    by = c("subjectkey","eventname", "sex", "interview_age", "site_name"))

#2.2 Merge the site-visit analysis data with family ID data
site_visit_analysis_data <- left_join(site_visit_analysis_data, abcd_family_id_data)

#2.3 Ensure column ordering is optimized for downstream analyses
site_visit_analysis_data <- site_visit_analysis_data %>% 
  dplyr::select(subjectkey, site_name, eventname, interview_age, sex, group, age_in_years, analysis_group, imgincl_rsfmri_include, rel_family_id, rsfmri_c_ngd_meanmotion, everything())


## Output ## 

#1. Write the site-visit analysis data as a csv file for use in downstream analyses
write.csv(site_visit_analysis_data, "./data_processed/main_analysis/site_visit_analysis_data.csv", row.names = FALSE)
