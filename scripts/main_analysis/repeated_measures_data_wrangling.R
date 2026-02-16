## Setup ##

# Load packages for loading, wrangling, mutating, and visualizing data + set environmental variables
library(dplyr)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
options(digits = 6, scipen = 999)

# Load in resting state fMRI data
qcd_rsfmri_data <- read.csv("./data_processed/main_analysis/subset_qcd_imaging_data.csv")

# Load in merged HC + GAD subject sample & exclude imaging data
merged_hc_gad_sample <- read.csv("./data_processed/main_analysis/site_visit_analysis_data.csv") %>% 
  dplyr::select(c(subjectkey, site_name, rel_family_id, sex, group, analysis_group))

## Data Wrangling ## 

#1. Clean and merge the sample data with the imaging data
#1.1 Change the name of the family ID column
merged_hc_gad_sample <- merged_hc_gad_sample %>% rename(family_id = rel_family_id)

#1.2 Join the subjects of interest with the cleaned + QC'd/QA's rsfMRI imaging data
repeated_measures_grouped_imaging_data <- left_join(merged_hc_gad_sample, qcd_rsfmri_data, by = c("subjectkey" = "src_subject_id", "sex" = "sex", "site_name" = "site_name"))

#1.3 Keep only the vars of interest 
#1.31 Create the age in years variable
repeated_measures_grouped_imaging_data$age_in_years <- floor(as.numeric(as.character(repeated_measures_grouped_imaging_data$interview_age)) / 12)

#1.32 Retain columns of interest
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>% 
  dplyr::select(c(subjectkey, eventname, group, interview_age, age_in_years, sex, site_name, family_id, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))

#1.4 Change the values of the "group" variable to represent their experimental identity
#1.41 HC
repeated_measures_grouped_imaging_data$group[repeated_measures_grouped_imaging_data$group == "control"] <- "Control"

#1.42 GAD Remitters
repeated_measures_grouped_imaging_data$group[repeated_measures_grouped_imaging_data$group == "baseline_GAD"] <- "GAD Remitter"

#1.43 GAD Converters
repeated_measures_grouped_imaging_data$group[repeated_measures_grouped_imaging_data$group == "followup_GAD"] <- "GAD Converter"

#1.44 Continuous GAD
repeated_measures_grouped_imaging_data$group[repeated_measures_grouped_imaging_data$group == "GAD_Both"] <- "Continuous GAD"


#1.5 Remove any rows with NA or empty values and retain subjects who have both baseline and followup scans only 
#1.51 Remove rows with NA's
repeated_measures_grouped_imaging_data <- na.omit(repeated_measures_grouped_imaging_data)

#1.52 Remove rows with empty values
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data[rowSums(repeated_measures_grouped_imaging_data == "") == 0, ]

#1.53 Only retain rows with compelte data at both baseline and follow up assessment timepoints
repeated_measures_grouped_imaging_data <-  repeated_measures_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") &
           any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()


## Output ##

#1. Write the complete, merged clinical with imaging data as a csv for use in repeated measures analysis
write.csv(repeated_measures_grouped_imaging_data, "./data_processed/main_analysis/repeated_measures_grouped_imaging_data.csv", row.names = FALSE)
