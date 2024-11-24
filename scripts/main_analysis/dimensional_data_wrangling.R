## Setup ##

# Load packages for loading, wrangling, mutating, and visualizing data and set environmental variables
library(dplyr)
options(digits = 6, scipen = 999)


# Read in GAD baseline and followup data
baseline_GAD <- read.csv("./data_raw/baseline_GAD.csv")
followup_GAD <- read.csv("./data_raw/followup_GAD.csv")
both_GAD_random_assignment <- read.csv("./data_raw/both_GAD_random_assignment.csv")

# Read in control sample data
control_sample <- read.csv("./data_raw/control_sample.csv")

# Read in rsfMRI imaging data
ABCD_rsfMRI_Data <- read.csv("./data_raw/ABCD_rsfMRI_Data.csv")
ABCD_rsfMRI_Data <- ABCD_rsfMRI_Data[-1,]

# Read in cbcl data
dimensional_analysis_cbcl_data <- read.csv("./data_raw/abcd_cbcls01.csv")
dimensional_analysis_cbcl_data <- dimensional_analysis_cbcl_data[-1, ]

# Read in site visit analysis data
dimensional_analysis_repeat_sample <- read.csv("./data_processed/main_analysis/site_visit_analysis_data.csv")

# Read in ABCD family ID data
abcd_family_id_data <- read.delim("./data_raw/acspsw03.txt")


## Data Wrangling ##

#1. cbcl 
#1.1 Create a dataframe containing GAD subjectkeys and their diagnosis status
dimensional_analysis_data <- full_join(baseline_GAD, followup_GAD)

#1.2 Select relevant columns from both GAD data and assign group label
both_GAD <- both_GAD_random_assignment %>%
  dplyr::select(-group)
both_GAD$group <- rep("GAD_Both")

#1.3 Add the "both GAD" group to the main dataset
dimensional_analysis_data <- full_join(dimensional_analysis_data, both_GAD)

#1.4 Remove unnecessary columns
dimensional_analysis_data <- dimensional_analysis_data %>%
  dplyr::select(-c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"))

#1.5 Write GAD subject keys to a CSV
write.csv(dimensional_analysis_data, "./data_processed/main_analysis/dimensional_analysis_gad_subjectkeys_groups.csv", row.names = FALSE)

#1.6 Extract and save control group data
dimensional_analysis_controls <- control_sample %>%
  dplyr::select(c(subjectkey, group))
write.csv(dimensional_analysis_controls, "./data_processed/main_analysis/dimensional_analysis_hc_subjectkeys_groups.csv", row.names = FALSE)

#1.7 Merge GAD and control subjectkeys into a single dataframe
dimensional_analysis_GAD_sub_groups <- read.csv("./data_processed/main_analysis/dimensional_analysis_GAD_subjectkeys_groups.csv")
dimensional_analysis_CN_sub_groups <- read.csv("./data_processed/main_analysis/dimensional_analysis_HC_subjectkeys_groups.csv")
dimensional_analysis_sample_groups <- full_join(dimensional_analysis_CN_sub_groups, dimensional_analysis_GAD_sub_groups)

#1.8 Merge the grouped subjectkeys with the rsfMRI imaging data
dimensional_analysis_grouped_imaging_data <- merge(dimensional_analysis_sample_groups, ABCD_rsfMRI_Data)

#1.9 Calculate age in years and keep relevant columns
dimensional_analysis_grouped_imaging_data$age_in_years <- floor(as.numeric(as.character(dimensional_analysis_grouped_imaging_data$interview_age)) / 12)
dimensional_analysis_grouped_imaging_data <- dimensional_analysis_grouped_imaging_data %>%
  dplyr::select(c(subjectkey, eventname, group, interview_age, age_in_years, sex, site_name,
                  rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))

#1.10 Convert numeric columns for rsfMRI data
dimensional_analysis_grouped_imaging_data <- dimensional_analysis_grouped_imaging_data %>%
  mutate_at(vars(8:13), as.numeric)

#1.11 Filter to include only subjects with both baseline and followup scans
dimensional_analysis_grouped_imaging_data <- dimensional_analysis_grouped_imaging_data %>%
  na.omit() %>%
  filter(rowSums(.) != "") %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") & any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()

#1.12 Clean and preprocess cbcl data
dimensional_analysis_cbcl_data <- dimensional_analysis_cbcl_data[-1, ]
dimensional_analysis_cbcl_data$interview_age <- as.numeric(dimensional_analysis_cbcl_data$interview_age)
dimensional_analysis_cbcl_data <- dimensional_analysis_cbcl_data %>%
  mutate_at(vars(10:90), as.numeric)

#1.13 Merge family ID data with analysis data
abcd_family_id_data <- abcd_family_id_data %>%
  filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::select(c(subjectkey, rel_family_id)) %>%
  rename(family_id = rel_family_id)
dimensional_analysis_repeat_sample <- left_join(dimensional_analysis_repeat_sample, abcd_family_id_data)

#1.14 Combine imaging, cbcl, and subject group data
cbcl_data <- left_join(dimensional_analysis_repeat_sample, dimensional_analysis_cbcl_data) %>%
  dplyr::select(c(subjectkey, eventname, family_id, analysis_group, group, interview_age,
                  age_in_years, sex, site_name, rsfmri_c_ngd_meanmotion,
                  rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta,
                  rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh,
                  cbcl_scr_dsm5_anxdisord_t, cbcl_scr_syn_anxdep_t, cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t))

#1.15 Convert factor and numeric columns as needed
cbcl_data <- cbcl_data %>%
  mutate(subjectkey = as.factor(subjectkey),
         family_id = as.factor(family_id),
         analysis_group = as.factor(analysis_group),
         eventname = as.factor(eventname),
         sex = as.factor(sex),
         site_name = as.factor(site_name),
         cbcl_scr_dsm5_anxdisord_t = as.numeric(cbcl_scr_dsm5_anxdisord_t),
         cbcl_scr_syn_anxdep_t = as.numeric(cbcl_scr_syn_anxdep_t),
         cbcl_scr_syn_internal_t = as.numeric(cbcl_scr_syn_internal_t),
         cbcl_scr_syn_external_t = as.numeric(cbcl_scr_syn_external_t))

#1.16 Exclude rows with NA or empty values in cbcl columns
#1.161 Create a cevtor containing the columns to check
cbcl_columns_to_check_for_empties <-
  c("cbcl_scr_dsm5_anxdisord_t",
    "cbcl_scr_syn_anxdep_t",
    "cbcl_scr_syn_internal_t",
    "cbcl_scr_syn_external_t")

#1.162 Check the columns of interest for 
for (col in cbcl_columns_to_check_for_empties) {
  cbcl_data <-
    cbcl_data[!(is.na(cbcl_data[[col]]) | cbcl_data[[col]] == ""),]
}

## Output ## 

#1. Write the CBCL data as a csv file
write.csv(cbcl_data, "./data_processed/main_analysis/dimensional_analysis_imaging_cbcl_data.csv", row.names = FALSE)