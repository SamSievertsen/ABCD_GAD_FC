## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(sampler)
library(sampling)
options(digits = 8, scipen = 999) 

# Read in required data 
# Healthy control subjects
healthy_control_group <- read.csv("./data_processed/healty_control_group.csv")

# GAD subjects
GAD_group <- read.csv("./data_processed/GAD_group.csv")

# Resting state fMRI data
rsfMRI_data <- read.csv("./data_raw/ABCD_rsfMRI_Data.csv")

# Recommended imaging data for QC/QA 
recommended_imaging <- read.csv("./data_raw/abcd_imgincl01.csv")

# Site Data
site_data <- read.csv("./data_raw/ABCD_sitename_data.csv") %>% 
  dplyr::select(c(subjectkey, eventname, interview_age, sex, site_name))


## Data Cleaning ## 

#1. Imaging Data
#1.11 Remove the non-data rows from the rsfMRI data
rsfMRI_data <- rsfMRI_data[-1, ]

#1.12 Remove the non-data rows from the recommended imaging data
recommended_imaging <- recommended_imaging[-1, ]

#1.21 Retain only the columns of interest in the rsfMRI data
rsfMRI_data <- rsfMRI_data %>% 
  dplyr::select(c(src_subject_id, interview_age, sex, eventname, rsfmri_c_ngd_meanmotion, rsfmri_c_ngd_cgc_ngd_cgc, rsfmri_c_ngd_cgc_ngd_dt, rsfmri_c_ngd_cgc_ngd_fo, rsfmri_c_ngd_cgc_ngd_vta, rsfmri_c_ngd_dla_ngd_cgc, rsfmri_c_ngd_dla_ngd_dla, rsfmri_c_ngd_dla_ngd_dt, rsfmri_c_ngd_dla_ngd_sa, rsfmri_c_ngd_dla_ngd_vta, rsfmri_c_ngd_dt_ngd_dt, rsfmri_c_ngd_dt_ngd_fo, rsfmri_c_ngd_dt_ngd_sa, rsfmri_c_ngd_dt_ngd_vta, rsfmri_c_ngd_fo_ngd_dt, rsfmri_c_ngd_fo_ngd_fo, rsfmri_c_ngd_fo_ngd_sa, rsfmri_c_ngd_fo_ngd_vta, rsfmri_c_ngd_sa_ngd_cgc, rsfmri_c_ngd_sa_ngd_sa, rsfmri_c_ngd_sa_ngd_vta, rsfmri_c_ngd_vta_ngd_dt, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_cerc_scs_aalh, rsfmri_cor_ngd_cerc_scs_aarh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_cor_ngd_cerc_scs_agrh, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_cderh, rsfmri_cor_ngd_cerc_scs_hplh, rsfmri_cor_ngd_cerc_scs_hprh, rsfmri_cor_ngd_cerc_scs_ptlh, rsfmri_cor_ngd_cerc_scs_ptrh, rsfmri_cor_ngd_cerc_scs_thplh, rsfmri_cor_ngd_cerc_scs_thprh, rsfmri_cor_ngd_df_scs_aalh, rsfmri_cor_ngd_df_scs_aarh, rsfmri_cor_ngd_df_scs_aglh, rsfmri_cor_ngd_df_scs_agrh, rsfmri_cor_ngd_df_scs_cdelh, rsfmri_cor_ngd_df_scs_cderh, rsfmri_cor_ngd_df_scs_hplh, rsfmri_cor_ngd_df_scs_hprh, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_df_scs_ptrh, rsfmri_cor_ngd_df_scs_thplh, rsfmri_cor_ngd_df_scs_thprh, rsfmri_cor_ngd_dsa_scs_aalh, rsfmri_cor_ngd_dsa_scs_aarh, rsfmri_cor_ngd_dsa_scs_aglh, rsfmri_cor_ngd_dsa_scs_agrh, rsfmri_cor_ngd_dsa_scs_cdelh, rsfmri_cor_ngd_dsa_scs_cderh, rsfmri_cor_ngd_dsa_scs_hplh, rsfmri_cor_ngd_dsa_scs_hprh, rsfmri_cor_ngd_dsa_scs_ptlh, rsfmri_cor_ngd_dsa_scs_ptrh, rsfmri_cor_ngd_dsa_scs_thplh, rsfmri_cor_ngd_dsa_scs_thprh, rsfmri_cor_ngd_fopa_scs_aalh, rsfmri_cor_ngd_fopa_scs_aarh, rsfmri_cor_ngd_fopa_scs_aglh, rsfmri_cor_ngd_fopa_scs_agrh, rsfmri_cor_ngd_fopa_scs_cdelh, rsfmri_cor_ngd_fopa_scs_cderh, rsfmri_cor_ngd_fopa_scs_hplh, rsfmri_cor_ngd_fopa_scs_hprh, rsfmri_cor_ngd_fopa_scs_ptlh, rsfmri_cor_ngd_fopa_scs_ptrh, rsfmri_cor_ngd_fopa_scs_thplh, rsfmri_cor_ngd_fopa_scs_thprh, rsfmri_cor_ngd_sa_scs_aalh, rsfmri_cor_ngd_sa_scs_aarh, rsfmri_cor_ngd_sa_scs_aglh, rsfmri_cor_ngd_sa_scs_agrh, rsfmri_cor_ngd_sa_scs_cdelh, rsfmri_cor_ngd_sa_scs_cderh, rsfmri_cor_ngd_sa_scs_hplh, rsfmri_cor_ngd_sa_scs_hprh, rsfmri_cor_ngd_sa_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptrh, rsfmri_cor_ngd_sa_scs_thplh, rsfmri_cor_ngd_sa_scs_thprh, rsfmri_cor_ngd_vta_scs_aalh, rsfmri_cor_ngd_vta_scs_aarh, rsfmri_cor_ngd_vta_scs_aglh, rsfmri_cor_ngd_vta_scs_agrh, rsfmri_cor_ngd_vta_scs_cdelh, rsfmri_cor_ngd_vta_scs_cderh, rsfmri_cor_ngd_vta_scs_hplh, rsfmri_cor_ngd_vta_scs_hprh, rsfmri_cor_ngd_vta_scs_ptlh, rsfmri_cor_ngd_vta_scs_ptrh, rsfmri_cor_ngd_vta_scs_thplh, rsfmri_cor_ngd_vta_scs_thprh))

#1.22 Retain only the columns of interest in the recommended imaging data
recommended_imaging <- recommended_imaging %>% 
  dplyr::select(c(src_subject_id, eventname, imgincl_rsfmri_include))

#1.31 Change all rsfMRI columns to numeric type 
rsfMRI_data <- rsfMRI_data %>%
  mutate_at(vars(5:99), as.numeric)

#1.32 Change the recommended imaging column to numeric type
recommended_imaging$imgincl_rsfmri_include <- as.numeric(recommended_imaging$imgincl_rsfmri_include)

#1.33 Change the interview age column in the imaging + site data to numeric type
rsfMRI_data$interview_age <- as.numeric(rsfMRI_data$interview_age)
site_data$interview_age <- as.numeric(site_data$interview_age)

#1.34 Change the subjectkey column in the site data to src_subject_id for merging purposes
colnames(site_data)[colnames(site_data) == "subjectkey"] <- "src_subject_id"

#1.41 Merge the recommended imaging and rsfMRI data 
imaging_data_qcd <- full_join(rsfMRI_data, recommended_imaging)

#1.42 Merge the site name data with the imaging data 
imaging_data_qcd <- left_join(imaging_data_qcd, site_data)

#1.5 Only retain subjects with complete, quality controlled data
#1.51 Remove any rows not containing values in the rsfMRI colmns of interest
imaging_data_qcd <- imaging_data_qcd %>%
  filter(across(5:99, ~ !is.na(.x) & .x != ""))

#1.52 Only retain subject rows whos data passed QC/QA
imaging_data_qcd <- imaging_data_qcd %>% 
  filter(imgincl_rsfmri_include == 1)

#1.6 Retain only columns of interest
imaging_data_qcd <- imaging_data_qcd %>% 
  dplyr::select(c(src_subject_id, interview_age, sex, eventname, site_name, everything())) %>% 
  dplyr::select(-imgincl_rsfmri_include)


## Data Wrangling + Harmonization ##

#1. Control Group
#1.1 Create a merged version of the control group + rsfMRI data
healthy_control_rsfMRI_data <- inner_join(healthy_control_group, imaging_data_qcd)

#2. GAD Group
#2.1 Create a merged version of the GAD group + rsfMRI data
GAD_group_rsfMRI_data <- left_join(GAD_group, imaging_data_qcd)

#2.2 Re-adjust randomly sampled continuous GAD subjects based upon the availability of their data
#2.21 Create a function to check for missing rsfMRI data in columns 10 to 104
is_complete_data <- function(row_data) {
  complete <- all(!is.na(row_data[10:104]) & row_data[10:104] != "")
  return(complete)
}

# 2.22 Filter GAD_timepoint = "both" and iterate through these subjects
GAD_group_rsfMRI_data_cleaned <- GAD_group_rsfMRI_data %>%
  group_by(src_subject_id) %>%
  mutate(flagged_for_review = FALSE) %>%
  group_modify(~ {
    
    # 2.221 If GAD_timepoint is not "both", return the original data
    if (!all(.x$GAD_timepoint == "both")) {
      return(.x)
    }
    
    # 2.222 For "both" GAD subjects, extract rows for baseline and followup
    baseline_row <- .x %>% filter(eventname == "baseline_year_1_arm_1")
    followup_row <- .x %>% filter(eventname == "2_year_follow_up_y_arm_1")
    
    # 2.223 Determine which event the subject was assigned to in the analysis_group
    assigned_event <- ifelse(.x$analysis_group[1] == "baseline", "baseline_year_1_arm_1", "2_year_follow_up_y_arm_1")
    
    # 2.224 Check if both baseline and follow-up rows have complete rsfMRI data
    if (is_complete_data(baseline_row) && is_complete_data(followup_row)) {
      # If both have complete data, select one based on a criterion (e.g., random or assigned event)
      if (assigned_event == "baseline_year_1_arm_1") {
        return(baseline_row)  # If assigned to baseline, keep the baseline row
      } else {
        return(followup_row)  # If assigned to followup, keep the followup row
      }
    }
    
    # 2.225 Check if the assigned event row has complete data
    if (assigned_event == "baseline_year_1_arm_1" && is_complete_data(baseline_row)) {
      return(baseline_row)  # Keep baseline row, drop followup row
    } else if (assigned_event == "2_year_follow_up_y_arm_1" && is_complete_data(followup_row)) {
      return(followup_row)  # Keep followup row, drop baseline row
    }
    
    # 2.226 If the assigned event row does not have complete data, check the other row
    if (assigned_event == "baseline_year_1_arm_1" && is_complete_data(followup_row)) {
      followup_row <- followup_row %>% mutate(analysis_group = "followup")  # Change assignment to followup
      return(followup_row)  # Keep followup row, drop baseline row
    } else if (assigned_event == "2_year_follow_up_y_arm_1" && is_complete_data(baseline_row)) {
      baseline_row <- baseline_row %>% mutate(analysis_group = "baseline")  # Change assignment to baseline
      return(baseline_row)  # Keep baseline row, drop followup row
    }
    
    # 2.227 If neither row has complete data, log for review
    .x <- .x %>% mutate(flagged_for_review = TRUE)
    return(.x)
  }) %>%
  ungroup()

# 2.3 Flag rows with incomplete data for non-"both" subjects
GAD_group_rsfMRI_data_cleaned <- GAD_group_rsfMRI_data_cleaned %>%
  rowwise() %>% 
  mutate(flagged_for_review = ifelse(
    GAD_timepoint != "both" & !is_complete_data(cur_data()),
    TRUE, flagged_for_review)) %>%
  ungroup() 

# 2.4 Extract the rows that are flagged for review into a separate dataframe
flagged_subjects <- GAD_group_rsfMRI_data_cleaned %>%
  filter(flagged_for_review == TRUE) %>%
  dplyr::select(src_subject_id, eventname)

# 2.5 Remove flagged rows and the flagged column
GAD_group_rsfMRI_data_cleaned <- GAD_group_rsfMRI_data_cleaned %>%
  filter(flagged_for_review == FALSE) %>%
  select(-flagged_for_review)

#3. Merged Data
#3.1 Full join the healthy control and GAD group data
merged_GAD_HC_sample_data <- full_join(GAD_group_rsfMRI_data_cleaned, healthy_control_rsfMRI_data)

#3.2 Add HC values to the analysis_group variable (+ retain current values for GAD group)
merged_GAD_HC_sample_data <- merged_GAD_HC_sample_data %>% 
  mutate(analysis_group = if_else(group == "HC", "HC", analysis_group))

#3.3 Retain only columns of interest
merged_GAD_HC_sample_data_cleaned <- merged_GAD_HC_sample_data %>% 
  dplyr::select(c(src_subject_id, interview_age, sex, site_name, eventname, group, GAD_timepoint, analysis_group, subgroup))


## Output ## 

#1. Write the unsubset imaging data as a csv for later use
write.csv(imaging_data_qcd, "./data_processed/rsfMRI_data_qcd_unsubset.csv", row.names = FALSE)

#2. Write the merged GAD + HC group data as a csv for use in resampling 
write.csv(merged_GAD_HC_sample_data_cleaned, "./data_processed/GAD_HC_subjects_not_resampled.csv", row.names = FALSE)
