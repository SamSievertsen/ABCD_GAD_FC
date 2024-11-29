## Set-Up ## 

# Load packages for loading, wrangling, mutating, and visualizing data
library(dplyr)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(stats)
library(lme4)
library(sampling)
library(rsample)
library(splitstackshape)
options(digits = 9, scipen = 999)

# Load in raw clinical data
KSADS_Data <- read.csv("./data_raw/abcd_ksad_data_amended_dx.csv")

# Remove the first row of the KSADS df to account for the descriptions of the variables in the pertinent diagnoses and raw data
KSADS_Data_no_desc <- KSADS_Data[-1,]

# Load in the rsfMRI data
ABCD_rsfMRI_Data <- read.csv("./data_raw/ABCD_rsfMRI_Data.csv")

# Remove the first row of the df to account for the descriptions of the variables in the pertinent diagnoses and raw data
ABCD_rsfMRI_Data_no_desc <- ABCD_rsfMRI_Data[-1,]

# Create the site_name variable containing the scanner and site info
ABCD_rsfMRI_Data_no_desc$site_name <- substr(ABCD_rsfMRI_Data_no_desc$rsfmri_c_ngd_visitid, 1, 4)

# Load in super healthy control subject data to create the control vector
control_subjects_all_timepoints <- read_xlsx("./data_processed/main_analysis/ABCD_control_subjects_all_timepoints.xlsx")

# Create a vector containing the super healthy control ID's; these will later be merged based on 
sh_controls_vector <- control_subjects_all_timepoints$subjectkey

# Load in the recommended imaging inclusion data
recommended_imaging <- read.csv("./data_raw/abcd_imgincl01.csv")

# Remove the first row of the recommended imaging df to account for the descriptions of the variables in the pertinent columns
recommended_imaging <- recommended_imaging[-1,]


## Data Wrangling ## 

#1.11 Create a cleaned version of the KSADS data containing vars of interest
KSADS_vars_with_NA <- subset(
  KSADS_Data_no_desc,
  select = c("subjectkey", "eventname",
    na.omit("ksads_10_869_p"), "interview_age", "sex"))

KSADS_cleaned <-
  KSADS_vars_with_NA[complete.cases(KSADS_vars_with_NA[, "ksads_10_869_p"]), ]

#1.12 Create a cleaned version of the rsfMRI data containing vars of interest
rsfMRI_cleaned <- ABCD_rsfMRI_Data_no_desc[,c("subjectkey", "eventname", "site_name", "rsfmri_c_ngd_meanmotion", "rsfmri_c_ngd_cgc_ngd_cgc", "rsfmri_c_ngd_cgc_ngd_dt", "rsfmri_c_ngd_cgc_ngd_fo", "rsfmri_c_ngd_cgc_ngd_vta", "rsfmri_c_ngd_dla_ngd_cgc", "rsfmri_c_ngd_dla_ngd_dla", "rsfmri_c_ngd_dla_ngd_dt", "rsfmri_c_ngd_dla_ngd_sa", "rsfmri_c_ngd_dla_ngd_vta", "rsfmri_c_ngd_dt_ngd_dt", "rsfmri_c_ngd_dt_ngd_fo", "rsfmri_c_ngd_dt_ngd_sa", "rsfmri_c_ngd_dt_ngd_vta", "rsfmri_c_ngd_fo_ngd_dt", "rsfmri_c_ngd_fo_ngd_fo", "rsfmri_c_ngd_fo_ngd_sa", "rsfmri_c_ngd_fo_ngd_vta", "rsfmri_c_ngd_sa_ngd_cgc", "rsfmri_c_ngd_sa_ngd_sa", "rsfmri_c_ngd_sa_ngd_vta", "rsfmri_c_ngd_vta_ngd_dt", "rsfmri_c_ngd_vta_ngd_vta", "rsfmri_cor_ngd_cerc_scs_aalh", "rsfmri_cor_ngd_cerc_scs_aarh", "rsfmri_cor_ngd_cerc_scs_aglh", "rsfmri_cor_ngd_cerc_scs_agrh", "rsfmri_cor_ngd_cerc_scs_cdelh", "rsfmri_cor_ngd_cerc_scs_cderh", "rsfmri_cor_ngd_cerc_scs_hplh", "rsfmri_cor_ngd_cerc_scs_hprh", "rsfmri_cor_ngd_cerc_scs_ptlh", "rsfmri_cor_ngd_cerc_scs_ptrh", "rsfmri_cor_ngd_cerc_scs_thplh", "rsfmri_cor_ngd_cerc_scs_thprh", "rsfmri_cor_ngd_df_scs_aalh", "rsfmri_cor_ngd_df_scs_aarh", "rsfmri_cor_ngd_df_scs_aglh", "rsfmri_cor_ngd_df_scs_agrh", "rsfmri_cor_ngd_df_scs_cdelh", "rsfmri_cor_ngd_df_scs_cderh", "rsfmri_cor_ngd_df_scs_hplh", "rsfmri_cor_ngd_df_scs_hprh", "rsfmri_cor_ngd_df_scs_ptlh", "rsfmri_cor_ngd_df_scs_ptrh", "rsfmri_cor_ngd_df_scs_thplh", "rsfmri_cor_ngd_df_scs_thprh", "rsfmri_cor_ngd_dsa_scs_aalh", "rsfmri_cor_ngd_dsa_scs_aarh", "rsfmri_cor_ngd_dsa_scs_aglh", "rsfmri_cor_ngd_dsa_scs_agrh", "rsfmri_cor_ngd_dsa_scs_cdelh", "rsfmri_cor_ngd_dsa_scs_cderh", "rsfmri_cor_ngd_dsa_scs_hplh", "rsfmri_cor_ngd_dsa_scs_hprh", "rsfmri_cor_ngd_dsa_scs_ptlh", "rsfmri_cor_ngd_dsa_scs_ptrh", "rsfmri_cor_ngd_dsa_scs_thplh", "rsfmri_cor_ngd_dsa_scs_thprh", "rsfmri_cor_ngd_fopa_scs_aalh", "rsfmri_cor_ngd_fopa_scs_aarh", "rsfmri_cor_ngd_fopa_scs_aglh", "rsfmri_cor_ngd_fopa_scs_agrh", "rsfmri_cor_ngd_fopa_scs_cdelh", "rsfmri_cor_ngd_fopa_scs_cderh", "rsfmri_cor_ngd_fopa_scs_hplh", "rsfmri_cor_ngd_fopa_scs_hprh", "rsfmri_cor_ngd_fopa_scs_ptlh", "rsfmri_cor_ngd_fopa_scs_ptrh", "rsfmri_cor_ngd_fopa_scs_thplh", "rsfmri_cor_ngd_fopa_scs_thprh", "rsfmri_cor_ngd_sa_scs_aalh", "rsfmri_cor_ngd_sa_scs_aarh", "rsfmri_cor_ngd_sa_scs_aglh", "rsfmri_cor_ngd_sa_scs_agrh", "rsfmri_cor_ngd_sa_scs_cdelh", "rsfmri_cor_ngd_sa_scs_cderh", "rsfmri_cor_ngd_sa_scs_hplh", "rsfmri_cor_ngd_sa_scs_hprh", "rsfmri_cor_ngd_sa_scs_ptlh", "rsfmri_cor_ngd_sa_scs_ptrh", "rsfmri_cor_ngd_sa_scs_thplh", "rsfmri_cor_ngd_sa_scs_thprh", "rsfmri_cor_ngd_vta_scs_aalh", "rsfmri_cor_ngd_vta_scs_aarh", "rsfmri_cor_ngd_vta_scs_aglh", "rsfmri_cor_ngd_vta_scs_agrh", "rsfmri_cor_ngd_vta_scs_cdelh", "rsfmri_cor_ngd_vta_scs_cderh", "rsfmri_cor_ngd_vta_scs_hplh", "rsfmri_cor_ngd_vta_scs_hprh", "rsfmri_cor_ngd_vta_scs_ptlh", "rsfmri_cor_ngd_vta_scs_ptrh", "rsfmri_cor_ngd_vta_scs_thplh", "rsfmri_cor_ngd_vta_scs_thprh")] 

#1.2 Merge the KSADS data with the rsfMRI data using the rsfMRI data as the reference (to retain only baseline and 2 year follow up data)
merged_data_NA <- merge(rsfMRI_cleaned, KSADS_cleaned, by = c("subjectkey","eventname"), all.x = TRUE, suffixes = c("_KSADS","_rsfMRI"))

#1.31 Remove any rows with NA values (no imaging data) for analysis purposes
imaging_and_clinical_data_no_NA <- na.omit(merged_data_NA)

#1.32 Remove rows with empty strings for the same reason 
clinical_imaging_data_no_empty_NA <- imaging_and_clinical_data_no_NA[rowSums(imaging_and_clinical_data_no_NA == "") == 0, ]

#1.331 Remove rows from the data that contain rsfMRI scan data not passing QC. First, merge the two data frames
merged_recommended_clinical_imaging_data <- merge(clinical_imaging_data_no_empty_NA, recommended_imaging, by = c("subjectkey", "eventname", "interview_age", "sex"))

#1.332 Filter out rows that did not pass rsfMRI QC and remove unneccessary columns
clinical_imaging_data <- subset(merged_recommended_clinical_imaging_data, imgincl_rsfmri_include == 1)

clinical_imaging_data <- clinical_imaging_data %>% 
  dplyr::select(-c(collection_id, abcd_imgincl01_id, dataset_id, src_subject_id, interview_date, visit, imgincl_t1w_include, imgincl_t2w_include, imgincl_dmri_include, imgincl_rsfmri_include, imgincl_mid_include, imgincl_nback_include, imgincl_sst_include, collection_title))


#2. Create a merged data frame with subject IDs and their diagnosis of generalized anxiety disorder at baseline and 2-year follow-up
wide_clinical_data <- clinical_imaging_data %>%
  filter(ksads_10_869_p == 1) %>%
  dplyr::select(subjectkey, eventname, ksads_10_869_p) %>%
  pivot_wider(names_from = eventname, values_from = ksads_10_869_p)

#2.1 Save the data columns as numeric to avoid issues related to numbers savedas characters
wide_clinical_data$baseline_year_1_arm_1 <- as.numeric(wide_clinical_data$baseline_year_1_arm_1)
wide_clinical_data$`2_year_follow_up_y_arm_1` <- as.numeric(wide_clinical_data$`2_year_follow_up_y_arm_1`)

#2.2 Replace NA values with a "999" value to allow the use of filtering arguments in the next steps
wide_clinical_data$baseline_year_1_arm_1 <- replace(as.numeric(wide_clinical_data$baseline_year_1_arm_1), is.na(as.numeric(wide_clinical_data$baseline_year_1_arm_1)), 999)
wide_clinical_data$`2_year_follow_up_y_arm_1` <- replace(as.numeric(wide_clinical_data$`2_year_follow_up_y_arm_1`), is.na(as.numeric(wide_clinical_data$`2_year_follow_up_y_arm_1`)), 999)

#3.11 Create a group of subjects who have a positive diagnosis of generalized anxiety disorder at baseline, but not at 2-year follow-up
baseline_GAD <- wide_clinical_data %>%
  filter(`baseline_year_1_arm_1` == 1 & `2_year_follow_up_y_arm_1` != 1) %>%
  mutate(group = "baseline_GAD")

#3.12 Create a group of subjects who have a positive diagnosis of generalized anxiety disorder at 2-year followup, but not at baseline
followup_GAD <- wide_clinical_data %>%
  filter(`2_year_follow_up_y_arm_1` == 1 & `baseline_year_1_arm_1` != 1) %>%
  mutate(group = "followup_GAD")

#3.2 Identify subjects with positive diagnosis of GAD at both time points and randomly assign them to either group
both_GAD_random_assignment <- wide_clinical_data %>%
  filter(baseline_year_1_arm_1 == 1 & `2_year_follow_up_y_arm_1` == 1) %>%
  mutate(group = case_when(runif(n()) > 0.5 ~ "baseline_GAD", TRUE ~ "followup_GAD"))

#3.3 Bind the three data frames together
GAD_groups <- bind_rows(baseline_GAD, followup_GAD, both_GAD_random_assignment) 

#3.41 Add in the relevant demographic data for determining distribution of controls
GAD_grouped_demo_info_raw <- GAD_groups %>%
  left_join(clinical_imaging_data %>%
      dplyr::select(subjectkey, sex, interview_age, site_name, eventname), by = "subjectkey") %>%
  mutate(interview_age = as.numeric(interview_age),
    age = case_when(
      group == "baseline_GAD" &
        eventname == "baseline_year_1_arm_1" ~ interview_age,
      group == "followup_GAD" &
        eventname == "2_year_follow_up_y_arm_1" ~ interview_age,
      TRUE ~ NA_real_)) %>%
  filter(!is.na(age)) 

#3.421 Create a cleaned version of the newly merged demo/dx data
GAD_grouped_demo <- GAD_grouped_demo_info_raw %>% 
  dplyr::select(subjectkey, sex, interview_age, site_name, group, eventname)

#3.422 Create a new column for age in years to get more descriptive distributions
GAD_grouped_demo$age_in_years <- floor((as.numeric(GAD_grouped_demo$interview_age)) / 12)

#3.43 Create a new dataframe containing the relevant information for the control subjects 
raw_control_data <- clinical_imaging_data %>%
  filter(subjectkey %in% sh_controls_vector) %>%
  dplyr::select(subjectkey, sex, site_name, eventname, interview_age)

#3.44 Ensure data is converted to the correct type
raw_control_data$age_in_years <- floor((as.numeric(raw_control_data$interview_age)) / 12)

#3.51 Determine if there are any subjects who have assessment data at separate sites
multi_site_subjects <- c("NDAR_INV6D02DABX", "NDAR_INV7YJ6UFDB", "NDAR_INVA2ZPZTJJ", "NDAR_INVNTAR3TAF", "NDAR_INVPP5NK3LX", "NDAR_INVUFF64VGJ", "NDAR_INVUMJVH4W6", "NDAR_INVWF7C1DEL", "NDAR_INVY92TEZW6")

#3.52 Store IDs of subjects with assessment data at multiple sites
multi_site_subjects <- c("NDAR_INV6D02DABX", "NDAR_INV7YJ6UFDB", "NDAR_INVA2ZPZTJJ", "NDAR_INVNTAR3TAF", "NDAR_INVPP5NK3LX", "NDAR_INVUFF64VGJ", "NDAR_INVUMJVH4W6", "NDAR_INVWF7C1DEL", "NDAR_INVY92TEZW6")

#3.532 Randomly select one data point for each subject with assessment data at multiple sites. First, initialize a new data frame to store the selected data points
selected_data <- data.frame()

# Loop through the multi-site subjects
for (subject in multi_site_subjects) {
  
  # Subset the data for the current subject
  subject_data <- subset(raw_control_data, subjectkey == subject)
  
  # Randomly select either the baseline or followup data point with a 50% chance for each
  if (runif(1) > 0.5) {
    selected_data <- rbind(selected_data, subject_data[subject_data$eventname == "baseline_year_1_arm_1", ])
  } else {
    selected_data <- rbind(selected_data, subject_data[subject_data$eventname == "2_year_follow_up_y_arm_1", ])
  }
}

#3.54 Remove the original data points for the multi-site subjects from the raw control data
raw_control_data <- subset(raw_control_data, !(subjectkey %in% multi_site_subjects))

#3.55 Append the selected data points to the raw control data
raw_control_data <- rbind(raw_control_data, selected_data)


#4.1 Summarize the count of GAD subjects for each site_name and eventname group and calculate the percentage of visits to each site/visit combo attributed to GAD. Them perform the same on the raw_control_data df and create a new df  that contains the count and percentage of visits to each site attributed to the control group (by site/event)
site_visit_GAD_count <- GAD_grouped_demo %>% 
  group_by(site_name, eventname) %>%
  summarize(GAD_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(GAD_site_n = sum(GAD_site_visit_n),
         GAD_site_pct = (GAD_site_visit_n / GAD_site_n))

#4.11 Round the site percentage to 2 digits for thresholding resampling accuracy to a viable level
site_visit_GAD_count$GAD_site_pct <- round(site_visit_GAD_count$GAD_site_pct, digits = 2)

#4.12 Generate the current control counts/percentages
site_visit_control_count <- raw_control_data %>% 
  group_by(site_name, eventname) %>%
  summarize(control_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(control_site_n = sum(control_site_visit_n),
         control_site_pct = (control_site_visit_n / control_site_n))

#4.2 Merge the pct distributions at each site-visit combination in the GAD data
site_visit_merged_pct <- merge(site_visit_GAD_count, site_visit_control_count, by = c("eventname", "site_name"), all.x =  TRUE)
site_visit_merged_pct <- site_visit_merged_pct %>% 
  group_by(eventname, site_name) %>%
  mutate(max_controls = floor(GAD_site_pct * control_site_n))


## G010 Resampling ##

#5.111 Filter the control data for just Site G010
G010_control_data <- raw_control_data %>%
  filter(site_name == "G010") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.112 Combine subject data at one site with n and pct values for that site
G010_merged_pct <- merge(site_visit_merged_pct, G010_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.113 Create a df for just baseline and a df for just followup data 
G010_control_baseline_only <- G010_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)

G010_control_followup_only <- G010_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.114 save number of baseline only subjects and followup only subjects
G010_baseline_only_sub_n <- as.numeric(nrow(G010_control_baseline_only))
G010_followup_only_sub_n <- as.numeric(nrow(G010_control_followup_only))

#5.115 save the subject IDs for the baseline only subjects and followup only subjects
G010_baseline_control_subs_only <- G010_control_baseline_only$subjectkey
G010_followup_control_subs_only <- G010_control_followup_only$subjectkey

#5.116 save number of unique control subjects at this site
G010_num_unique_control_subs <- as.numeric(length(unique(G010_merged_pct$subjectkey)))

#5.117 save subsetted dataframe for only G010 GAD subjects
G010_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G010") 

#5.118 calculate total number of subjects at site G010
G010_total_GAD_sub_n <- as.numeric(length(unique(G010_GAD_subs$subjectkey)))

#5.119 further subset the dataframe to only contain baseline scans or followup GAD scans
G010_baseline_GAD_data <- G010_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
G010_followup_GAD_data <- G010_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1110 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
G010_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G010" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
G010_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G010" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1111 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
G010_both_timepoint_sub_list <- G010_control_data$subjectkey[duplicated(G010_control_data$subjectkey) ]
G010_both_timepoint_sub_list_n <- as.numeric(length(unique(G010_both_timepoint_sub_list)))

#5.1112 Subset a dataframe of those subjects
G010_both_timepoint_subs <- G010_control_data[G010_control_data$subjectkey %in% G010_both_timepoint_sub_list, ]

#5.1113 Create vectors for the max number of baseline and max number of followup scans 
G010_max_baseline_controls <- as.numeric(G010_baseline_only_sub_n + G010_both_timepoint_sub_list_n)
G010_max_followup_controls <- as.numeric(G010_followup_only_sub_n + G010_both_timepoint_sub_list_n)

#5.1114 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
G010_baseline_controls_desired <- as.numeric(G010_num_unique_control_subs * G010_baseline_GAD_site_visit_percent)
G010_followup_controls_desired <- as.numeric(G010_num_unique_control_subs * G010_followup_GAD_site_visit_percent)

#5.1115 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
G010_sufficient_baseline_control_n <- if_else(G010_max_baseline_controls > G010_baseline_controls_desired, TRUE, FALSE)
G010_sufficient_followup_control_n <- if_else(G010_max_followup_controls > G010_followup_controls_desired, TRUE, FALSE)

#5.1116 Create an adjusted total control N to sample in the event either of the above checks are failed 
G010_adjusted_total_control_n_if_insufficient_baseline <- floor(G010_max_baseline_controls / G010_baseline_GAD_site_visit_percent)
G010_adjusted_total_control_n_if_insufficient_followup <- floor(G010_max_followup_controls / G010_followup_GAD_site_visit_percent)

#5.1117 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
G010_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1118 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!G010_sufficient_baseline_control_n &
    G010_sufficient_followup_control_n) {
  G010_updated_control_counts$num_baseline_controls <-
    G010_max_baseline_controls
  G010_updated_control_counts$num_followup_controls <-
    floor(
      G010_adjusted_total_control_n_if_insufficient_baseline * G010_followup_GAD_site_visit_percent
    )
} else if (G010_sufficient_baseline_control_n &
           !G010_sufficient_followup_control_n) {
  G010_updated_control_counts$num_baseline_controls <-
    floor(
      G010_adjusted_total_control_n_if_insufficient_followup * G010_baseline_GAD_site_visit_percent
    )
  G010_updated_control_counts$num_followup_controls <-
    G010_max_followup_controls
} else {
  G010_updated_control_counts$num_baseline_controls <-
    G010_baseline_controls_desired
  G010_updated_control_counts$num_followup_controls <-
    G010_followup_controls_desired
}

#5.1119 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
G010_baseline_control_sample <- G010_control_data[G010_control_data$subjectkey %in% G010_baseline_control_subs_only, ]
G010_followup_control_sample <- G010_control_data[G010_control_data$subjectkey %in% G010_followup_control_subs_only, ] 

#5.1120 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
G010_num_baseline_controls_to_sample <- as.numeric(G010_updated_control_counts$num_baseline_controls)
G010_current_baseline_control_n <- as.numeric(length(unique(G010_baseline_control_sample$subjectkey)))

#5.1121 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (G010_current_baseline_control_n == G010_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G010_current_baseline_control_n > G010_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    G010_baseline_control_sample <-
      G010_baseline_control_sample[sample(nrow(G010_baseline_control_sample),
                                          G010_num_baseline_controls_to_sample),]
  } else if (G010_current_baseline_control_n < G010_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    G010_baseline_control_sample_to_add <-
      G010_control_data[G010_control_data$subjectkey %in% G010_both_timepoint_sub_list &
                          G010_control_data$eventname == "baseline_year_1_arm_1",]
    G010_baseline_control_sample_to_add <-
      G010_baseline_control_sample_to_add[sample(
        nrow(G010_baseline_control_sample_to_add),
        G010_num_baseline_controls_to_sample - G010_current_baseline_control_n),]
    G010_baseline_control_sample <-
      rbind(G010_baseline_control_sample,
            G010_baseline_control_sample_to_add)
  }
}

#5.1122 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
G010_num_followup_controls_to_sample <- as.numeric(G010_updated_control_counts$num_followup_controls)
G010_current_followup_control_n <- as.numeric(length(unique(G010_followup_control_sample$subjectkey)))
G010_subs_already_sampled <- G010_baseline_control_sample$subjectkey

#5.1123 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (G010_current_followup_control_n == G010_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G010_current_followup_control_n > G010_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    G010_followup_control_sample <-
      G010_followup_control_sample[sample(nrow(G010_followup_control_sample),
                                          G010_num_followup_controls_to_sample),]
  } else if (G010_current_followup_control_n < G010_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    G010_followup_control_sample_to_add <-
      G010_control_data[G010_control_data$subjectkey %in% G010_both_timepoint_sub_list &
                          !G010_control_data$subjectkey %in% G010_subs_already_sampled &
                          G010_control_data$eventname == "2_year_follow_up_y_arm_1",]
    G010_followup_control_sample_to_add <-
      G010_followup_control_sample_to_add[sample(
        nrow(G010_followup_control_sample_to_add),
        G010_num_followup_controls_to_sample - G010_current_followup_control_n),]
    G010_followup_control_sample <-
      rbind(G010_followup_control_sample,
            G010_followup_control_sample_to_add)
  }
}

#5.1124 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G010_baseline_control_sample$group <- rep("control")
G010_followup_control_sample$group <- rep("control")

#5.1125 Merge the baseline and followup G010 samples 
G010_control_sample <- merge(G010_baseline_control_sample, G010_followup_control_sample, all = TRUE)

#5.1126 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G010_cn_sample <- ifelse(length(unique(G010_control_sample$subjectkey)) == length(G010_control_sample$subjectkey), TRUE, FALSE)

#5.1127 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- G010_control_sample


## G031 Resampling ##

#5.121 Filter the control data for just Site G031
G031_control_data <- raw_control_data %>%
  filter(site_name == "G031") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.122 Combine subject data at one site with n and pct values for that site
G031_merged_pct <- merge(site_visit_merged_pct, G031_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.123 Create a df for just baseline and a df for just followup data 
G031_control_baseline_only <- G031_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
G031_control_followup_only <- G031_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.124 save number of baseline only subjects and followup only subjects
G031_baseline_only_sub_n <- as.numeric(nrow(G031_control_baseline_only))
G031_followup_only_sub_n <- as.numeric(nrow(G031_control_followup_only))

#5.125 save the subject IDs for the baseline only subjects and followup only subjects
G031_baseline_control_subs_only <- G031_control_baseline_only$subjectkey
G031_followup_control_subs_only <- G031_control_followup_only$subjectkey

#5.126 save number of unique control subjects at this site
G031_num_unique_control_subs <- as.numeric(length(unique(G031_merged_pct$subjectkey)))

#5.127 save subsetted dataframe for only G031 GAD subjects
G031_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G031") 

#5.128 calculate total number of subjects at site G031
G031_total_GAD_sub_n <- as.numeric(length(unique(G031_GAD_subs$subjectkey)))

#5.129 further subset the dataframe to only contain baseline scans or followup GAD scans
G031_baseline_GAD_data <- G031_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
G031_followup_GAD_data <- G031_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1210 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
G031_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G031" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
G031_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G031" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1211 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
G031_both_timepoint_sub_list <- G031_control_data$subjectkey[duplicated(G031_control_data$subjectkey) ]
G031_both_timepoint_sub_list_n <- as.numeric(length(unique(G031_both_timepoint_sub_list)))

#5.1212 Subset a dataframe of those subjects
G031_both_timepoint_subs <- G031_control_data[G031_control_data$subjectkey %in% G031_both_timepoint_sub_list, ]

#5.1213 Create vectors for the max number of baseline and max number of followup scans 
G031_max_baseline_controls <- as.numeric(G031_baseline_only_sub_n + G031_both_timepoint_sub_list_n)
G031_max_followup_controls <- as.numeric(G031_followup_only_sub_n + G031_both_timepoint_sub_list_n)

#5.1214 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
G031_baseline_controls_desired <- as.numeric(G031_num_unique_control_subs * G031_baseline_GAD_site_visit_percent)
G031_followup_controls_desired <- as.numeric(G031_num_unique_control_subs * G031_followup_GAD_site_visit_percent)

#5.1215 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
G031_sufficient_baseline_control_n <- if_else(G031_max_baseline_controls > G031_baseline_controls_desired, TRUE, FALSE)
G031_sufficient_followup_control_n <- if_else(G031_max_followup_controls > G031_followup_controls_desired, TRUE, FALSE)

#5.1216 Create an adjusted total control N to sample in the event either of the above checks are failed 
G031_adjusted_total_control_n_if_insufficient_baseline <- floor(G031_max_baseline_controls / G031_baseline_GAD_site_visit_percent)
G031_adjusted_total_control_n_if_insufficient_followup <- floor(G031_max_followup_controls / G031_followup_GAD_site_visit_percent)

#5.1217 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
G031_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1218 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!G031_sufficient_baseline_control_n &
    G031_sufficient_followup_control_n) {
  G031_updated_control_counts$num_baseline_controls <-
    G031_max_baseline_controls
  G031_updated_control_counts$num_followup_controls <-
    floor(
      G031_adjusted_total_control_n_if_insufficient_baseline * G031_followup_GAD_site_visit_percent
    )
} else if (G031_sufficient_baseline_control_n &
           !G031_sufficient_followup_control_n) {
  G031_updated_control_counts$num_baseline_controls <-
    floor(
      G031_adjusted_total_control_n_if_insufficient_followup * G031_baseline_GAD_site_visit_percent
    )
  G031_updated_control_counts$num_followup_controls <-
    G031_max_followup_controls
} else {
  G031_updated_control_counts$num_baseline_controls <-
    G031_baseline_controls_desired
  G031_updated_control_counts$num_followup_controls <-
    G031_followup_controls_desired
}

#5.1219 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
G031_baseline_control_sample <- G031_control_data[G031_control_data$subjectkey %in% G031_baseline_control_subs_only, ]
G031_followup_control_sample <- G031_control_data[G031_control_data$subjectkey %in% G031_followup_control_subs_only, ] 

#5.1220 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
G031_num_baseline_controls_to_sample <- as.numeric(G031_updated_control_counts$num_baseline_controls)
G031_current_baseline_control_n <- as.numeric(length(unique(G031_baseline_control_sample$subjectkey)))

#5.1221 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (G031_current_baseline_control_n == G031_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G031_current_baseline_control_n > G031_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    G031_baseline_control_sample <-
      G031_baseline_control_sample[sample(nrow(G031_baseline_control_sample),
                                          G031_num_baseline_controls_to_sample),]
  } else if (G031_current_baseline_control_n < G031_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    G031_baseline_control_sample_to_add <-
      G031_control_data[G031_control_data$subjectkey %in% G031_both_timepoint_sub_list &
                          G031_control_data$eventname == "baseline_year_1_arm_1",]
    G031_baseline_control_sample_to_add <-
      G031_baseline_control_sample_to_add[sample(
        nrow(G031_baseline_control_sample_to_add),
        G031_num_baseline_controls_to_sample - G031_current_baseline_control_n),]
    G031_baseline_control_sample <-
      rbind(G031_baseline_control_sample,
            G031_baseline_control_sample_to_add)
  }
}

#5.1222 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
G031_num_followup_controls_to_sample <- as.numeric(G031_updated_control_counts$num_followup_controls)
G031_current_followup_control_n <- as.numeric(length(unique(G031_followup_control_sample$subjectkey)))
G031_subs_already_sampled <- G031_baseline_control_sample$subjectkey

#5.1223 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (G031_current_followup_control_n == G031_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G031_current_followup_control_n > G031_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    G031_followup_control_sample <-
      G031_followup_control_sample[sample(nrow(G031_followup_control_sample),
                                          G031_num_followup_controls_to_sample),]
  } else if (G031_current_followup_control_n < G031_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    G031_followup_control_sample_to_add <-
      G031_control_data[G031_control_data$subjectkey %in% G031_both_timepoint_sub_list &
                          !G031_control_data$subjectkey %in% G031_subs_already_sampled &
                          G031_control_data$eventname == "2_year_follow_up_y_arm_1",]
    G031_followup_control_sample_to_add <-
      G031_followup_control_sample_to_add[sample(
        nrow(G031_followup_control_sample_to_add),
        G031_num_followup_controls_to_sample - G031_current_followup_control_n
      ),]
    G031_followup_control_sample <-
      rbind(G031_followup_control_sample,
            G031_followup_control_sample_to_add)
  }
}

#5.1224 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G031_baseline_control_sample$group <- rep("control")
G031_followup_control_sample$group <- rep("control")

#5.1225 Merge the baseline and followup G031 samples 
G031_control_sample <- merge(G031_baseline_control_sample, G031_followup_control_sample, all = TRUE)

#5.1226 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G031_cn_sample <- ifelse(length(unique(G031_control_sample$subjectkey)) == length(G031_control_sample$subjectkey), TRUE, FALSE)

#5.1227 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, G031_control_sample, all = TRUE)


## G032 Resampling ##

#5.131 Filter the control data for just Site G032
G032_control_data <- raw_control_data %>%
  filter(site_name == "G032") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.132 Combine subject data at one site with n and pct values for that site
G032_merged_pct <- merge(site_visit_merged_pct, G032_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.133 Create a df for just baseline and a df for just followup data 
G032_control_baseline_only <- G032_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
G032_control_followup_only <- G032_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.134 save number of baseline only subjects and followup only subjects
G032_baseline_only_sub_n <- as.numeric(nrow(G032_control_baseline_only))
G032_followup_only_sub_n <- as.numeric(nrow(G032_control_followup_only))

#5.135 save the subject IDs for the baseline only subjects and followup only subjects
G032_baseline_control_subs_only <- G032_control_baseline_only$subjectkey
G032_followup_control_subs_only <- G032_control_followup_only$subjectkey

#5.136 save number of unique control subjects at this site
G032_num_unique_control_subs <- as.numeric(length(unique(G032_merged_pct$subjectkey)))

#5.137 save subsetted dataframe for only G032 GAD subjects
G032_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G032") 

#5.138 calculate total number of subjects at site G032
G032_total_GAD_sub_n <- as.numeric(length(unique(G032_GAD_subs$subjectkey)))

#5.139 further subset the dataframe to only contain baseline scans or followup GAD scans
G032_baseline_GAD_data <- G032_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
G032_followup_GAD_data <- G032_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1310 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
G032_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G032" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
G032_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G032" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1311 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
G032_both_timepoint_sub_list <- G032_control_data$subjectkey[duplicated(G032_control_data$subjectkey) ]
G032_both_timepoint_sub_list_n <- as.numeric(length(unique(G032_both_timepoint_sub_list)))

#5.1312 Subset a dataframe of those subjects
G032_both_timepoint_subs <- G032_control_data[G032_control_data$subjectkey %in% G032_both_timepoint_sub_list, ]

#5.1313 Create vectors for the max number of baseline and max number of followup scans 
G032_max_baseline_controls <- as.numeric(G032_baseline_only_sub_n + G032_both_timepoint_sub_list_n)
G032_max_followup_controls <- as.numeric(G032_followup_only_sub_n + G032_both_timepoint_sub_list_n)

#5.1314 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
G032_baseline_controls_desired <- as.numeric(G032_num_unique_control_subs * G032_baseline_GAD_site_visit_percent)
G032_followup_controls_desired <- as.numeric(G032_num_unique_control_subs * G032_followup_GAD_site_visit_percent)

#5.1315 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
G032_sufficient_baseline_control_n <- if_else(G032_max_baseline_controls > G032_baseline_controls_desired, TRUE, FALSE)
G032_sufficient_followup_control_n <- if_else(G032_max_followup_controls > G032_followup_controls_desired, TRUE, FALSE)

#5.1316 Create an adjusted total control N to sample in the event either of the above checks are failed 
G032_adjusted_total_control_n_if_insufficient_baseline <- floor(G032_max_baseline_controls / G032_baseline_GAD_site_visit_percent)
G032_adjusted_total_control_n_if_insufficient_followup <- floor(G032_max_followup_controls / G032_followup_GAD_site_visit_percent)

#5.1317 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
G032_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1318 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!G032_sufficient_baseline_control_n &
    G032_sufficient_followup_control_n) {
  G032_updated_control_counts$num_baseline_controls <-
    G032_max_baseline_controls
  G032_updated_control_counts$num_followup_controls <-
    floor(
      G032_adjusted_total_control_n_if_insufficient_baseline * G032_followup_GAD_site_visit_percent
    )
} else if (G032_sufficient_baseline_control_n &
           !G032_sufficient_followup_control_n) {
  G032_updated_control_counts$num_baseline_controls <-
    floor(
      G032_adjusted_total_control_n_if_insufficient_followup * G032_baseline_GAD_site_visit_percent
    )
  G032_updated_control_counts$num_followup_controls <-
    G032_max_followup_controls
} else {
  G032_updated_control_counts$num_baseline_controls <-
    G032_baseline_controls_desired
  G032_updated_control_counts$num_followup_controls <-
    G032_followup_controls_desired
}

#5.1319 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
G032_baseline_control_sample <- G032_control_data[G032_control_data$subjectkey %in% G032_baseline_control_subs_only, ]
G032_followup_control_sample <- G032_control_data[G032_control_data$subjectkey %in% G032_followup_control_subs_only, ] 

#5.1320 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
G032_num_baseline_controls_to_sample <- as.numeric(G032_updated_control_counts$num_baseline_controls)
G032_current_baseline_control_n <- as.numeric(length(unique(G032_baseline_control_sample$subjectkey)))

#5.1321 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (G032_current_baseline_control_n == G032_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G032_current_baseline_control_n > G032_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    G032_baseline_control_sample <-
      G032_baseline_control_sample[sample(nrow(G032_baseline_control_sample),
                                          G032_num_baseline_controls_to_sample),]
  } else if (G032_current_baseline_control_n < G032_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    G032_baseline_control_sample_to_add <-
      G032_control_data[G032_control_data$subjectkey %in% G032_both_timepoint_sub_list &
                          G032_control_data$eventname == "baseline_year_1_arm_1",]
    G032_baseline_control_sample_to_add <-
      G032_baseline_control_sample_to_add[sample(
        nrow(G032_baseline_control_sample_to_add),
        G032_num_baseline_controls_to_sample - G032_current_baseline_control_n),]
    G032_baseline_control_sample <-
      rbind(G032_baseline_control_sample,
            G032_baseline_control_sample_to_add)
  }
}

#5.1322 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
G032_num_followup_controls_to_sample <- as.numeric(G032_updated_control_counts$num_followup_controls)
G032_current_followup_control_n <- as.numeric(length(unique(G032_followup_control_sample$subjectkey)))
G032_subs_already_sampled <- G032_baseline_control_sample$subjectkey

#5.1323 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (G032_current_followup_control_n == G032_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G032_current_followup_control_n > G032_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    G032_followup_control_sample <-
      G032_followup_control_sample[sample(nrow(G032_followup_control_sample),
                                          G032_num_followup_controls_to_sample),]
  } else if (G032_current_followup_control_n < G032_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    G032_followup_control_sample_to_add <-
      G032_control_data[G032_control_data$subjectkey %in% G032_both_timepoint_sub_list &
                          !G032_control_data$subjectkey %in% G032_subs_already_sampled &
                          G032_control_data$eventname == "2_year_follow_up_y_arm_1",]
    G032_followup_control_sample_to_add <-
      G032_followup_control_sample_to_add[sample(
        nrow(G032_followup_control_sample_to_add),
        G032_num_followup_controls_to_sample - G032_current_followup_control_n
      ),]
    G032_followup_control_sample <-
      rbind(G032_followup_control_sample,
            G032_followup_control_sample_to_add)
  }
}

#5.1324 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G032_baseline_control_sample$group <- rep("control")
G032_followup_control_sample$group <- rep("control")

#5.1325 Merge the baseline and followup G032 samples 
G032_control_sample <- merge(G032_baseline_control_sample, G032_followup_control_sample, all = TRUE)

#5.1326 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G032_cn_sample <- ifelse(length(unique(G032_control_sample$subjectkey)) == length(G032_control_sample$subjectkey), TRUE, FALSE)

#5.1327 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, G032_control_sample, all = TRUE)


## G054 Resampling ##

#5.141 Filter the control data for just Site G054
G054_control_data <- raw_control_data %>%
  filter(site_name == "G054") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.142 Combine subject data at one site with n and pct values for that site
G054_merged_pct <- merge(site_visit_merged_pct, G054_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.143 Create a df for just baseline and a df for just followup data 
G054_control_baseline_only <- G054_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)

#5.144 save number of baseline only subjects and followup only subjects
G054_baseline_only_sub_n <- as.numeric(nrow(G054_control_baseline_only))

#5.145 save the subject IDs for the baseline only subjects and followup only subjects
G054_baseline_control_subs_only <- G054_control_baseline_only$subjectkey

#5.146 save number of unique control subjects at this site
G054_num_unique_control_subs <- as.numeric(length(unique(G054_merged_pct$subjectkey)))

#5.147 save subsetted dataframe for only G054 GAD subjects
G054_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G054") 

#5.148 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G054_control_baseline_only$group <- rep("control")

#5.149 Merge the baseline and followup G054 samples 
G054_control_sample <- G054_control_baseline_only %>% dplyr::select(-c(GAD_site_visit_n, GAD_site_n, GAD_site_pct, control_site_n, control_site_visit_n, control_site_pct, max_controls))

#5.1410 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G054_cn_sample <- ifelse(length(unique(G054_control_sample$subjectkey)) == length(G054_control_sample$subjectkey), TRUE, FALSE)

#5.1411 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, G054_control_sample, all = TRUE)


## G075 Resampling ##

#5.151 Filter the control data for just Site G075
G075_control_data <- raw_control_data %>%
  filter(site_name == "G075") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.152 Combine subject data at one site with n and pct values for that site
G075_merged_pct <- merge(site_visit_merged_pct, G075_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.153 Create a df for just baseline and a df for just followup data 
G075_control_baseline_only <- G075_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)

#5.154 save number of baseline only subjects and followup only subjects
G075_baseline_only_sub_n <- as.numeric(nrow(G075_control_baseline_only))

#5.155 save the subject IDs for the baseline only subjects and followup only subjects
G075_baseline_control_subs_only <- G075_control_baseline_only$subjectkey

#5.156 save number of unique control subjects at this site
G075_num_unique_control_subs <- as.numeric(length(unique(G075_merged_pct$subjectkey)))

#5.157 save subsetted dataframe for only G075 GAD subjects
G075_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G075") 

#5.158 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G075_control_baseline_only$group <- rep("control")

#5.159 Merge the baseline and followup G075 samples 
G075_control_sample <- G075_control_baseline_only %>% dplyr::select(-c(GAD_site_visit_n, GAD_site_n, GAD_site_pct, control_site_n, control_site_visit_n, control_site_pct, max_controls))

#5.1510 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G075_cn_sample <- ifelse(length(unique(G075_control_sample$subjectkey)) == length(G075_control_sample$subjectkey), TRUE, FALSE)

#5.1511 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, G075_control_sample, all = TRUE)


## G087 Resampling ##

#5.161 Filter the control data for just Site G087
G087_control_data <- raw_control_data %>%
  filter(site_name == "G087") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.162 Combine subject data at one site with n and pct values for that site
G087_merged_pct <- merge(site_visit_merged_pct, G087_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.163 Create a df for just baseline and a df for just followup data 
G087_control_baseline_only <- G087_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
G087_control_followup_only <- G087_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.164 save number of baseline only subjects and followup only subjects
G087_baseline_only_sub_n <- as.numeric(nrow(G087_control_baseline_only))
G087_followup_only_sub_n <- as.numeric(nrow(G087_control_followup_only))

#5.165 save the subject IDs for the baseline only subjects and followup only subjects
G087_baseline_control_subs_only <- G087_control_baseline_only$subjectkey
G087_followup_control_subs_only <- G087_control_followup_only$subjectkey

#5.166 save number of unique control subjects at this site
G087_num_unique_control_subs <- as.numeric(length(unique(G087_merged_pct$subjectkey)))

#5.167 save subsetted dataframe for only G087 GAD subjects
G087_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "G087") 

#5.168 calculate total number of subjects at site G087
G087_total_GAD_sub_n <- as.numeric(length(unique(G087_GAD_subs$subjectkey)))

#5.169 further subset the dataframe to only contain baseline scans or followup GAD scans
G087_baseline_GAD_data <- G087_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
G087_followup_GAD_data <- G087_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1610 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
G087_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G087" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
G087_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "G087" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1611 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
G087_both_timepoint_sub_list <- G087_control_data$subjectkey[duplicated(G087_control_data$subjectkey) ]
G087_both_timepoint_sub_list_n <- as.numeric(length(unique(G087_both_timepoint_sub_list)))

#5.1612 Subset a dataframe of those subjects
G087_both_timepoint_subs <- G087_control_data[G087_control_data$subjectkey %in% G087_both_timepoint_sub_list, ]

#5.1613 Create vectors for the max number of baseline and max number of followup scans 
G087_max_baseline_controls <- as.numeric(G087_baseline_only_sub_n + G087_both_timepoint_sub_list_n)
G087_max_followup_controls <- as.numeric(G087_followup_only_sub_n + G087_both_timepoint_sub_list_n)

#5.1614 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
G087_baseline_controls_desired <- as.numeric(G087_num_unique_control_subs * G087_baseline_GAD_site_visit_percent)
G087_followup_controls_desired <- as.numeric(G087_num_unique_control_subs * G087_followup_GAD_site_visit_percent)

#5.1615 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
G087_sufficient_baseline_control_n <- if_else(G087_max_baseline_controls > G087_baseline_controls_desired, TRUE, FALSE)
G087_sufficient_followup_control_n <- if_else(G087_max_followup_controls > G087_followup_controls_desired, TRUE, FALSE)

#5.1616 Create an adjusted total control N to sample in the event either of the above checks are failed 
G087_adjusted_total_control_n_if_insufficient_baseline <- floor(G087_max_baseline_controls / G087_baseline_GAD_site_visit_percent)
G087_adjusted_total_control_n_if_insufficient_followup <- floor(G087_max_followup_controls / G087_followup_GAD_site_visit_percent)

#5.1617 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
G087_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1618 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!G087_sufficient_baseline_control_n &
    G087_sufficient_followup_control_n) {
  G087_updated_control_counts$num_baseline_controls <-
    G087_max_baseline_controls
  G087_updated_control_counts$num_followup_controls <-
    floor(
      G087_adjusted_total_control_n_if_insufficient_baseline * G087_followup_GAD_site_visit_percent
    )
} else if (G087_sufficient_baseline_control_n &
           !G087_sufficient_followup_control_n) {
  G087_updated_control_counts$num_baseline_controls <-
    floor(
      G087_adjusted_total_control_n_if_insufficient_followup * G087_baseline_GAD_site_visit_percent
    )
  G087_updated_control_counts$num_followup_controls <-
    G087_max_followup_controls
} else {
  G087_updated_control_counts$num_baseline_controls <-
    G087_baseline_controls_desired
  G087_updated_control_counts$num_followup_controls <-
    G087_followup_controls_desired
}

#5.1619 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
G087_baseline_control_sample <- G087_control_data[G087_control_data$subjectkey %in% G087_baseline_control_subs_only, ]
G087_followup_control_sample <- G087_control_data[G087_control_data$subjectkey %in% G087_followup_control_subs_only, ] 

#5.1620 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
G087_num_baseline_controls_to_sample <- as.numeric(G087_updated_control_counts$num_baseline_controls)
G087_current_baseline_control_n <- as.numeric(length(unique(G087_baseline_control_sample$subjectkey)))

#5.1621 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (G087_current_baseline_control_n == G087_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G087_current_baseline_control_n > G087_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    G087_baseline_control_sample <-
      G087_baseline_control_sample[sample(nrow(G087_baseline_control_sample),
                                          G087_num_baseline_controls_to_sample),]
  } else if (G087_current_baseline_control_n < G087_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    G087_baseline_control_sample_to_add <-
      G087_control_data[G087_control_data$subjectkey %in% G087_both_timepoint_sub_list &
                          G087_control_data$eventname == "baseline_year_1_arm_1",]
    G087_baseline_control_sample_to_add <-
      G087_baseline_control_sample_to_add[sample(
        nrow(G087_baseline_control_sample_to_add),
        G087_num_baseline_controls_to_sample - G087_current_baseline_control_n),]
    G087_baseline_control_sample <-
      rbind(G087_baseline_control_sample,
            G087_baseline_control_sample_to_add)
  }
}

#5.1622 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
G087_num_followup_controls_to_sample <- as.numeric(G087_updated_control_counts$num_followup_controls)
G087_current_followup_control_n <- as.numeric(length(unique(G087_followup_control_sample$subjectkey)))
G087_subs_already_sampled <- G087_baseline_control_sample$subjectkey

#5.1623 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (G087_current_followup_control_n == G087_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (G087_current_followup_control_n > G087_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    G087_followup_control_sample <-
      G087_followup_control_sample[sample(nrow(G087_followup_control_sample),
                                          G087_num_followup_controls_to_sample),]
  } else if (G087_current_followup_control_n < G087_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    G087_followup_control_sample_to_add <-
      G087_control_data[G087_control_data$subjectkey %in% G087_both_timepoint_sub_list &
                          !G087_control_data$subjectkey %in% G087_subs_already_sampled &
                          G087_control_data$eventname == "2_year_follow_up_y_arm_1",]
    G087_followup_control_sample_to_add <-
      G087_followup_control_sample_to_add[sample(
        nrow(G087_followup_control_sample_to_add),
        G087_num_followup_controls_to_sample - G087_current_followup_control_n),]
    G087_followup_control_sample <-
      rbind(G087_followup_control_sample,
            G087_followup_control_sample_to_add)
  }
}

#5.1624 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
G087_baseline_control_sample$group <- rep("control")
G087_followup_control_sample$group <- rep("control")

#5.1625 Merge the baseline and followup G087 samples 
G087_control_sample <- merge(G087_baseline_control_sample, G087_followup_control_sample, all = TRUE)

#5.1626 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_G087_cn_sample <- ifelse(length(unique(G087_control_sample$subjectkey)) == length(G087_control_sample$subjectkey), TRUE, FALSE)

#5.1627 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, G087_control_sample, all = TRUE)


## P023 Resampling ##

#5.171 Filter the control data for just Site P023
P023_control_data <- raw_control_data %>%
  filter(site_name == "P023") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.172 Combine subject data at one site with n and pct values for that site
P023_merged_pct <- merge(site_visit_merged_pct, P023_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.173 Create a df for just baseline and a df for just followup data 
P023_control_baseline_only <- P023_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
P023_control_followup_only <- P023_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.174 save number of baseline only subjects and followup only subjects
P023_baseline_only_sub_n <- as.numeric(nrow(P023_control_baseline_only))
P023_followup_only_sub_n <- as.numeric(nrow(P023_control_followup_only))

#5.175 save the subject IDs for the baseline only subjects and followup only subjects
P023_baseline_control_subs_only <- P023_control_baseline_only$subjectkey
P023_followup_control_subs_only <- P023_control_followup_only$subjectkey

#5.176 save number of unique control subjects at this site
P023_num_unique_control_subs <- as.numeric(length(unique(P023_merged_pct$subjectkey)))

#5.177 save subsetted dataframe for only P023 GAD subjects
P023_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "P023") 

#5.178 calculate total number of subjects at site P023
P023_total_GAD_sub_n <- as.numeric(length(unique(P023_GAD_subs$subjectkey)))

#5.179 further subset the dataframe to only contain baseline scans or followup GAD scans
P023_baseline_GAD_data <- P023_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
P023_followup_GAD_data <- P023_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1710 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
P023_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P023" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
P023_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P023" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1711 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
P023_both_timepoint_sub_list <- P023_control_data$subjectkey[duplicated(P023_control_data$subjectkey) ]
P023_both_timepoint_sub_list_n <- as.numeric(length(unique(P023_both_timepoint_sub_list)))

#5.1712 Subset a dataframe of those subjects
P023_both_timepoint_subs <- P023_control_data[P023_control_data$subjectkey %in% P023_both_timepoint_sub_list, ]

#5.1713 Create vectors for the max number of baseline and max number of followup scans 
P023_max_baseline_controls <- as.numeric(P023_baseline_only_sub_n + P023_both_timepoint_sub_list_n)
P023_max_followup_controls <- as.numeric(P023_followup_only_sub_n + P023_both_timepoint_sub_list_n)

#5.1714 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
P023_baseline_controls_desired <- as.numeric(P023_num_unique_control_subs * P023_baseline_GAD_site_visit_percent)
P023_followup_controls_desired <- as.numeric(P023_num_unique_control_subs * P023_followup_GAD_site_visit_percent)

#5.1715 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
P023_sufficient_baseline_control_n <- if_else(P023_max_baseline_controls > P023_baseline_controls_desired, TRUE, FALSE)
P023_sufficient_followup_control_n <- if_else(P023_max_followup_controls > P023_followup_controls_desired, TRUE, FALSE)

#5.1716 Create an adjusted total control N to sample in the event either of the above checks are failed 
P023_adjusted_total_control_n_if_insufficient_baseline <- floor(P023_max_baseline_controls / P023_baseline_GAD_site_visit_percent)
P023_adjusted_total_control_n_if_insufficient_followup <- floor(P023_max_followup_controls / P023_followup_GAD_site_visit_percent)

#5.1717 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
P023_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1718 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!P023_sufficient_baseline_control_n &
    P023_sufficient_followup_control_n) {
  P023_updated_control_counts$num_baseline_controls <-
    P023_max_baseline_controls
  P023_updated_control_counts$num_followup_controls <-
    floor(
      P023_adjusted_total_control_n_if_insufficient_baseline * P023_followup_GAD_site_visit_percent
    )
} else if (P023_sufficient_baseline_control_n &
           !P023_sufficient_followup_control_n) {
  P023_updated_control_counts$num_baseline_controls <-
    floor(
      P023_adjusted_total_control_n_if_insufficient_followup * P023_baseline_GAD_site_visit_percent
    )
  P023_updated_control_counts$num_followup_controls <-
    P023_max_followup_controls
} else {
  P023_updated_control_counts$num_baseline_controls <-
    P023_baseline_controls_desired
  P023_updated_control_counts$num_followup_controls <-
    P023_followup_controls_desired
}

#5.1719 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
P023_baseline_control_sample <- P023_control_data[P023_control_data$subjectkey %in% P023_baseline_control_subs_only, ]
P023_followup_control_sample <- P023_control_data[P023_control_data$subjectkey %in% P023_followup_control_subs_only, ] 

#5.1720 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
P023_num_baseline_controls_to_sample <- as.numeric(P023_updated_control_counts$num_baseline_controls)
P023_current_baseline_control_n <- as.numeric(length(unique(P023_baseline_control_sample$subjectkey)))

#5.1721 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (P023_current_baseline_control_n == P023_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P023_current_baseline_control_n > P023_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    P023_baseline_control_sample <-
      P023_baseline_control_sample[sample(nrow(P023_baseline_control_sample),
                                          P023_num_baseline_controls_to_sample),]
  } else if (P023_current_baseline_control_n < P023_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    P023_baseline_control_sample_to_add <-
      P023_control_data[P023_control_data$subjectkey %in% P023_both_timepoint_sub_list &
                          P023_control_data$eventname == "baseline_year_1_arm_1",]
    P023_baseline_control_sample_to_add <-
      P023_baseline_control_sample_to_add[sample(
        nrow(P023_baseline_control_sample_to_add),
        P023_num_baseline_controls_to_sample - P023_current_baseline_control_n),]
    P023_baseline_control_sample <-
      rbind(P023_baseline_control_sample,
            P023_baseline_control_sample_to_add)
  }
}

#5.1722 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
P023_num_followup_controls_to_sample <- as.numeric(P023_updated_control_counts$num_followup_controls)
P023_current_followup_control_n <- as.numeric(length(unique(P023_followup_control_sample$subjectkey)))
P023_subs_already_sampled <- P023_baseline_control_sample$subjectkey

#5.1723 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (P023_current_followup_control_n == P023_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P023_current_followup_control_n > P023_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    P023_followup_control_sample <-
      P023_followup_control_sample[sample(nrow(P023_followup_control_sample),
                                          P023_num_followup_controls_to_sample),]
  } else if (P023_current_followup_control_n < P023_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    P023_followup_control_sample_to_add <-
      P023_control_data[P023_control_data$subjectkey %in% P023_both_timepoint_sub_list &
                          !P023_control_data$subjectkey %in% P023_subs_already_sampled &
                          P023_control_data$eventname == "2_year_follow_up_y_arm_1",]
    P023_followup_control_sample_to_add <-
      P023_followup_control_sample_to_add[sample(
        nrow(P023_followup_control_sample_to_add),
        P023_num_followup_controls_to_sample - P023_current_followup_control_n),]
    P023_followup_control_sample <-
      rbind(P023_followup_control_sample,
            P023_followup_control_sample_to_add)
  }
}

#5.1724 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
P023_baseline_control_sample$group <- rep("control")
P023_followup_control_sample$group <- rep("control")

#5.1725 Merge the baseline and followup P023 samples 
P023_control_sample <- merge(P023_baseline_control_sample, P023_followup_control_sample, all = TRUE)

#5.1726 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_P023_cn_sample <- ifelse(length(unique(P023_control_sample$subjectkey)) == length(P023_control_sample$subjectkey), TRUE, FALSE)

#5.1727 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, P023_control_sample, all = TRUE)


## P043 Resampling ##

#5.181 Filter the control data for just Site P043
P043_control_data <- raw_control_data %>%
  filter(site_name == "P043") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.182 Combine subject data at one site with n and pct values for that site
P043_merged_pct <- merge(site_visit_merged_pct, P043_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.183 Create a df for just baseline and a df for just followup data 
P043_control_baseline_only <- P043_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
P043_control_followup_only <- P043_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.184 save number of baseline only subjects and followup only subjects
P043_baseline_only_sub_n <- as.numeric(nrow(P043_control_baseline_only))
P043_followup_only_sub_n <- as.numeric(nrow(P043_control_followup_only))

#5.185 save the subject IDs for the baseline only subjects and followup only subjects
P043_baseline_control_subs_only <- P043_control_baseline_only$subjectkey
P043_followup_control_subs_only <- P043_control_followup_only$subjectkey

#5.186 save number of unique control subjects at this site
P043_num_unique_control_subs <- as.numeric(length(unique(P043_merged_pct$subjectkey)))

#5.187 save subsetted dataframe for only P043 GAD subjects
P043_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "P043") 

#5.188 calculate total number of subjects at site P043
P043_total_GAD_sub_n <- as.numeric(length(unique(P043_GAD_subs$subjectkey)))

#5.189 further subset the dataframe to only contain baseline scans or followup GAD scans
P043_baseline_GAD_data <- P043_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
P043_followup_GAD_data <- P043_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1810 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
P043_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P043" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
P043_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P043" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1811 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
P043_both_timepoint_sub_list <- P043_control_data$subjectkey[duplicated(P043_control_data$subjectkey) ]
P043_both_timepoint_sub_list_n <- as.numeric(length(unique(P043_both_timepoint_sub_list)))

#5.1812 Subset a dataframe of those subjects
P043_both_timepoint_subs <- P043_control_data[P043_control_data$subjectkey %in% P043_both_timepoint_sub_list, ]

#5.1813 Create vectors for the max number of baseline and max number of followup scans 
P043_max_baseline_controls <- as.numeric(P043_baseline_only_sub_n + P043_both_timepoint_sub_list_n)
P043_max_followup_controls <- as.numeric(P043_followup_only_sub_n + P043_both_timepoint_sub_list_n)

#5.1814 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
P043_baseline_controls_desired <- as.numeric(P043_num_unique_control_subs * P043_baseline_GAD_site_visit_percent)
P043_followup_controls_desired <- as.numeric(P043_num_unique_control_subs * P043_followup_GAD_site_visit_percent)

#5.1815 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
P043_sufficient_baseline_control_n <- if_else(P043_max_baseline_controls > P043_baseline_controls_desired, TRUE, FALSE)
P043_sufficient_followup_control_n <- if_else(P043_max_followup_controls > P043_followup_controls_desired, TRUE, FALSE)

#5.1816 Create an adjusted total control N to sample in the event either of the above checks are failed 
P043_adjusted_total_control_n_if_insufficient_baseline <- floor(P043_max_baseline_controls / P043_baseline_GAD_site_visit_percent)
P043_adjusted_total_control_n_if_insufficient_followup <- floor(P043_max_followup_controls / P043_followup_GAD_site_visit_percent)

#5.1817 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
P043_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1818 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!P043_sufficient_baseline_control_n &
    P043_sufficient_followup_control_n) {
  P043_updated_control_counts$num_baseline_controls <-
    P043_max_baseline_controls
  P043_updated_control_counts$num_followup_controls <-
    floor(
      P043_adjusted_total_control_n_if_insufficient_baseline * P043_followup_GAD_site_visit_percent
    )
} else if (P043_sufficient_baseline_control_n &
           !P043_sufficient_followup_control_n) {
  P043_updated_control_counts$num_baseline_controls <-
    floor(
      P043_adjusted_total_control_n_if_insufficient_followup * P043_baseline_GAD_site_visit_percent
    )
  P043_updated_control_counts$num_followup_controls <-
    P043_max_followup_controls
} else {
  P043_updated_control_counts$num_baseline_controls <-
    P043_baseline_controls_desired
  P043_updated_control_counts$num_followup_controls <-
    P043_followup_controls_desired
}

#5.1819 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
P043_baseline_control_sample <- P043_control_data[P043_control_data$subjectkey %in% P043_baseline_control_subs_only, ]
P043_followup_control_sample <- P043_control_data[P043_control_data$subjectkey %in% P043_followup_control_subs_only, ] 

#5.1820 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
P043_num_baseline_controls_to_sample <- as.numeric(P043_updated_control_counts$num_baseline_controls)
P043_current_baseline_control_n <- as.numeric(length(unique(P043_baseline_control_sample$subjectkey)))

#5.1821 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (P043_current_baseline_control_n == P043_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P043_current_baseline_control_n > P043_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    P043_baseline_control_sample <-
      P043_baseline_control_sample[sample(nrow(P043_baseline_control_sample),
                                          P043_num_baseline_controls_to_sample),]
  } else if (P043_current_baseline_control_n < P043_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    P043_baseline_control_sample_to_add <-
      P043_control_data[P043_control_data$subjectkey %in% P043_both_timepoint_sub_list &
                          P043_control_data$eventname == "baseline_year_1_arm_1",]
    P043_baseline_control_sample_to_add <-
      P043_baseline_control_sample_to_add[sample(
        nrow(P043_baseline_control_sample_to_add),
        P043_num_baseline_controls_to_sample - P043_current_baseline_control_n
      ),]
    P043_baseline_control_sample <-
      rbind(P043_baseline_control_sample,
            P043_baseline_control_sample_to_add)
  }
}

#5.1822 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
P043_num_followup_controls_to_sample <- as.numeric(P043_updated_control_counts$num_followup_controls)
P043_current_followup_control_n <- as.numeric(length(unique(P043_followup_control_sample$subjectkey)))
P043_subs_already_sampled <- P043_baseline_control_sample$subjectkey

#5.1823 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (P043_current_followup_control_n == P043_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P043_current_followup_control_n > P043_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    P043_followup_control_sample <-
      P043_followup_control_sample[sample(nrow(P043_followup_control_sample),
                                          P043_num_followup_controls_to_sample),]
  } else if (P043_current_followup_control_n < P043_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    P043_followup_control_sample_to_add <-
      P043_control_data[P043_control_data$subjectkey %in% P043_both_timepoint_sub_list &
                          !P043_control_data$subjectkey %in% P043_subs_already_sampled &
                          P043_control_data$eventname == "2_year_follow_up_y_arm_1",]
    P043_followup_control_sample_to_add <-
      P043_followup_control_sample_to_add[sample(
        nrow(P043_followup_control_sample_to_add),
        P043_num_followup_controls_to_sample - P043_current_followup_control_n),]
    P043_followup_control_sample <-
      rbind(P043_followup_control_sample,
            P043_followup_control_sample_to_add)
  }
}

#5.1824 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
P043_baseline_control_sample$group <- rep("control")
P043_followup_control_sample$group <- rep("control")

#5.1825 Merge the baseline and followup P043 samples 
P043_control_sample <- merge(P043_baseline_control_sample, P043_followup_control_sample, all = TRUE)

#5.1826 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_P043_cn_sample <- ifelse(length(unique(P043_control_sample$subjectkey)) == length(P043_control_sample$subjectkey), TRUE, FALSE)

#5.1827 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, P043_control_sample, all = TRUE)


## P064 Resampling ##

#5.191 Filter the control data for just Site P064
P064_control_data <- raw_control_data %>%
  filter(site_name == "P064") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.192 Combine subject data at one site with n and pct values for that site
P064_merged_pct <- merge(site_visit_merged_pct, P064_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.193 Create a df for just baseline and a df for just followup data 
P064_control_baseline_only <- P064_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
P064_control_followup_only <- P064_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.194 save number of baseline only subjects and followup only subjects
P064_baseline_only_sub_n <- as.numeric(nrow(P064_control_baseline_only))
P064_followup_only_sub_n <- as.numeric(nrow(P064_control_followup_only))

#5.195 save the subject IDs for the baseline only subjects and followup only subjects
P064_baseline_control_subs_only <- P064_control_baseline_only$subjectkey
P064_followup_control_subs_only <- P064_control_followup_only$subjectkey

#5.196 save number of unique control subjects at this site
P064_num_unique_control_subs <- as.numeric(length(unique(P064_merged_pct$subjectkey)))

#5.197 save subsetted dataframe for only P064 GAD subjects
P064_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "P064") 

#5.198 calculate total number of subjects at site P064
P064_total_GAD_sub_n <- as.numeric(length(unique(P064_GAD_subs$subjectkey)))

#5.199 further subset the dataframe to only contain baseline scans or followup GAD scans
P064_baseline_GAD_data <- P064_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
P064_followup_GAD_data <- P064_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.1910 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
P064_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P064" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
P064_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "P064" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.1911 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
P064_both_timepoint_sub_list <- P064_control_data$subjectkey[duplicated(P064_control_data$subjectkey) ]
P064_both_timepoint_sub_list_n <- as.numeric(length(unique(P064_both_timepoint_sub_list)))

#5.1912 Subset a dataframe of those subjects
P064_both_timepoint_subs <- P064_control_data[P064_control_data$subjectkey %in% P064_both_timepoint_sub_list, ]

#5.1913 Create vectors for the max number of baseline and max number of followup scans 
P064_max_baseline_controls <- as.numeric(P064_baseline_only_sub_n + P064_both_timepoint_sub_list_n)
P064_max_followup_controls <- as.numeric(P064_followup_only_sub_n + P064_both_timepoint_sub_list_n)

#5.1914 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
P064_baseline_controls_desired <- as.numeric(P064_num_unique_control_subs * P064_baseline_GAD_site_visit_percent)
P064_followup_controls_desired <- as.numeric(P064_num_unique_control_subs * P064_followup_GAD_site_visit_percent)

#5.1915 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
P064_sufficient_baseline_control_n <- if_else(P064_max_baseline_controls > P064_baseline_controls_desired, TRUE, FALSE)
P064_sufficient_followup_control_n <- if_else(P064_max_followup_controls > P064_followup_controls_desired, TRUE, FALSE)

#5.1916 Create an adjusted total control N to sample in the event either of the above checks are failed 
P064_adjusted_total_control_n_if_insufficient_baseline <- floor(P064_max_baseline_controls / P064_baseline_GAD_site_visit_percent)
P064_adjusted_total_control_n_if_insufficient_followup <- floor(P064_max_followup_controls / P064_followup_GAD_site_visit_percent)

#5.1917 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
P064_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.1918 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!P064_sufficient_baseline_control_n &
    P064_sufficient_followup_control_n) {
  P064_updated_control_counts$num_baseline_controls <-
    P064_max_baseline_controls
  P064_updated_control_counts$num_followup_controls <-
    floor(
      P064_adjusted_total_control_n_if_insufficient_baseline * P064_followup_GAD_site_visit_percent
    )
} else if (P064_sufficient_baseline_control_n &
           !P064_sufficient_followup_control_n) {
  P064_updated_control_counts$num_baseline_controls <-
    floor(
      P064_adjusted_total_control_n_if_insufficient_followup * P064_baseline_GAD_site_visit_percent
    )
  P064_updated_control_counts$num_followup_controls <-
    P064_max_followup_controls
} else {
  P064_updated_control_counts$num_baseline_controls <-
    P064_baseline_controls_desired
  P064_updated_control_counts$num_followup_controls <-
    P064_followup_controls_desired
}

#5.1919 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
P064_baseline_control_sample <- P064_control_data[P064_control_data$subjectkey %in% P064_baseline_control_subs_only, ]
P064_followup_control_sample <- P064_control_data[P064_control_data$subjectkey %in% P064_followup_control_subs_only, ] 

#5.1920 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
P064_num_baseline_controls_to_sample <- as.numeric(P064_updated_control_counts$num_baseline_controls)
P064_current_baseline_control_n <- as.numeric(length(unique(P064_baseline_control_sample$subjectkey)))

#5.1921 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (P064_current_baseline_control_n == P064_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P064_current_baseline_control_n > P064_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    P064_baseline_control_sample <-
      P064_baseline_control_sample[sample(nrow(P064_baseline_control_sample),
                                          P064_num_baseline_controls_to_sample),]
  } else if (P064_current_baseline_control_n < P064_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    P064_baseline_control_sample_to_add <-
      P064_control_data[P064_control_data$subjectkey %in% P064_both_timepoint_sub_list &
                          P064_control_data$eventname == "baseline_year_1_arm_1",]
    P064_baseline_control_sample_to_add <-
      P064_baseline_control_sample_to_add[sample(
        nrow(P064_baseline_control_sample_to_add),
        P064_num_baseline_controls_to_sample - P064_current_baseline_control_n),]
    P064_baseline_control_sample <-
      rbind(P064_baseline_control_sample,
            P064_baseline_control_sample_to_add)
  }
}

#5.1922 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
P064_num_followup_controls_to_sample <- as.numeric(P064_updated_control_counts$num_followup_controls)
P064_current_followup_control_n <- as.numeric(length(unique(P064_followup_control_sample$subjectkey)))
P064_subs_already_sampled <- P064_baseline_control_sample$subjectkey

#5.1923 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (P064_current_followup_control_n == P064_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (P064_current_followup_control_n > P064_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    P064_followup_control_sample <-
      P064_followup_control_sample[sample(nrow(P064_followup_control_sample),
                                          P064_num_followup_controls_to_sample),]
  } else if (P064_current_followup_control_n < P064_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    P064_followup_control_sample_to_add <-
      P064_control_data[P064_control_data$subjectkey %in% P064_both_timepoint_sub_list &
                          !P064_control_data$subjectkey %in% P064_subs_already_sampled &
                          P064_control_data$eventname == "2_year_follow_up_y_arm_1",]
    P064_followup_control_sample_to_add <-
      P064_followup_control_sample_to_add[sample(
        nrow(P064_followup_control_sample_to_add),
        P064_num_followup_controls_to_sample - P064_current_followup_control_n),]
    P064_followup_control_sample <-
      rbind(P064_followup_control_sample,
            P064_followup_control_sample_to_add)
  }
}

#5.1924 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
P064_baseline_control_sample$group <- rep("control")
P064_followup_control_sample$group <- rep("control")

#5.1925 Merge the baseline and followup P064 samples 
P064_control_sample <- merge(P064_baseline_control_sample, P064_followup_control_sample, all = TRUE)

#5.1926 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_P064_cn_sample <- ifelse(length(unique(P064_control_sample$subjectkey)) == length(P064_control_sample$subjectkey), TRUE, FALSE)

#5.1927 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, P064_control_sample, all = TRUE)


## S011 Resampling ##

#5.201 Filter the control data for just Site S011
S011_control_data <- raw_control_data %>%
  filter(site_name == "S011") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.202 Combine subject data at one site with n and pct values for that site
S011_merged_pct <- merge(site_visit_merged_pct, S011_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.203 Create a df for just baseline and a df for just followup data 
S011_control_baseline_only <- S011_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S011_control_followup_only <- S011_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.204 save number of baseline only subjects and followup only subjects
S011_baseline_only_sub_n <- as.numeric(nrow(S011_control_baseline_only))
S011_followup_only_sub_n <- as.numeric(nrow(S011_control_followup_only))

#5.205 save the subject IDs for the baseline only subjects and followup only subjects
S011_baseline_control_subs_only <- S011_control_baseline_only$subjectkey
S011_followup_control_subs_only <- S011_control_followup_only$subjectkey

#5.206 save number of unique control subjects at this site
S011_num_unique_control_subs <- as.numeric(length(unique(S011_merged_pct$subjectkey)))

#5.207 save subsetted dataframe for only S011 GAD subjects
S011_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S011") 

#5.208 calculate total number of subjects at site S011
S011_total_GAD_sub_n <- as.numeric(length(unique(S011_GAD_subs$subjectkey)))

#5.209 further subset the dataframe to only contain baseline scans or followup GAD scans
S011_baseline_GAD_data <- S011_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S011_followup_GAD_data <- S011_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2010 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S011_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S011" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S011_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S011" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2011 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S011_both_timepoint_sub_list <- S011_control_data$subjectkey[duplicated(S011_control_data$subjectkey) ]
S011_both_timepoint_sub_list_n <- as.numeric(length(unique(S011_both_timepoint_sub_list)))

#5.2012 Subset a dataframe of those subjects
S011_both_timepoint_subs <- S011_control_data[S011_control_data$subjectkey %in% S011_both_timepoint_sub_list, ]

#5.2013 Create vectors for the max number of baseline and max number of followup scans 
S011_max_baseline_controls <- as.numeric(S011_baseline_only_sub_n + S011_both_timepoint_sub_list_n)
S011_max_followup_controls <- as.numeric(S011_followup_only_sub_n + S011_both_timepoint_sub_list_n)

#5.2014 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S011_baseline_controls_desired <- as.numeric(S011_num_unique_control_subs * S011_baseline_GAD_site_visit_percent)
S011_followup_controls_desired <- as.numeric(S011_num_unique_control_subs * S011_followup_GAD_site_visit_percent)

#5.2015 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S011_sufficient_baseline_control_n <- if_else(S011_max_baseline_controls > S011_baseline_controls_desired, TRUE, FALSE)
S011_sufficient_followup_control_n <- if_else(S011_max_followup_controls > S011_followup_controls_desired, TRUE, FALSE)

#5.2016 Create an adjusted total control N to sample in the event either of the above checks are failed 
S011_adjusted_total_control_n_if_insufficient_baseline <- floor(S011_max_baseline_controls / S011_baseline_GAD_site_visit_percent)
S011_adjusted_total_control_n_if_insufficient_followup <- floor(S011_max_followup_controls / S011_followup_GAD_site_visit_percent)

#5.2017 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S011_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2018 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S011_sufficient_baseline_control_n &
    S011_sufficient_followup_control_n) {
  S011_updated_control_counts$num_baseline_controls <-
    S011_max_baseline_controls
  S011_updated_control_counts$num_followup_controls <-
    floor(
      S011_adjusted_total_control_n_if_insufficient_baseline * S011_followup_GAD_site_visit_percent
    )
} else if (S011_sufficient_baseline_control_n &
           !S011_sufficient_followup_control_n) {
  S011_updated_control_counts$num_baseline_controls <-
    floor(
      S011_adjusted_total_control_n_if_insufficient_followup * S011_baseline_GAD_site_visit_percent
    )
  S011_updated_control_counts$num_followup_controls <-
    S011_max_followup_controls
} else {
  S011_updated_control_counts$num_baseline_controls <-
    S011_baseline_controls_desired
  S011_updated_control_counts$num_followup_controls <-
    S011_followup_controls_desired
}

#5.2019 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S011_baseline_control_sample <- S011_control_data[S011_control_data$subjectkey %in% S011_baseline_control_subs_only, ]
S011_followup_control_sample <- S011_control_data[S011_control_data$subjectkey %in% S011_followup_control_subs_only, ] 

#5.2020 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S011_num_baseline_controls_to_sample <- as.numeric(S011_updated_control_counts$num_baseline_controls)
S011_current_baseline_control_n <- as.numeric(length(unique(S011_baseline_control_sample$subjectkey)))

#5.2021 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S011_current_baseline_control_n == S011_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S011_current_baseline_control_n > S011_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S011_baseline_control_sample <-
      S011_baseline_control_sample[sample(nrow(S011_baseline_control_sample),
                                          S011_num_baseline_controls_to_sample),]
  } else if (S011_current_baseline_control_n < S011_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S011_baseline_control_sample_to_add <-
      S011_control_data[S011_control_data$subjectkey %in% S011_both_timepoint_sub_list &
                          S011_control_data$eventname == "baseline_year_1_arm_1",]
    S011_baseline_control_sample_to_add <-
      S011_baseline_control_sample_to_add[sample(
        nrow(S011_baseline_control_sample_to_add),
        S011_num_baseline_controls_to_sample - S011_current_baseline_control_n),]
    S011_baseline_control_sample <-
      rbind(S011_baseline_control_sample,
            S011_baseline_control_sample_to_add)
  }
}

#5.2022 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S011_num_followup_controls_to_sample <- as.numeric(S011_updated_control_counts$num_followup_controls)
S011_current_followup_control_n <- as.numeric(length(unique(S011_followup_control_sample$subjectkey)))
S011_subs_already_sampled <- S011_baseline_control_sample$subjectkey

#5.2023 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S011_current_followup_control_n == S011_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S011_current_followup_control_n > S011_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S011_followup_control_sample <-
      S011_followup_control_sample[sample(nrow(S011_followup_control_sample),
                                          S011_num_followup_controls_to_sample),]
  } else if (S011_current_followup_control_n < S011_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S011_followup_control_sample_to_add <-
      S011_control_data[S011_control_data$subjectkey %in% S011_both_timepoint_sub_list &
                          !S011_control_data$subjectkey %in% S011_subs_already_sampled &
                          S011_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S011_followup_control_sample_to_add <-
      S011_followup_control_sample_to_add[sample(
        nrow(S011_followup_control_sample_to_add),
        S011_num_followup_controls_to_sample - S011_current_followup_control_n),]
    S011_followup_control_sample <-
      rbind(S011_followup_control_sample,
            S011_followup_control_sample_to_add)
  }
}

#5.2024 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S011_baseline_control_sample$group <- rep("control")
S011_followup_control_sample$group <- rep("control")

#5.2025 Merge the baseline and followup S011 samples 
S011_control_sample <- merge(S011_baseline_control_sample, S011_followup_control_sample, all = TRUE)

#5.2026 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S011_cn_sample <- ifelse(length(unique(S011_control_sample$subjectkey)) == length(S011_control_sample$subjectkey), TRUE, FALSE)

#5.2027 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S011_control_sample, all = TRUE)


## S012 Resampling ##

#5.211 Filter the control data for just Site S012
S012_control_data <- raw_control_data %>%
  filter(site_name == "S012") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.212 Combine subject data at one site with n and pct values for that site
S012_merged_pct <- merge(site_visit_merged_pct, S012_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.213 Create a df for just baseline and a df for just followup data 
S012_control_baseline_only <- S012_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S012_control_followup_only <- S012_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.214 save number of baseline only subjects and followup only subjects
S012_baseline_only_sub_n <- as.numeric(nrow(S012_control_baseline_only))
S012_followup_only_sub_n <- as.numeric(nrow(S012_control_followup_only))

#5.215 save the subject IDs for the baseline only subjects and followup only subjects
S012_baseline_control_subs_only <- S012_control_baseline_only$subjectkey
S012_followup_control_subs_only <- S012_control_followup_only$subjectkey

#5.216 save number of unique control subjects at this site
S012_num_unique_control_subs <- as.numeric(length(unique(S012_merged_pct$subjectkey)))

#5.217 save subsetted dataframe for only S012 GAD subjects
S012_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S012") 

#5.218 calculate total number of subjects at site S012
S012_total_GAD_sub_n <- as.numeric(length(unique(S012_GAD_subs$subjectkey)))

#5.219 further subset the dataframe to only contain baseline scans or followup GAD scans
S012_baseline_GAD_data <- S012_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S012_followup_GAD_data <- S012_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2110 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S012_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S012" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S012_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S012" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2111 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S012_both_timepoint_sub_list <- S012_control_data$subjectkey[duplicated(S012_control_data$subjectkey) ]
S012_both_timepoint_sub_list_n <- as.numeric(length(unique(S012_both_timepoint_sub_list)))

#5.2112 Subset a dataframe of those subjects
S012_both_timepoint_subs <- S012_control_data[S012_control_data$subjectkey %in% S012_both_timepoint_sub_list, ]

#5.2113 Create vectors for the max number of baseline and max number of followup scans 
S012_max_baseline_controls <- as.numeric(S012_baseline_only_sub_n + S012_both_timepoint_sub_list_n)
S012_max_followup_controls <- as.numeric(S012_followup_only_sub_n + S012_both_timepoint_sub_list_n)

#5.2114 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S012_baseline_controls_desired <- as.numeric(S012_num_unique_control_subs * S012_baseline_GAD_site_visit_percent)
S012_followup_controls_desired <- as.numeric(S012_num_unique_control_subs * S012_followup_GAD_site_visit_percent)

#5.2115 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S012_sufficient_baseline_control_n <- if_else(S012_max_baseline_controls > S012_baseline_controls_desired, TRUE, FALSE)
S012_sufficient_followup_control_n <- if_else(S012_max_followup_controls > S012_followup_controls_desired, TRUE, FALSE)

#5.2116 Create an adjusted total control N to sample in the event either of the above checks are failed 
S012_adjusted_total_control_n_if_insufficient_baseline <- floor(S012_max_baseline_controls / S012_baseline_GAD_site_visit_percent)
S012_adjusted_total_control_n_if_insufficient_followup <- floor(S012_max_followup_controls / S012_followup_GAD_site_visit_percent)

#5.2117 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S012_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2118 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S012_sufficient_baseline_control_n &
    S012_sufficient_followup_control_n) {
  S012_updated_control_counts$num_baseline_controls <-
    S012_max_baseline_controls
  S012_updated_control_counts$num_followup_controls <-
    floor(
      S012_adjusted_total_control_n_if_insufficient_baseline * S012_followup_GAD_site_visit_percent
    )
} else if (S012_sufficient_baseline_control_n &
           !S012_sufficient_followup_control_n) {
  S012_updated_control_counts$num_baseline_controls <-
    floor(
      S012_adjusted_total_control_n_if_insufficient_followup * S012_baseline_GAD_site_visit_percent
    )
  S012_updated_control_counts$num_followup_controls <-
    S012_max_followup_controls
} else {
  S012_updated_control_counts$num_baseline_controls <-
    S012_baseline_controls_desired
  S012_updated_control_counts$num_followup_controls <-
    S012_followup_controls_desired
}

#5.2119 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S012_baseline_control_sample <- S012_control_data[S012_control_data$subjectkey %in% S012_baseline_control_subs_only, ]
S012_followup_control_sample <- S012_control_data[S012_control_data$subjectkey %in% S012_followup_control_subs_only, ] 

#5.2120 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S012_num_baseline_controls_to_sample <- as.numeric(S012_updated_control_counts$num_baseline_controls)
S012_current_baseline_control_n <- as.numeric(length(unique(S012_baseline_control_sample$subjectkey)))

#5.2121 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S012_current_baseline_control_n == S012_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S012_current_baseline_control_n > S012_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S012_baseline_control_sample <-
      S012_baseline_control_sample[sample(nrow(S012_baseline_control_sample),
                                          S012_num_baseline_controls_to_sample),]
  } else if (S012_current_baseline_control_n < S012_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S012_baseline_control_sample_to_add <-
      S012_control_data[S012_control_data$subjectkey %in% S012_both_timepoint_sub_list &
                          S012_control_data$eventname == "baseline_year_1_arm_1",]
    S012_baseline_control_sample_to_add <-
      S012_baseline_control_sample_to_add[sample(
        nrow(S012_baseline_control_sample_to_add),
        S012_num_baseline_controls_to_sample - S012_current_baseline_control_n),]
    S012_baseline_control_sample <-
      rbind(S012_baseline_control_sample,
            S012_baseline_control_sample_to_add)
  }
}

#5.2122 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S012_num_followup_controls_to_sample <- as.numeric(S012_updated_control_counts$num_followup_controls)
S012_current_followup_control_n <- as.numeric(length(unique(S012_followup_control_sample$subjectkey)))
S012_subs_already_sampled <- S012_baseline_control_sample$subjectkey

#5.2123 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S012_current_followup_control_n == S012_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S012_current_followup_control_n > S012_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S012_followup_control_sample <-
      S012_followup_control_sample[sample(nrow(S012_followup_control_sample),
                                          S012_num_followup_controls_to_sample),]
  } else if (S012_current_followup_control_n < S012_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S012_followup_control_sample_to_add <-
      S012_control_data[S012_control_data$subjectkey %in% S012_both_timepoint_sub_list &
                          !S012_control_data$subjectkey %in% S012_subs_already_sampled &
                          S012_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S012_followup_control_sample_to_add <-
      S012_followup_control_sample_to_add[sample(
        nrow(S012_followup_control_sample_to_add),
        S012_num_followup_controls_to_sample - S012_current_followup_control_n),]
    S012_followup_control_sample <-
      rbind(S012_followup_control_sample,
            S012_followup_control_sample_to_add)
  }
}

#5.2124 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S012_baseline_control_sample$group <- rep("control")
S012_followup_control_sample$group <- rep("control")

#5.2125 Merge the baseline and followup S012 samples 
S012_control_sample <- merge(S012_baseline_control_sample, S012_followup_control_sample, all = TRUE)

#5.2126 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S012_cn_sample <- ifelse(length(unique(S012_control_sample$subjectkey)) == length(S012_control_sample$subjectkey), TRUE, FALSE)

#5.2127 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S012_control_sample, all = TRUE)


## S013 Resampling ##

#5.221 Filter the control data for just Site S013
S013_control_data <- raw_control_data %>%
  filter(site_name == "S013") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.222 Combine subject data at one site with n and pct values for that site
S013_merged_pct <- merge(site_visit_merged_pct, S013_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.223 Create a df for just baseline and a df for just followup data 
S013_control_baseline_only <- S013_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S013_control_followup_only <- S013_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.224 save number of baseline only subjects and followup only subjects
S013_baseline_only_sub_n <- as.numeric(nrow(S013_control_baseline_only))
S013_followup_only_sub_n <- as.numeric(nrow(S013_control_followup_only))

#5.225 save the subject IDs for the baseline only subjects and followup only subjects
S013_baseline_control_subs_only <- S013_control_baseline_only$subjectkey
S013_followup_control_subs_only <- S013_control_followup_only$subjectkey

#5.226 save number of unique control subjects at this site
S013_num_unique_control_subs <- as.numeric(length(unique(S013_merged_pct$subjectkey)))

#5.227 save subsetted dataframe for only S013 GAD subjects
S013_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S013") 

#5.228 calculate total number of subjects at site S013
S013_total_GAD_sub_n <- as.numeric(length(unique(S013_GAD_subs$subjectkey)))

#5.229 further subset the dataframe to only contain baseline scans or followup GAD scans
S013_baseline_GAD_data <- S013_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S013_followup_GAD_data <- S013_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2210 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S013_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S013" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S013_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S013" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2211 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S013_both_timepoint_sub_list <- S013_control_data$subjectkey[duplicated(S013_control_data$subjectkey) ]
S013_both_timepoint_sub_list_n <- as.numeric(length(unique(S013_both_timepoint_sub_list)))

#5.2212 Subset a dataframe of those subjects
S013_both_timepoint_subs <- S013_control_data[S013_control_data$subjectkey %in% S013_both_timepoint_sub_list, ]

#5.2213 Create vectors for the max number of baseline and max number of followup scans 
S013_max_baseline_controls <- as.numeric(S013_baseline_only_sub_n + S013_both_timepoint_sub_list_n)
S013_max_followup_controls <- as.numeric(S013_followup_only_sub_n + S013_both_timepoint_sub_list_n)

#5.2214 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S013_baseline_controls_desired <- as.numeric(S013_num_unique_control_subs * S013_baseline_GAD_site_visit_percent)
S013_followup_controls_desired <- as.numeric(S013_num_unique_control_subs * S013_followup_GAD_site_visit_percent)

#5.2215 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S013_sufficient_baseline_control_n <- if_else(S013_max_baseline_controls > S013_baseline_controls_desired, TRUE, FALSE)
S013_sufficient_followup_control_n <- if_else(S013_max_followup_controls > S013_followup_controls_desired, TRUE, FALSE)

#5.2216 Create an adjusted total control N to sample in the event either of the above checks are failed 
S013_adjusted_total_control_n_if_insufficient_baseline <- floor(S013_max_baseline_controls / S013_baseline_GAD_site_visit_percent)
S013_adjusted_total_control_n_if_insufficient_followup <- floor(S013_max_followup_controls / S013_followup_GAD_site_visit_percent)

#5.2217 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S013_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2218 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S013_sufficient_baseline_control_n &
    S013_sufficient_followup_control_n) {
  S013_updated_control_counts$num_baseline_controls <-
    S013_max_baseline_controls
  S013_updated_control_counts$num_followup_controls <-
    floor(
      S013_adjusted_total_control_n_if_insufficient_baseline * S013_followup_GAD_site_visit_percent
    )
} else if (S013_sufficient_baseline_control_n &
           !S013_sufficient_followup_control_n) {
  S013_updated_control_counts$num_baseline_controls <-
    floor(
      S013_adjusted_total_control_n_if_insufficient_followup * S013_baseline_GAD_site_visit_percent
    )
  S013_updated_control_counts$num_followup_controls <-
    S013_max_followup_controls
} else {
  S013_updated_control_counts$num_baseline_controls <-
    S013_baseline_controls_desired
  S013_updated_control_counts$num_followup_controls <-
    S013_followup_controls_desired
}

#5.2219 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S013_baseline_control_sample <- S013_control_data[S013_control_data$subjectkey %in% S013_baseline_control_subs_only, ]
S013_followup_control_sample <- S013_control_data[S013_control_data$subjectkey %in% S013_followup_control_subs_only, ] 

#5.2220 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S013_num_baseline_controls_to_sample <- as.numeric(S013_updated_control_counts$num_baseline_controls)
S013_current_baseline_control_n <- as.numeric(length(unique(S013_baseline_control_sample$subjectkey)))

#5.2221 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S013_current_baseline_control_n == S013_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S013_current_baseline_control_n > S013_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S013_baseline_control_sample <-
      S013_baseline_control_sample[sample(nrow(S013_baseline_control_sample),
                                          S013_num_baseline_controls_to_sample),]
  } else if (S013_current_baseline_control_n < S013_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S013_baseline_control_sample_to_add <-
      S013_control_data[S013_control_data$subjectkey %in% S013_both_timepoint_sub_list &
                          S013_control_data$eventname == "baseline_year_1_arm_1",]
    S013_baseline_control_sample_to_add <-
      S013_baseline_control_sample_to_add[sample(
        nrow(S013_baseline_control_sample_to_add),
        S013_num_baseline_controls_to_sample - S013_current_baseline_control_n),]
    S013_baseline_control_sample <-
      rbind(S013_baseline_control_sample,
            S013_baseline_control_sample_to_add)
  }
}

#5.2222 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S013_num_followup_controls_to_sample <- as.numeric(S013_updated_control_counts$num_followup_controls)
S013_current_followup_control_n <- as.numeric(length(unique(S013_followup_control_sample$subjectkey)))
S013_subs_already_sampled <- S013_baseline_control_sample$subjectkey

#5.2223 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S013_current_followup_control_n == S013_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S013_current_followup_control_n > S013_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S013_followup_control_sample <-
      S013_followup_control_sample[sample(nrow(S013_followup_control_sample),
                                          S013_num_followup_controls_to_sample),]
  } else if (S013_current_followup_control_n < S013_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S013_followup_control_sample_to_add <-
      S013_control_data[S013_control_data$subjectkey %in% S013_both_timepoint_sub_list &
                          !S013_control_data$subjectkey %in% S013_subs_already_sampled &
                          S013_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S013_followup_control_sample_to_add <-
      S013_followup_control_sample_to_add[sample(
        nrow(S013_followup_control_sample_to_add),
        S013_num_followup_controls_to_sample - S013_current_followup_control_n),]
    S013_followup_control_sample <-
      rbind(S013_followup_control_sample,
            S013_followup_control_sample_to_add)
  }
}

#5.2224 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S013_baseline_control_sample$group <- rep("control")
S013_followup_control_sample$group <- rep("control")

#5.2225 Merge the baseline and followup S013 samples 
S013_control_sample <- merge(S013_baseline_control_sample, S013_followup_control_sample, all = TRUE)

#5.2226 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S013_cn_sample <- ifelse(length(unique(S013_control_sample$subjectkey)) == length(S013_control_sample$subjectkey), TRUE, FALSE)

#5.2227 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S013_control_sample, all = TRUE)


## S014 Resampling ##

#5.231 Filter the control data for just Site S014
S014_control_data <- raw_control_data %>%
  filter(site_name == "S014") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.232 Combine subject data at one site with n and pct values for that site
S014_merged_pct <- merge(site_visit_merged_pct, S014_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.233 Create a df for just baseline and a df for just followup data 
S014_control_baseline_only <- S014_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S014_control_followup_only <- S014_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.234 save number of baseline only subjects and followup only subjects
S014_baseline_only_sub_n <- as.numeric(nrow(S014_control_baseline_only))
S014_followup_only_sub_n <- as.numeric(nrow(S014_control_followup_only))

#5.235 save the subject IDs for the baseline only subjects and followup only subjects
S014_baseline_control_subs_only <- S014_control_baseline_only$subjectkey
S014_followup_control_subs_only <- S014_control_followup_only$subjectkey

#5.236 save number of unique control subjects at this site
S014_num_unique_control_subs <- as.numeric(length(unique(S014_merged_pct$subjectkey)))

#5.237 save subsetted dataframe for only S014 GAD subjects
S014_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S014") 

#5.238 calculate total number of subjects at site S014
S014_total_GAD_sub_n <- as.numeric(length(unique(S014_GAD_subs$subjectkey)))

#5.239 further subset the dataframe to only contain baseline scans or followup GAD scans
S014_baseline_GAD_data <- S014_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S014_followup_GAD_data <- S014_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2310 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S014_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S014" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S014_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S014" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2311 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S014_both_timepoint_sub_list <- S014_control_data$subjectkey[duplicated(S014_control_data$subjectkey) ]
S014_both_timepoint_sub_list_n <- as.numeric(length(unique(S014_both_timepoint_sub_list)))

#5.2312 Subset a dataframe of those subjects
S014_both_timepoint_subs <- S014_control_data[S014_control_data$subjectkey %in% S014_both_timepoint_sub_list, ]

#5.2313 Create vectors for the max number of baseline and max number of followup scans 
S014_max_baseline_controls <- as.numeric(S014_baseline_only_sub_n + S014_both_timepoint_sub_list_n)
S014_max_followup_controls <- as.numeric(S014_followup_only_sub_n + S014_both_timepoint_sub_list_n)

#5.2314 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S014_baseline_controls_desired <- as.numeric(S014_num_unique_control_subs * S014_baseline_GAD_site_visit_percent)
S014_followup_controls_desired <- as.numeric(S014_num_unique_control_subs * S014_followup_GAD_site_visit_percent)

#5.2315 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S014_sufficient_baseline_control_n <- if_else(S014_max_baseline_controls > S014_baseline_controls_desired, TRUE, FALSE)
S014_sufficient_followup_control_n <- if_else(S014_max_followup_controls > S014_followup_controls_desired, TRUE, FALSE)

#5.2316 Create an adjusted total control N to sample in the event either of the above checks are failed 
S014_adjusted_total_control_n_if_insufficient_baseline <- floor(S014_max_baseline_controls / S014_baseline_GAD_site_visit_percent)
S014_adjusted_total_control_n_if_insufficient_followup <- floor(S014_max_followup_controls / S014_followup_GAD_site_visit_percent)

#5.2317 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S014_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2318 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S014_sufficient_baseline_control_n &
    S014_sufficient_followup_control_n) {
  S014_updated_control_counts$num_baseline_controls <-
    S014_max_baseline_controls
  S014_updated_control_counts$num_followup_controls <-
    floor(
      S014_adjusted_total_control_n_if_insufficient_baseline * S014_followup_GAD_site_visit_percent
    )
} else if (S014_sufficient_baseline_control_n &
           !S014_sufficient_followup_control_n) {
  S014_updated_control_counts$num_baseline_controls <-
    floor(
      S014_adjusted_total_control_n_if_insufficient_followup * S014_baseline_GAD_site_visit_percent
    )
  S014_updated_control_counts$num_followup_controls <-
    S014_max_followup_controls
} else {
  S014_updated_control_counts$num_baseline_controls <-
    S014_baseline_controls_desired
  S014_updated_control_counts$num_followup_controls <-
    S014_followup_controls_desired
}

#5.2319 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S014_baseline_control_sample <- S014_control_data[S014_control_data$subjectkey %in% S014_baseline_control_subs_only, ]
S014_followup_control_sample <- S014_control_data[S014_control_data$subjectkey %in% S014_followup_control_subs_only, ] 

#5.2320 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S014_num_baseline_controls_to_sample <- as.numeric(S014_updated_control_counts$num_baseline_controls)
S014_current_baseline_control_n <- as.numeric(length(unique(S014_baseline_control_sample$subjectkey)))

#5.2321 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S014_current_baseline_control_n == S014_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S014_current_baseline_control_n > S014_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S014_baseline_control_sample <-
      S014_baseline_control_sample[sample(nrow(S014_baseline_control_sample),
                                          S014_num_baseline_controls_to_sample),]
  } else if (S014_current_baseline_control_n < S014_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S014_baseline_control_sample_to_add <-
      S014_control_data[S014_control_data$subjectkey %in% S014_both_timepoint_sub_list &
                          S014_control_data$eventname == "baseline_year_1_arm_1",]
    S014_baseline_control_sample_to_add <-
      S014_baseline_control_sample_to_add[sample(
        nrow(S014_baseline_control_sample_to_add),
        S014_num_baseline_controls_to_sample - S014_current_baseline_control_n),]
    S014_baseline_control_sample <-
      rbind(S014_baseline_control_sample,
            S014_baseline_control_sample_to_add)
  }
}

#5.2322 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S014_num_followup_controls_to_sample <- as.numeric(S014_updated_control_counts$num_followup_controls)
S014_current_followup_control_n <- as.numeric(length(unique(S014_followup_control_sample$subjectkey)))
S014_subs_already_sampled <- S014_baseline_control_sample$subjectkey

#5.2323 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S014_current_followup_control_n == S014_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S014_current_followup_control_n > S014_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S014_followup_control_sample <-
      S014_followup_control_sample[sample(nrow(S014_followup_control_sample),
                                          S014_num_followup_controls_to_sample),]
  } else if (S014_current_followup_control_n < S014_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S014_followup_control_sample_to_add <-
      S014_control_data[S014_control_data$subjectkey %in% S014_both_timepoint_sub_list &
                          !S014_control_data$subjectkey %in% S014_subs_already_sampled &
                          S014_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S014_followup_control_sample_to_add <-
      S014_followup_control_sample_to_add[sample(
        nrow(S014_followup_control_sample_to_add),
        S014_num_followup_controls_to_sample - S014_current_followup_control_n),]
    S014_followup_control_sample <-
      rbind(S014_followup_control_sample,
            S014_followup_control_sample_to_add)
  }
}

#5.2324 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S014_baseline_control_sample$group <- rep("control")
S014_followup_control_sample$group <- rep("control")

#5.2325 Merge the baseline and followup S014 samples 
S014_control_sample <- merge(S014_baseline_control_sample, S014_followup_control_sample, all = TRUE)

#5.2326 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S014_cn_sample <- ifelse(length(unique(S014_control_sample$subjectkey)) == length(S014_control_sample$subjectkey), TRUE, FALSE)

#5.2327 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S014_control_sample, all = TRUE)


## S020 Resampling ##

#5.241 Filter the control data for just Site S020
S020_control_data <- raw_control_data %>%
  filter(site_name == "S020") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.242 Combine subject data at one site with n and pct values for that site
S020_merged_pct <- merge(site_visit_merged_pct, S020_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.243 Create a df for just baseline and a df for just followup data 
S020_control_baseline_only <- S020_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S020_control_followup_only <- S020_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.244 save number of baseline only subjects and followup only subjects
S020_baseline_only_sub_n <- as.numeric(nrow(S020_control_baseline_only))
S020_followup_only_sub_n <- as.numeric(nrow(S020_control_followup_only))

#5.245 save the subject IDs for the baseline only subjects and followup only subjects
S020_baseline_control_subs_only <- S020_control_baseline_only$subjectkey
S020_followup_control_subs_only <- S020_control_followup_only$subjectkey

#5.246 save number of unique control subjects at this site
S020_num_unique_control_subs <- as.numeric(length(unique(S020_merged_pct$subjectkey)))

#5.247 save subsetted dataframe for only S020 GAD subjects
S020_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S020") 

#5.248 calculate total number of subjects at site S020
S020_total_GAD_sub_n <- as.numeric(length(unique(S020_GAD_subs$subjectkey)))

#5.249 further subset the dataframe to only contain baseline scans or followup GAD scans
S020_baseline_GAD_data <- S020_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S020_followup_GAD_data <- S020_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2410 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S020_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S020" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S020_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S020" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2411 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S020_both_timepoint_sub_list <- S020_control_data$subjectkey[duplicated(S020_control_data$subjectkey) ]
S020_both_timepoint_sub_list_n <- as.numeric(length(unique(S020_both_timepoint_sub_list)))

#5.2412 Subset a dataframe of those subjects
S020_both_timepoint_subs <- S020_control_data[S020_control_data$subjectkey %in% S020_both_timepoint_sub_list, ]

#5.2413 Create vectors for the max number of baseline and max number of followup scans 
S020_max_baseline_controls <- as.numeric(S020_baseline_only_sub_n + S020_both_timepoint_sub_list_n)
S020_max_followup_controls <- as.numeric(S020_followup_only_sub_n + S020_both_timepoint_sub_list_n)

#5.2414 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S020_baseline_controls_desired <- as.numeric(S020_num_unique_control_subs * S020_baseline_GAD_site_visit_percent)
S020_followup_controls_desired <- as.numeric(S020_num_unique_control_subs * S020_followup_GAD_site_visit_percent)

#5.2415 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S020_sufficient_baseline_control_n <- if_else(S020_max_baseline_controls > S020_baseline_controls_desired, TRUE, FALSE)
S020_sufficient_followup_control_n <- if_else(S020_max_followup_controls > S020_followup_controls_desired, TRUE, FALSE)

#5.2416 Create an adjusted total control N to sample in the event either of the above checks are failed 
S020_adjusted_total_control_n_if_insufficient_baseline <- floor(S020_max_baseline_controls / S020_baseline_GAD_site_visit_percent)
S020_adjusted_total_control_n_if_insufficient_followup <- floor(S020_max_followup_controls / S020_followup_GAD_site_visit_percent)

#5.2417 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S020_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2418 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S020_sufficient_baseline_control_n &
    S020_sufficient_followup_control_n) {
  S020_updated_control_counts$num_baseline_controls <-
    S020_max_baseline_controls
  S020_updated_control_counts$num_followup_controls <-
    floor(
      S020_adjusted_total_control_n_if_insufficient_baseline * S020_followup_GAD_site_visit_percent
    )
} else if (S020_sufficient_baseline_control_n &
           !S020_sufficient_followup_control_n) {
  S020_updated_control_counts$num_baseline_controls <-
    floor(
      S020_adjusted_total_control_n_if_insufficient_followup * S020_baseline_GAD_site_visit_percent
    )
  S020_updated_control_counts$num_followup_controls <-
    S020_max_followup_controls
} else {
  S020_updated_control_counts$num_baseline_controls <-
    S020_baseline_controls_desired
  S020_updated_control_counts$num_followup_controls <-
    S020_followup_controls_desired
}

#5.2419 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S020_baseline_control_sample <- S020_control_data[S020_control_data$subjectkey %in% S020_baseline_control_subs_only, ]
S020_followup_control_sample <- S020_control_data[S020_control_data$subjectkey %in% S020_followup_control_subs_only, ] 

#5.2420 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S020_num_baseline_controls_to_sample <- as.numeric(S020_updated_control_counts$num_baseline_controls)
S020_current_baseline_control_n <- as.numeric(length(unique(S020_baseline_control_sample$subjectkey)))

#5.2421 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S020_current_baseline_control_n == S020_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S020_current_baseline_control_n > S020_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S020_baseline_control_sample <-
      S020_baseline_control_sample[sample(nrow(S020_baseline_control_sample),
                                          S020_num_baseline_controls_to_sample),]
  } else if (S020_current_baseline_control_n < S020_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S020_baseline_control_sample_to_add <-
      S020_control_data[S020_control_data$subjectkey %in% S020_both_timepoint_sub_list &
                          S020_control_data$eventname == "baseline_year_1_arm_1",]
    S020_baseline_control_sample_to_add <-
      S020_baseline_control_sample_to_add[sample(
        nrow(S020_baseline_control_sample_to_add),
        S020_num_baseline_controls_to_sample - S020_current_baseline_control_n),]
    S020_baseline_control_sample <-
      rbind(S020_baseline_control_sample,
            S020_baseline_control_sample_to_add)
  }
}

#5.2422 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S020_num_followup_controls_to_sample <- as.numeric(S020_updated_control_counts$num_followup_controls)
S020_current_followup_control_n <- as.numeric(length(unique(S020_followup_control_sample$subjectkey)))
S020_subs_already_sampled <- S020_baseline_control_sample$subjectkey

#5.2423 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S020_current_followup_control_n == S020_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S020_current_followup_control_n > S020_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S020_followup_control_sample <-
      S020_followup_control_sample[sample(nrow(S020_followup_control_sample),
                                          S020_num_followup_controls_to_sample),]
  } else if (S020_current_followup_control_n < S020_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S020_followup_control_sample_to_add <-
      S020_control_data[S020_control_data$subjectkey %in% S020_both_timepoint_sub_list &
                          !S020_control_data$subjectkey %in% S020_subs_already_sampled &
                          S020_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S020_followup_control_sample_to_add <-
      S020_followup_control_sample_to_add[sample(
        nrow(S020_followup_control_sample_to_add),
        S020_num_followup_controls_to_sample - S020_current_followup_control_n),]
    S020_followup_control_sample <-
      rbind(S020_followup_control_sample,
            S020_followup_control_sample_to_add)
  }
}

#5.2424 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S020_baseline_control_sample$group <- rep("control")
S020_followup_control_sample$group <- rep("control")

#5.2425 Merge the baseline and followup S020 samples 
S020_control_sample <- merge(S020_baseline_control_sample, S020_followup_control_sample, all = TRUE)

#5.2426 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S020_cn_sample <- ifelse(length(unique(S020_control_sample$subjectkey)) == length(S020_control_sample$subjectkey), TRUE, FALSE)

#5.2427 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S020_control_sample, all = TRUE)


## S021 Resampling ##

#5.251 Filter the control data for just Site S021
S021_control_data <- raw_control_data %>%
  filter(site_name == "S021") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.252 Combine subject data at one site with n and pct values for that site
S021_merged_pct <- merge(site_visit_merged_pct, S021_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.253 Create a df for just baseline and a df for just followup data 
S021_control_baseline_only <- S021_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S021_control_followup_only <- S021_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.254 save number of baseline only subjects and followup only subjects
S021_baseline_only_sub_n <- as.numeric(nrow(S021_control_baseline_only))
S021_followup_only_sub_n <- as.numeric(nrow(S021_control_followup_only))

#5.255 save the subject IDs for the baseline only subjects and followup only subjects
S021_baseline_control_subs_only <- S021_control_baseline_only$subjectkey
S021_followup_control_subs_only <- S021_control_followup_only$subjectkey

#5.256 save number of unique control subjects at this site
S021_num_unique_control_subs <- as.numeric(length(unique(S021_merged_pct$subjectkey)))

#5.257 save subsetted dataframe for only S021 GAD subjects
S021_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S021") 

#5.258 calculate total number of subjects at site S021
S021_total_GAD_sub_n <- as.numeric(length(unique(S021_GAD_subs$subjectkey)))

#5.259 further subset the dataframe to only contain baseline scans or followup GAD scans
S021_baseline_GAD_data <- S021_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S021_followup_GAD_data <- S021_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2510 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S021_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S021" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S021_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S021" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2511 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S021_both_timepoint_sub_list <- S021_control_data$subjectkey[duplicated(S021_control_data$subjectkey) ]
S021_both_timepoint_sub_list_n <- as.numeric(length(unique(S021_both_timepoint_sub_list)))

#5.2512 Subset a dataframe of those subjects
S021_both_timepoint_subs <- S021_control_data[S021_control_data$subjectkey %in% S021_both_timepoint_sub_list, ]

#5.2513 Create vectors for the max number of baseline and max number of followup scans 
S021_max_baseline_controls <- as.numeric(S021_baseline_only_sub_n + S021_both_timepoint_sub_list_n)
S021_max_followup_controls <- as.numeric(S021_followup_only_sub_n + S021_both_timepoint_sub_list_n)

#5.2514 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S021_baseline_controls_desired <- as.numeric(S021_num_unique_control_subs * S021_baseline_GAD_site_visit_percent)
S021_followup_controls_desired <- as.numeric(S021_num_unique_control_subs * S021_followup_GAD_site_visit_percent)

#5.2515 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S021_sufficient_baseline_control_n <- if_else(S021_max_baseline_controls > S021_baseline_controls_desired, TRUE, FALSE)
S021_sufficient_followup_control_n <- if_else(S021_max_followup_controls > S021_followup_controls_desired, TRUE, FALSE)

#5.2516 Create an adjusted total control N to sample in the event either of the above checks are failed 
S021_adjusted_total_control_n_if_insufficient_baseline <- floor(S021_max_baseline_controls / S021_baseline_GAD_site_visit_percent)
S021_adjusted_total_control_n_if_insufficient_followup <- floor(S021_max_followup_controls / S021_followup_GAD_site_visit_percent)

#5.2517 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S021_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2518 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S021_sufficient_baseline_control_n &
    S021_sufficient_followup_control_n) {
  S021_updated_control_counts$num_baseline_controls <-
    S021_max_baseline_controls
  S021_updated_control_counts$num_followup_controls <-
    floor(
      S021_adjusted_total_control_n_if_insufficient_baseline * S021_followup_GAD_site_visit_percent
    )
} else if (S021_sufficient_baseline_control_n &
           !S021_sufficient_followup_control_n) {
  S021_updated_control_counts$num_baseline_controls <-
    floor(
      S021_adjusted_total_control_n_if_insufficient_followup * S021_baseline_GAD_site_visit_percent
    )
  S021_updated_control_counts$num_followup_controls <-
    S021_max_followup_controls
} else {
  S021_updated_control_counts$num_baseline_controls <-
    S021_baseline_controls_desired
  S021_updated_control_counts$num_followup_controls <-
    S021_followup_controls_desired
}

#5.2519 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S021_baseline_control_sample <- S021_control_data[S021_control_data$subjectkey %in% S021_baseline_control_subs_only, ]
S021_followup_control_sample <- S021_control_data[S021_control_data$subjectkey %in% S021_followup_control_subs_only, ]

#5.2520 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S021_num_baseline_controls_to_sample <- as.numeric(S021_updated_control_counts$num_baseline_controls)
S021_current_baseline_control_n <- as.numeric(length(unique(S021_baseline_control_sample$subjectkey)))

#5.2521 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S021_current_baseline_control_n == S021_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S021_current_baseline_control_n > S021_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S021_baseline_control_sample <-
      S021_baseline_control_sample[sample(nrow(S021_baseline_control_sample),
                                          S021_num_baseline_controls_to_sample),]
  } else if (S021_current_baseline_control_n < S021_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S021_baseline_control_sample_to_add <-
      S021_control_data[S021_control_data$subjectkey %in% S021_both_timepoint_sub_list &
                          S021_control_data$eventname == "baseline_year_1_arm_1",]
    S021_baseline_control_sample_to_add <-
      S021_baseline_control_sample_to_add[sample(
        nrow(S021_baseline_control_sample_to_add),
        S021_num_baseline_controls_to_sample - S021_current_baseline_control_n),]
    S021_baseline_control_sample <-
      rbind(S021_baseline_control_sample,
            S021_baseline_control_sample_to_add)
  }
}

#5.2522 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S021_num_followup_controls_to_sample <- as.numeric(S021_updated_control_counts$num_followup_controls)
S021_current_followup_control_n <- as.numeric(length(unique(S021_followup_control_sample$subjectkey)))
S021_subs_already_sampled <- S021_baseline_control_sample$subjectkey

#5.2523 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S021_current_followup_control_n == S021_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S021_current_followup_control_n > S021_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S021_followup_control_sample <-
      S021_followup_control_sample[sample(nrow(S021_followup_control_sample),
                                          S021_num_followup_controls_to_sample),]
  } else if (S021_current_followup_control_n < S021_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S021_followup_control_sample_to_add <-
      S021_control_data[S021_control_data$subjectkey %in% S021_both_timepoint_sub_list &
                          !S021_control_data$subjectkey %in% S021_subs_already_sampled &
                          S021_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S021_followup_control_sample_to_add <-
      S021_followup_control_sample_to_add[sample(
        nrow(S021_followup_control_sample_to_add),
        S021_num_followup_controls_to_sample - S021_current_followup_control_n),]
    S021_followup_control_sample <-
      rbind(S021_followup_control_sample,
            S021_followup_control_sample_to_add)
  }
}

#5.2524 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S021_baseline_control_sample$group <- rep("control")
S021_followup_control_sample$group <- rep("control")

#5.2525 Merge the baseline and followup S021 samples 
S021_control_sample <- merge(S021_baseline_control_sample, S021_followup_control_sample, all = TRUE)

#5.2526 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S021_cn_sample <- ifelse(length(unique(S021_control_sample$subjectkey)) == length(S021_control_sample$subjectkey), TRUE, FALSE)

#5.2527 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S021_control_sample, all = TRUE)


## S022 Resampling ##

#5.261 Filter the control data for just Site S022
S022_control_data <- raw_control_data %>%
  filter(site_name == "S022") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.262 Combine subject data at one site with n and pct values for that site
S022_merged_pct <- merge(site_visit_merged_pct, S022_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.263 Create a df for just baseline and a df for just followup data 
S022_control_baseline_only <- S022_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S022_control_followup_only <- S022_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.264 save number of baseline only subjects and followup only subjects
S022_baseline_only_sub_n <- as.numeric(nrow(S022_control_baseline_only))
S022_followup_only_sub_n <- as.numeric(nrow(S022_control_followup_only))

#5.265 save the subject IDs for the baseline only subjects and followup only subjects
S022_baseline_control_subs_only <- S022_control_baseline_only$subjectkey
S022_followup_control_subs_only <- S022_control_followup_only$subjectkey

#5.266 save number of unique control subjects at this site
S022_num_unique_control_subs <- as.numeric(length(unique(S022_merged_pct$subjectkey)))

#5.267 save subsetted dataframe for only S022 GAD subjects
S022_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S022") 

#5.268 calculate total number of subjects at site S022
S022_total_GAD_sub_n <- as.numeric(length(unique(S022_GAD_subs$subjectkey)))

#5.269 further subset the dataframe to only contain baseline scans or followup GAD scans
S022_baseline_GAD_data <- S022_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S022_followup_GAD_data <- S022_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2610 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S022_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S022" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S022_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S022" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2611 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S022_both_timepoint_sub_list <- S022_control_data$subjectkey[duplicated(S022_control_data$subjectkey) ]
S022_both_timepoint_sub_list_n <- as.numeric(length(unique(S022_both_timepoint_sub_list)))

#5.2612 Subset a dataframe of those subjects
S022_both_timepoint_subs <- S022_control_data[S022_control_data$subjectkey %in% S022_both_timepoint_sub_list, ]

#5.2613 Create vectors for the max number of baseline and max number of followup scans 
S022_max_baseline_controls <- as.numeric(S022_baseline_only_sub_n + S022_both_timepoint_sub_list_n)
S022_max_followup_controls <- as.numeric(S022_followup_only_sub_n + S022_both_timepoint_sub_list_n)

#5.2614 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S022_baseline_controls_desired <- as.numeric(S022_num_unique_control_subs * S022_baseline_GAD_site_visit_percent)
S022_followup_controls_desired <- as.numeric(S022_num_unique_control_subs * S022_followup_GAD_site_visit_percent)

#5.2615 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S022_sufficient_baseline_control_n <- if_else(S022_max_baseline_controls > S022_baseline_controls_desired, TRUE, FALSE)
S022_sufficient_followup_control_n <- if_else(S022_max_followup_controls > S022_followup_controls_desired, TRUE, FALSE)

#5.2616 Create an adjusted total control N to sample in the event either of the above checks are failed 
S022_adjusted_total_control_n_if_insufficient_baseline <- floor(S022_max_baseline_controls / S022_baseline_GAD_site_visit_percent)
S022_adjusted_total_control_n_if_insufficient_followup <- floor(S022_max_followup_controls / S022_followup_GAD_site_visit_percent)

#5.2617 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S022_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2618 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S022_sufficient_baseline_control_n &
    S022_sufficient_followup_control_n) {
  S022_updated_control_counts$num_baseline_controls <-
    S022_max_baseline_controls
  S022_updated_control_counts$num_followup_controls <-
    floor(
      S022_adjusted_total_control_n_if_insufficient_baseline * S022_followup_GAD_site_visit_percent
    )
} else if (S022_sufficient_baseline_control_n &
           !S022_sufficient_followup_control_n) {
  S022_updated_control_counts$num_baseline_controls <-
    floor(
      S022_adjusted_total_control_n_if_insufficient_followup * S022_baseline_GAD_site_visit_percent
    )
  S022_updated_control_counts$num_followup_controls <-
    S022_max_followup_controls
} else {
  S022_updated_control_counts$num_baseline_controls <-
    S022_baseline_controls_desired
  S022_updated_control_counts$num_followup_controls <-
    S022_followup_controls_desired
}

#5.2619 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S022_baseline_control_sample <- S022_control_data[S022_control_data$subjectkey %in% S022_baseline_control_subs_only, ]
S022_followup_control_sample <- S022_control_data[S022_control_data$subjectkey %in% S022_followup_control_subs_only, ]

#5.2620 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S022_num_baseline_controls_to_sample <- as.numeric(S022_updated_control_counts$num_baseline_controls)
S022_current_baseline_control_n <- as.numeric(length(unique(S022_baseline_control_sample$subjectkey)))

#5.2621 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S022_current_baseline_control_n == S022_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S022_current_baseline_control_n > S022_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S022_baseline_control_sample <-
      S022_baseline_control_sample[sample(nrow(S022_baseline_control_sample),
                                          S022_num_baseline_controls_to_sample),]
  } else if (S022_current_baseline_control_n < S022_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S022_baseline_control_sample_to_add <-
      S022_control_data[S022_control_data$subjectkey %in% S022_both_timepoint_sub_list &
                          S022_control_data$eventname == "baseline_year_1_arm_1",]
    S022_baseline_control_sample_to_add <-
      S022_baseline_control_sample_to_add[sample(
        nrow(S022_baseline_control_sample_to_add),
        S022_num_baseline_controls_to_sample - S022_current_baseline_control_n),]
    S022_baseline_control_sample <-
      rbind(S022_baseline_control_sample,
            S022_baseline_control_sample_to_add)
  }
}

#5.2622 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S022_num_followup_controls_to_sample <- as.numeric(S022_updated_control_counts$num_followup_controls)
S022_current_followup_control_n <- as.numeric(length(unique(S022_followup_control_sample$subjectkey)))
S022_subs_already_sampled <- S022_baseline_control_sample$subjectkey

#5.2623 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S022_current_followup_control_n == S022_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S022_current_followup_control_n > S022_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S022_followup_control_sample <-
      S022_followup_control_sample[sample(nrow(S022_followup_control_sample),
                                          S022_num_followup_controls_to_sample),]
  } else if (S022_current_followup_control_n < S022_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S022_followup_control_sample_to_add <-
      S022_control_data[S022_control_data$subjectkey %in% S022_both_timepoint_sub_list &
                          !S022_control_data$subjectkey %in% S022_subs_already_sampled &
                          S022_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S022_followup_control_sample_to_add <-
      S022_followup_control_sample_to_add[sample(
        nrow(S022_followup_control_sample_to_add),
        S022_num_followup_controls_to_sample - S022_current_followup_control_n),]
    S022_followup_control_sample <-
      rbind(S022_followup_control_sample,
            S022_followup_control_sample_to_add)
  }
}

#5.2624 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S022_baseline_control_sample$group <- rep("control")
S022_followup_control_sample$group <- rep("control")

#5.2625 Merge the baseline and followup S022 samples 
S022_control_sample <- merge(S022_baseline_control_sample, S022_followup_control_sample, all = TRUE)

#5.2626 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S022_cn_sample <- ifelse(length(unique(S022_control_sample$subjectkey)) == length(S022_control_sample$subjectkey), TRUE, FALSE)

#5.2627 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S022_control_sample, all = TRUE)


## S042 Resampling ##

#5.271 Filter the control data for just Site S042
S042_control_data <- raw_control_data %>%
  filter(site_name == "S042") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.272 Combine subject data at one site with n and pct values for that site
S042_merged_pct <- merge(site_visit_merged_pct, S042_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.273 Create a df for just baseline and a df for just followup data 
S042_control_baseline_only <- S042_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S042_control_followup_only <- S042_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.274 save number of baseline only subjects and followup only subjects
S042_baseline_only_sub_n <- as.numeric(nrow(S042_control_baseline_only))
S042_followup_only_sub_n <- as.numeric(nrow(S042_control_followup_only))

#5.275 save the subject IDs for the baseline only subjects and followup only subjects
S042_baseline_control_subs_only <- S042_control_baseline_only$subjectkey
S042_followup_control_subs_only <- S042_control_followup_only$subjectkey

#5.276 save number of unique control subjects at this site
S042_num_unique_control_subs <- as.numeric(length(unique(S042_merged_pct$subjectkey)))

#5.277 save subsetted dataframe for only S042 GAD subjects
S042_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S042") 

#5.278 calculate total number of subjects at site S042
S042_total_GAD_sub_n <- as.numeric(length(unique(S042_GAD_subs$subjectkey)))

#5.279 further subset the dataframe to only contain baseline scans or followup GAD scans
S042_baseline_GAD_data <- S042_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S042_followup_GAD_data <- S042_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.2710 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S042_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S042" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S042_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S042" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.2711 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S042_both_timepoint_sub_list <- S042_control_data$subjectkey[duplicated(S042_control_data$subjectkey) ]
S042_both_timepoint_sub_list_n <- as.numeric(length(unique(S042_both_timepoint_sub_list)))

#5.2712 Subset a dataframe of those subjects
S042_both_timepoint_subs <- S042_control_data[S042_control_data$subjectkey %in% S042_both_timepoint_sub_list, ]

#5.2713 Create vectors for the max number of baseline and max number of followup scans 
S042_max_baseline_controls <- as.numeric(S042_baseline_only_sub_n + S042_both_timepoint_sub_list_n)
S042_max_followup_controls <- as.numeric(S042_followup_only_sub_n + S042_both_timepoint_sub_list_n)

#5.2714 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S042_baseline_controls_desired <- as.numeric(S042_num_unique_control_subs * S042_baseline_GAD_site_visit_percent)
S042_followup_controls_desired <- as.numeric(S042_num_unique_control_subs * S042_followup_GAD_site_visit_percent)

#5.2715 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S042_sufficient_baseline_control_n <- if_else(S042_max_baseline_controls > S042_baseline_controls_desired, TRUE, FALSE)
S042_sufficient_followup_control_n <- if_else(S042_max_followup_controls > S042_followup_controls_desired, TRUE, FALSE)

#5.2716 Create an adjusted total control N to sample in the event either of the above checks are failed 
S042_adjusted_total_control_n_if_insufficient_baseline <- floor(S042_max_baseline_controls / S042_baseline_GAD_site_visit_percent)
S042_adjusted_total_control_n_if_insufficient_followup <- floor(S042_max_followup_controls / S042_followup_GAD_site_visit_percent)

#5.2717 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S042_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.2718 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S042_sufficient_baseline_control_n &
    S042_sufficient_followup_control_n) {
  S042_updated_control_counts$num_baseline_controls <-
    S042_max_baseline_controls
  S042_updated_control_counts$num_followup_controls <-
    floor(
      S042_adjusted_total_control_n_if_insufficient_baseline * S042_followup_GAD_site_visit_percent
    )
} else if (S042_sufficient_baseline_control_n &
           !S042_sufficient_followup_control_n) {
  S042_updated_control_counts$num_baseline_controls <-
    floor(
      S042_adjusted_total_control_n_if_insufficient_followup * S042_baseline_GAD_site_visit_percent
    )
  S042_updated_control_counts$num_followup_controls <-
    S042_max_followup_controls
} else {
  S042_updated_control_counts$num_baseline_controls <-
    S042_baseline_controls_desired
  S042_updated_control_counts$num_followup_controls <-
    S042_followup_controls_desired
}

#5.2719 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S042_baseline_control_sample <- S042_control_data[S042_control_data$subjectkey %in% S042_baseline_control_subs_only, ]
S042_followup_control_sample <- S042_control_data[S042_control_data$subjectkey %in% S042_followup_control_subs_only, ] 

#5.2720 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S042_num_baseline_controls_to_sample <- as.numeric(S042_updated_control_counts$num_baseline_controls)
S042_current_baseline_control_n <- as.numeric(length(unique(S042_baseline_control_sample$subjectkey)))

#5.2721 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S042_current_baseline_control_n == S042_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S042_current_baseline_control_n > S042_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S042_baseline_control_sample <-
      S042_baseline_control_sample[sample(nrow(S042_baseline_control_sample),
                                          S042_num_baseline_controls_to_sample),]
  } else if (S042_current_baseline_control_n < S042_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S042_baseline_control_sample_to_add <-
      S042_control_data[S042_control_data$subjectkey %in% S042_both_timepoint_sub_list &
                          S042_control_data$eventname == "baseline_year_1_arm_1",]
    S042_baseline_control_sample_to_add <-
      S042_baseline_control_sample_to_add[sample(
        nrow(S042_baseline_control_sample_to_add),
        S042_num_baseline_controls_to_sample - S042_current_baseline_control_n),]
    S042_baseline_control_sample <-
      rbind(S042_baseline_control_sample,
            S042_baseline_control_sample_to_add)
  }
}

#5.2722 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S042_num_followup_controls_to_sample <- as.numeric(S042_updated_control_counts$num_followup_controls)
S042_current_followup_control_n <- as.numeric(length(unique(S042_followup_control_sample$subjectkey)))
S042_subs_already_sampled <- S042_baseline_control_sample$subjectkey

#5.2723 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S042_current_followup_control_n == S042_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S042_current_followup_control_n > S042_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S042_followup_control_sample <-
      S042_followup_control_sample[sample(nrow(S042_followup_control_sample),
                                          S042_num_followup_controls_to_sample),]
  } else if (S042_current_followup_control_n < S042_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S042_followup_control_sample_to_add <-
      S042_control_data[S042_control_data$subjectkey %in% S042_both_timepoint_sub_list &
                          !S042_control_data$subjectkey %in% S042_subs_already_sampled &
                          S042_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S042_followup_control_sample_to_add <-
      S042_followup_control_sample_to_add[sample(
        nrow(S042_followup_control_sample_to_add),
        S042_num_followup_controls_to_sample - S042_current_followup_control_n),]
    S042_followup_control_sample <-
      rbind(S042_followup_control_sample,
            S042_followup_control_sample_to_add)
  }
}

#5.2724 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S042_baseline_control_sample$group <- rep("control")
S042_followup_control_sample$group <- rep("control")

#5.2725 Merge the baseline and followup S042 samples 
S042_control_sample <- merge(S042_baseline_control_sample, S042_followup_control_sample, all = TRUE)

#5.2726 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S042_cn_sample <- ifelse(length(unique(S042_control_sample$subjectkey)) == length(S042_control_sample$subjectkey), TRUE, FALSE)

#5.2727 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S042_control_sample, all = TRUE)


## S053 Resampling ##

#5.281 Filter the control data for just Site S053
S053_control_data <- raw_control_data %>%
  filter(site_name == "S053") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.282 Combine subject data at one site with n and pct values for that site
S053_merged_pct <- merge(site_visit_merged_pct, S053_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.283 Create a df for just baseline and a df for just followup data 
S053_control_baseline_only <- S053_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)

#5.284 save number of baseline only subjects and followup only subjects
S053_baseline_only_sub_n <- as.numeric(nrow(S053_control_baseline_only))

#5.285 save the subject IDs for the baseline only subjects and followup only subjects
S053_baseline_control_subs_only <- S053_control_baseline_only$subjectkey

#5.286 save number of unique control subjects at this site
S053_num_unique_control_subs <- as.numeric(length(unique(S053_merged_pct$subjectkey)))

#5.287 save subsetted dataframe for only S053 GAD subjects
S053_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S053") 

#5.288 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S053_control_baseline_only$group <- rep("control")

#5.289 Merge the baseline and followup S053 samples 
S053_control_sample <- S053_control_baseline_only %>% dplyr::select(-c(GAD_site_visit_n, GAD_site_n, GAD_site_pct, control_site_n, control_site_visit_n, control_site_pct, max_controls))

#5.2810 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S053_cn_sample <- ifelse(length(unique(S053_control_sample$subjectkey)) == length(S053_control_sample$subjectkey), TRUE, FALSE)

#5.2811 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S053_control_sample, all = TRUE)


## S065 Resampling ##

#5.291 Filter the control data for just Site S065
S065_control_data <- raw_control_data %>%
  filter(site_name == "S065") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.292 Combine subject data at one site with n and pct values for that site
S065_merged_pct <- merge(site_visit_merged_pct, S065_control_data, by = c("eventname","site_name"), all.X = TRUE)

#5.293 Create a df for just baseline and a df for just followup data 
S065_control_baseline_only <- S065_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)

#5.294 save number of baseline only subjects and followup only subjects
S065_baseline_only_sub_n <- as.numeric(nrow(S065_control_baseline_only))

#5.295 save the subject IDs for the baseline only subjects and followup only subjects
S065_baseline_control_subs_only <- S065_control_baseline_only$subjectkey

#5.296 save number of unique control subjects at this site
S065_num_unique_control_subs <- as.numeric(length(unique(S065_merged_pct$subjectkey)))

#5.297 save subsetted dataframe for only S065 GAD subjects
S065_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S065") 

#5.298 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S065_control_baseline_only$group <- rep("control")

#5.299 Merge the baseline and followup S065 samples 
S065_control_sample <- S065_control_baseline_only %>% dplyr::select(-c(GAD_site_visit_n, GAD_site_n, GAD_site_pct, control_site_n, control_site_visit_n, control_site_pct, max_controls))

#5.2910 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S065_cn_sample <- ifelse(length(unique(S065_control_sample$subjectkey)) == length(S065_control_sample$subjectkey), TRUE, FALSE)

#5.2911 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S065_control_sample, all = TRUE)


## S076 Resampling ##

#5.301 Filter the control data for just Site S076
S076_control_data <- raw_control_data %>%
  filter(site_name == "S076") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.302 Combine subject data at one site with n and pct values for that site
S076_merged_pct <- merge(site_visit_merged_pct, S076_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.303 Create a df for just baseline and a df for just followup data 
S076_control_baseline_only <- S076_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S076_control_followup_only <- S076_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.304 save number of baseline only subjects and followup only subjects
S076_baseline_only_sub_n <- as.numeric(nrow(S076_control_baseline_only))
S076_followup_only_sub_n <- as.numeric(nrow(S076_control_followup_only))

#5.305 save the subject IDs for the baseline only subjects and followup only subjects
S076_baseline_control_subs_only <- S076_control_baseline_only$subjectkey
S076_followup_control_subs_only <- S076_control_followup_only$subjectkey

#5.306 save number of unique control subjects at this site
S076_num_unique_control_subs <- as.numeric(length(unique(S076_merged_pct$subjectkey)))

#5.307 save subsetted dataframe for only S076 GAD subjects
S076_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S076") 

#5.308 calculate total number of subjects at site S076
S076_total_GAD_sub_n <- as.numeric(length(unique(S076_GAD_subs$subjectkey)))

#5.309 further subset the dataframe to only contain baseline scans or followup GAD scans
S076_baseline_GAD_data <- S076_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S076_followup_GAD_data <- S076_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.3010 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S076_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S076" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S076_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S076" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.3011 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S076_both_timepoint_sub_list <- S076_control_data$subjectkey[duplicated(S076_control_data$subjectkey) ]
S076_both_timepoint_sub_list_n <- as.numeric(length(unique(S076_both_timepoint_sub_list)))

#5.3012 Subset a dataframe of those subjects
S076_both_timepoint_subs <- S076_control_data[S076_control_data$subjectkey %in% S076_both_timepoint_sub_list, ]

#5.3013 Create vectors for the max number of baseline and max number of followup scans 
S076_max_baseline_controls <- as.numeric(S076_baseline_only_sub_n + S076_both_timepoint_sub_list_n)
S076_max_followup_controls <- as.numeric(S076_followup_only_sub_n + S076_both_timepoint_sub_list_n)

#5.3014 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S076_baseline_controls_desired <- as.numeric(S076_num_unique_control_subs * S076_baseline_GAD_site_visit_percent)
S076_followup_controls_desired <- as.numeric(S076_num_unique_control_subs * S076_followup_GAD_site_visit_percent)

#5.3015 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S076_sufficient_baseline_control_n <- if_else(S076_max_baseline_controls > S076_baseline_controls_desired, TRUE, FALSE)
S076_sufficient_followup_control_n <- if_else(S076_max_followup_controls > S076_followup_controls_desired, TRUE, FALSE)

#5.3016 Create an adjusted total control N to sample in the event either of the above checks are failed 
S076_adjusted_total_control_n_if_insufficient_baseline <- floor(S076_max_baseline_controls / S076_baseline_GAD_site_visit_percent)
S076_adjusted_total_control_n_if_insufficient_followup <- floor(S076_max_followup_controls / S076_followup_GAD_site_visit_percent)

#5.3017 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S076_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.3018 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S076_sufficient_baseline_control_n &
    S076_sufficient_followup_control_n) {
  S076_updated_control_counts$num_baseline_controls <-
    S076_max_baseline_controls
  S076_updated_control_counts$num_followup_controls <-
    floor(
      S076_adjusted_total_control_n_if_insufficient_baseline * S076_followup_GAD_site_visit_percent
    )
} else if (S076_sufficient_baseline_control_n &
           !S076_sufficient_followup_control_n) {
  S076_updated_control_counts$num_baseline_controls <-
    floor(
      S076_adjusted_total_control_n_if_insufficient_followup * S076_baseline_GAD_site_visit_percent
    )
  S076_updated_control_counts$num_followup_controls <-
    S076_max_followup_controls
} else {
  S076_updated_control_counts$num_baseline_controls <-
    S076_baseline_controls_desired
  S076_updated_control_counts$num_followup_controls <-
    S076_followup_controls_desired
}

#5.3019 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S076_baseline_control_sample <- S076_control_data[S076_control_data$subjectkey %in% S076_baseline_control_subs_only, ]
S076_followup_control_sample <- S076_control_data[S076_control_data$subjectkey %in% S076_followup_control_subs_only, ] 

#5.3020 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S076_num_baseline_controls_to_sample <- as.numeric(S076_updated_control_counts$num_baseline_controls)
S076_current_baseline_control_n <- as.numeric(length(unique(S076_baseline_control_sample$subjectkey)))

#5.3021 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S076_current_baseline_control_n == S076_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S076_current_baseline_control_n > S076_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S076_baseline_control_sample <-
      S076_baseline_control_sample[sample(nrow(S076_baseline_control_sample),
                                          S076_num_baseline_controls_to_sample),]
  } else if (S076_current_baseline_control_n < S076_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S076_baseline_control_sample_to_add <-
      S076_control_data[S076_control_data$subjectkey %in% S076_both_timepoint_sub_list &
                          S076_control_data$eventname == "baseline_year_1_arm_1",]
    S076_baseline_control_sample_to_add <-
      S076_baseline_control_sample_to_add[sample(
        nrow(S076_baseline_control_sample_to_add),
        S076_num_baseline_controls_to_sample - S076_current_baseline_control_n),]
    S076_baseline_control_sample <-
      rbind(S076_baseline_control_sample,
            S076_baseline_control_sample_to_add)
  }
}

#5.3022 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S076_num_followup_controls_to_sample <- as.numeric(S076_updated_control_counts$num_followup_controls)
S076_current_followup_control_n <- as.numeric(length(unique(S076_followup_control_sample$subjectkey)))
S076_subs_already_sampled <- S076_baseline_control_sample$subjectkey

#5.3023 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S076_current_followup_control_n == S076_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S076_current_followup_control_n > S076_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S076_followup_control_sample <-
      S076_followup_control_sample[sample(nrow(S076_followup_control_sample),
                                          S076_num_followup_controls_to_sample),]
  } else if (S076_current_followup_control_n < S076_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S076_followup_control_sample_to_add <-
      S076_control_data[S076_control_data$subjectkey %in% S076_both_timepoint_sub_list &
                          !S076_control_data$subjectkey %in% S076_subs_already_sampled &
                          S076_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S076_followup_control_sample_to_add <-
      S076_followup_control_sample_to_add[sample(
        nrow(S076_followup_control_sample_to_add),
        S076_num_followup_controls_to_sample - S076_current_followup_control_n),]
    S076_followup_control_sample <-
      rbind(S076_followup_control_sample,
            S076_followup_control_sample_to_add)
  }
}

#5.3024 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S076_baseline_control_sample$group <- rep("control")
S076_followup_control_sample$group <- rep("control")

#5.3025 Merge the baseline and followup S076 samples 
S076_control_sample <- merge(S076_baseline_control_sample, S076_followup_control_sample, all = TRUE)

#5.3026 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S076_cn_sample <- ifelse(length(unique(S076_control_sample$subjectkey)) == length(S076_control_sample$subjectkey), TRUE, FALSE)

#5.3027 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S076_control_sample, all = TRUE)


## S086 Resampling ##

#5.311 Filter the control data for just Site S086
S086_control_data <- raw_control_data %>%
  filter(site_name == "S086") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.312 Combine subject data at one site with n and pct values for that site
S086_merged_pct <- merge(site_visit_merged_pct, S086_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.313 Create a df for just baseline and a df for just followup data 
S086_control_baseline_only <- S086_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S086_control_followup_only <- S086_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.314 save number of baseline only subjects and followup only subjects
S086_baseline_only_sub_n <- as.numeric(nrow(S086_control_baseline_only))
S086_followup_only_sub_n <- as.numeric(nrow(S086_control_followup_only))

#5.315 save the subject IDs for the baseline only subjects and followup only subjects
S086_baseline_control_subs_only <- S086_control_baseline_only$subjectkey
S086_followup_control_subs_only <- S086_control_followup_only$subjectkey

#5.316 save number of unique control subjects at this site
S086_num_unique_control_subs <- as.numeric(length(unique(S086_merged_pct$subjectkey)))

#5.317 save subsetted dataframe for only S086 GAD subjects
S086_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S086") 

#5.318 calculate total number of subjects at site S086
S086_total_GAD_sub_n <- as.numeric(length(unique(S086_GAD_subs$subjectkey)))

#5.319 further subset the dataframe to only contain baseline scans or followup GAD scans
S086_baseline_GAD_data <- S086_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S086_followup_GAD_data <- S086_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.3110 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S086_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S086" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S086_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S086" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.3111 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S086_both_timepoint_sub_list <- S086_control_data$subjectkey[duplicated(S086_control_data$subjectkey) ]
S086_both_timepoint_sub_list_n <- as.numeric(length(unique(S086_both_timepoint_sub_list)))

#5.3112 Subset a dataframe of those subjects
S086_both_timepoint_subs <- S086_control_data[S086_control_data$subjectkey %in% S086_both_timepoint_sub_list, ]

#5.3113 Create vectors for the max number of baseline and max number of followup scans 
S086_max_baseline_controls <- as.numeric(S086_baseline_only_sub_n + S086_both_timepoint_sub_list_n)
S086_max_followup_controls <- as.numeric(S086_followup_only_sub_n + S086_both_timepoint_sub_list_n)

#5.3114 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S086_baseline_controls_desired <- as.numeric(S086_num_unique_control_subs * S086_baseline_GAD_site_visit_percent)
S086_followup_controls_desired <- as.numeric(S086_num_unique_control_subs * S086_followup_GAD_site_visit_percent)

#5.3115 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S086_sufficient_baseline_control_n <- if_else(S086_max_baseline_controls > S086_baseline_controls_desired, TRUE, FALSE)
S086_sufficient_followup_control_n <- if_else(S086_max_followup_controls > S086_followup_controls_desired, TRUE, FALSE)

#5.3116 Create an adjusted total control N to sample in the event either of the above checks are failed 
S086_adjusted_total_control_n_if_insufficient_baseline <- floor(S086_max_baseline_controls / S086_baseline_GAD_site_visit_percent)
S086_adjusted_total_control_n_if_insufficient_followup <- floor(S086_max_followup_controls / S086_followup_GAD_site_visit_percent)

#5.3117 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S086_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.3118 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S086_sufficient_baseline_control_n &
    S086_sufficient_followup_control_n) {
  S086_updated_control_counts$num_baseline_controls <-
    S086_max_baseline_controls
  S086_updated_control_counts$num_followup_controls <-
    floor(
      S086_adjusted_total_control_n_if_insufficient_baseline * S086_followup_GAD_site_visit_percent
    )
} else if (S086_sufficient_baseline_control_n &
           !S086_sufficient_followup_control_n) {
  S086_updated_control_counts$num_baseline_controls <-
    floor(
      S086_adjusted_total_control_n_if_insufficient_followup * S086_baseline_GAD_site_visit_percent
    )
  S086_updated_control_counts$num_followup_controls <-
    S086_max_followup_controls
} else {
  S086_updated_control_counts$num_baseline_controls <-
    S086_baseline_controls_desired
  S086_updated_control_counts$num_followup_controls <-
    S086_followup_controls_desired
}

#5.3119 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S086_baseline_control_sample <- S086_control_data[S086_control_data$subjectkey %in% S086_baseline_control_subs_only, ]
S086_followup_control_sample <- S086_control_data[S086_control_data$subjectkey %in% S086_followup_control_subs_only, ] 

#5.3120 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S086_num_baseline_controls_to_sample <- as.numeric(S086_updated_control_counts$num_baseline_controls)
S086_current_baseline_control_n <- as.numeric(length(unique(S086_baseline_control_sample$subjectkey)))

#5.3121 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S086_current_baseline_control_n == S086_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S086_current_baseline_control_n > S086_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S086_baseline_control_sample <-
      S086_baseline_control_sample[sample(nrow(S086_baseline_control_sample),
                                          S086_num_baseline_controls_to_sample),]
  } else if (S086_current_baseline_control_n < S086_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S086_baseline_control_sample_to_add <-
      S086_control_data[S086_control_data$subjectkey %in% S086_both_timepoint_sub_list &
                          S086_control_data$eventname == "baseline_year_1_arm_1",]
    S086_baseline_control_sample_to_add <-
      S086_baseline_control_sample_to_add[sample(
        nrow(S086_baseline_control_sample_to_add),
        S086_num_baseline_controls_to_sample - S086_current_baseline_control_n),]
    S086_baseline_control_sample <-
      rbind(S086_baseline_control_sample,
            S086_baseline_control_sample_to_add)
  }
}

#5.3122 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S086_num_followup_controls_to_sample <- as.numeric(S086_updated_control_counts$num_followup_controls)
S086_current_followup_control_n <- as.numeric(length(unique(S086_followup_control_sample$subjectkey)))
S086_subs_already_sampled <- S086_baseline_control_sample$subjectkey

#5.3123 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S086_current_followup_control_n == S086_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S086_current_followup_control_n > S086_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S086_followup_control_sample <-
      S086_followup_control_sample[sample(nrow(S086_followup_control_sample),
                                          S086_num_followup_controls_to_sample),]
  } else if (S086_current_followup_control_n < S086_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S086_followup_control_sample_to_add <-
      S086_control_data[S086_control_data$subjectkey %in% S086_both_timepoint_sub_list &
                          !S086_control_data$subjectkey %in% S086_subs_already_sampled &
                          S086_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S086_followup_control_sample_to_add <-
      S086_followup_control_sample_to_add[sample(
        nrow(S086_followup_control_sample_to_add),
        S086_num_followup_controls_to_sample - S086_current_followup_control_n),]
    S086_followup_control_sample <-
      rbind(S086_followup_control_sample,
            S086_followup_control_sample_to_add)
  }
}

#5.3124 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S086_baseline_control_sample$group <- rep("control")
S086_followup_control_sample$group <- rep("control")

#5.3125 Merge the baseline and followup S086 samples 
S086_control_sample <- merge(S086_baseline_control_sample, S086_followup_control_sample, all = TRUE)

#5.3126 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S086_cn_sample <- ifelse(length(unique(S086_control_sample$subjectkey)) == length(S086_control_sample$subjectkey), TRUE, FALSE)

#5.3127 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S086_control_sample, all = TRUE)


## S090 Resampling ##

#5.321 Filter the control data for just Site S090
S090_control_data <- raw_control_data %>%
  filter(site_name == "S090") %>%
  dplyr::select(subjectkey, eventname, site_name)

#5.322 Combine subject data at one site with n and pct values for that site
S090_merged_pct <- merge(site_visit_merged_pct, S090_control_data, by = c("eventname","site_name"), all.X = TRUE) 

#5.323 Create a df for just baseline and a df for just followup data 
S090_control_baseline_only <- S090_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
S090_control_followup_only <- S090_merged_pct %>%
  group_by(subjectkey) %>%
  filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)

#5.324 save number of baseline only subjects and followup only subjects
S090_baseline_only_sub_n <- as.numeric(nrow(S090_control_baseline_only))
S090_followup_only_sub_n <- as.numeric(nrow(S090_control_followup_only))

#5.325 save the subject IDs for the baseline only subjects and followup only subjects
S090_baseline_control_subs_only <- S090_control_baseline_only$subjectkey
S090_followup_control_subs_only <- S090_control_followup_only$subjectkey

#5.326 save number of unique control subjects at this site
S090_num_unique_control_subs <- as.numeric(length(unique(S090_merged_pct$subjectkey)))

#5.327 save subsetted dataframe for only S090 GAD subjects
S090_GAD_subs <- GAD_grouped_demo %>% 
  filter(site_name == "S090") 

#5.328 calculate total number of subjects at site S090
S090_total_GAD_sub_n <- as.numeric(length(unique(S090_GAD_subs$subjectkey)))

#5.329 further subset the dataframe to only contain baseline scans or followup GAD scans
S090_baseline_GAD_data <- S090_GAD_subs %>% 
  filter(eventname == "baseline_year_1_arm_1")
S090_followup_GAD_data <- S090_GAD_subs %>% 
  filter(eventname == "2_year_follow_up_y_arm_1")

#5.3210 save variable containing followup GAD vs followup percent distribution, and the controls needed based on those percents
S090_baseline_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S090" & eventname == "baseline_year_1_arm_1") %>% pull(GAD_site_pct)
S090_followup_GAD_site_visit_percent <- site_visit_merged_pct %>% filter(site_name == "S090" & eventname == "2_year_follow_up_y_arm_1") %>% pull(GAD_site_pct)

#5.3211 Create sub list of subjects that have both time points of data and a vector of the length of this list to know how many subs are shared and can be pulled from later on
S090_both_timepoint_sub_list <- S090_control_data$subjectkey[duplicated(S090_control_data$subjectkey) ]
S090_both_timepoint_sub_list_n <- as.numeric(length(unique(S090_both_timepoint_sub_list)))

#5.3212 Subset a dataframe of those subjects
S090_both_timepoint_subs <- S090_control_data[S090_control_data$subjectkey %in% S090_both_timepoint_sub_list, ]

#5.3213 Create vectors for the max number of baseline and max number of followup scans 
S090_max_baseline_controls <- as.numeric(S090_baseline_only_sub_n + S090_both_timepoint_sub_list_n)
S090_max_followup_controls <- as.numeric(S090_followup_only_sub_n + S090_both_timepoint_sub_list_n)

#5.3214 Create vectors for the desired number of baseline and followup control subjects based on the GAD distribution and the number of unique control subjects 
S090_baseline_controls_desired <- as.numeric(S090_num_unique_control_subs * S090_baseline_GAD_site_visit_percent)
S090_followup_controls_desired <- as.numeric(S090_num_unique_control_subs * S090_followup_GAD_site_visit_percent)

#5.3215 Check to see if there are sufficient baseline subjects or sufficient followup subjects to match the desired number of controls sampled (and determine the rate limiting factor)
S090_sufficient_baseline_control_n <- if_else(S090_max_baseline_controls > S090_baseline_controls_desired, TRUE, FALSE)
S090_sufficient_followup_control_n <- if_else(S090_max_followup_controls > S090_followup_controls_desired, TRUE, FALSE)

#5.3216 Create an adjusted total control N to sample in the event either of the above checks are failed 
S090_adjusted_total_control_n_if_insufficient_baseline <- floor(S090_max_baseline_controls / S090_baseline_GAD_site_visit_percent)
S090_adjusted_total_control_n_if_insufficient_followup <- floor(S090_max_followup_controls / S090_followup_GAD_site_visit_percent)

#5.3217 Create a data frame to store the number of baseline and followup control subjects to sample based upon the parameters above
S090_updated_control_counts <- data.frame(
  num_baseline_controls = 0,
  num_followup_controls = 0)

#5.3218 Create a loop that checks which of the rate limiting factors is present (insufficient baseline subjects or insufficient followup subjects to match the desired number of controls sampled) and adjusts the number of controls to sample according to the rate limiting factor & updated total control N to sample
if (!S090_sufficient_baseline_control_n &
    S090_sufficient_followup_control_n) {
  S090_updated_control_counts$num_baseline_controls <-
    S090_max_baseline_controls
  S090_updated_control_counts$num_followup_controls <-
    floor(
      S090_adjusted_total_control_n_if_insufficient_baseline * S090_followup_GAD_site_visit_percent
    )
} else if (S090_sufficient_baseline_control_n &
           !S090_sufficient_followup_control_n) {
  S090_updated_control_counts$num_baseline_controls <-
    floor(
      S090_adjusted_total_control_n_if_insufficient_followup * S090_baseline_GAD_site_visit_percent
    )
  S090_updated_control_counts$num_followup_controls <-
    S090_max_followup_controls
} else {
  S090_updated_control_counts$num_baseline_controls <-
    S090_baseline_controls_desired
  S090_updated_control_counts$num_followup_controls <-
    S090_followup_controls_desired
}

#5.3219 First sample the subjects who only have one timepoint of data. If the number of subjects sampled is greater than the required number, randomly drop the number of subjects required to match the expected number of controls. If the number of sampled subjects is less than the number expected, randomly sample the number of subjects required to meet the expected number of controls. If the number of controls sampled matches the required number of controls, exit the loop. 
S090_baseline_control_sample <- S090_control_data[S090_control_data$subjectkey %in% S090_baseline_control_subs_only, ]
S090_followup_control_sample <- S090_control_data[S090_control_data$subjectkey %in% S090_followup_control_subs_only, ] 

#5.3220 Now that we have the count of how many subjects to sample, randomly sample the required number of subjects from the baseline event. First, create a new subset of data containing only rows where `subjectkey` matches the baseline only controls vector, then set two variables to the expected and current number of rows in this subset.
S090_num_baseline_controls_to_sample <- as.numeric(S090_updated_control_counts$num_baseline_controls)
S090_current_baseline_control_n <- as.numeric(length(unique(S090_baseline_control_sample$subjectkey)))

#5.3221 Match the expected number of baseline controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of baseline controls is greater or less than the expected number
if (S090_current_baseline_control_n == S090_num_baseline_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of baseline controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S090_current_baseline_control_n > S090_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is greater than the expected number, randomly drop subjects
    S090_baseline_control_sample <-
      S090_baseline_control_sample[sample(nrow(S090_baseline_control_sample),
                                          S090_num_baseline_controls_to_sample),]
  } else if (S090_current_baseline_control_n < S090_num_baseline_controls_to_sample) {
    # If the current number of baseline controls is less than the expected number, randomly sample additional subjects
    S090_baseline_control_sample_to_add <-
      S090_control_data[S090_control_data$subjectkey %in% S090_both_timepoint_sub_list &
                          S090_control_data$eventname == "baseline_year_1_arm_1",]
    S090_baseline_control_sample_to_add <-
      S090_baseline_control_sample_to_add[sample(
        nrow(S090_baseline_control_sample_to_add),
        S090_num_baseline_controls_to_sample - S090_current_baseline_control_n),]
    S090_baseline_control_sample <-
      rbind(S090_baseline_control_sample,
            S090_baseline_control_sample_to_add)
  }
}

#5.3222 Now sample the followup controls. First, create a new subset of data containing only rows where `subjectkey` matches the followup only controls vector, then set two variables to the expected and current number of rows in this subset. Also, create a variable that logs the subjectkeys already sampled from the baseline group in the event both groups are sufficiently powered to be sampled (avoiding re-sampling subjects already selected in loop to come)
S090_num_followup_controls_to_sample <- as.numeric(S090_updated_control_counts$num_followup_controls)
S090_current_followup_control_n <- as.numeric(length(unique(S090_followup_control_sample$subjectkey)))
S090_subs_already_sampled <- S090_baseline_control_sample$subjectkey

#5.3223 Match the expected number of followup controls; print a message if they match. If they do not match, the code randomly samples or drops subjects from the data based on whether the current number of followup controls is greater or less than the expected number.
if (S090_current_followup_control_n == S090_num_followup_controls_to_sample) {
  # Exit if the number of controls matches the expected number
  message("Number of followup controls matched expected number")
} else {
  # If the number of controls does not match the expected number, randomly sample or drop subjects
  if (S090_current_followup_control_n > S090_num_followup_controls_to_sample) {
    # If the current number of followup controls is greater than the expected number, randomly drop subjects
    S090_followup_control_sample <-
      S090_followup_control_sample[sample(nrow(S090_followup_control_sample),
                                          S090_num_followup_controls_to_sample),]
  } else if (S090_current_followup_control_n < S090_num_followup_controls_to_sample) {
    # If the current number of followup controls is less than the expected number, randomly sample additional subjects
    S090_followup_control_sample_to_add <-
      S090_control_data[S090_control_data$subjectkey %in% S090_both_timepoint_sub_list &
                          !S090_control_data$subjectkey %in% S090_subs_already_sampled &
                          S090_control_data$eventname == "2_year_follow_up_y_arm_1",]
    S090_followup_control_sample_to_add <-
      S090_followup_control_sample_to_add[sample(
        nrow(S090_followup_control_sample_to_add),
        S090_num_followup_controls_to_sample - S090_current_followup_control_n),]
    S090_followup_control_sample <-
      rbind(S090_followup_control_sample,
            S090_followup_control_sample_to_add)
  }
}

#5.3224 Add a group column to the baseline and followup control data frames that labels those subjects as controls 
S090_baseline_control_sample$group <- rep("control")
S090_followup_control_sample$group <- rep("control")

#5.3225 Merge the baseline and followup S090 samples 
S090_control_sample <- merge(S090_baseline_control_sample, S090_followup_control_sample, all = TRUE)

#5.3226 Verify that only unique control subjects are in the final sample (should return TRUE)
unique_S090_cn_sample <- ifelse(length(unique(S090_control_sample$subjectkey)) == length(S090_control_sample$subjectkey), TRUE, FALSE)

#5.3227 Create a Dataframe containing the site-specific control sample to be merged with other sites
control_sample <- merge(control_sample, S090_control_sample, all = TRUE)


## Data Integrity Checks ##
#1. Create a dataframe containing the distribution of  sampled controls at each site-event combination (both N and % distribution)
final_site_visit_control_count <- control_sample %>% 
  group_by(site_name, eventname) %>%
  summarize(control_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(control_site_n = sum(control_site_visit_n),
         control_site_pct = (control_site_visit_n / control_site_n))

#2. Create a dataframe comparing the distribution of GAD subjects and sampled controls at each site-event combination (both N and % distribution), and add a column that tests whether the GAD distribution matches the sampled control subject distribution
final_comparitive_count <- merge(site_visit_GAD_count, final_site_visit_control_count, by = c("site_name", "eventname"), all = TRUE)
final_comparitive_count$control_site_pct <- round(final_comparitive_count$control_site_pct, digits = 2)
final_comparitive_count$distribution_converged <- ifelse(final_comparitive_count$control_site_pct == final_comparitive_count$GAD_site_pct, "Success", "Debug")

#3. Check whether the sampled control distributions at each site-event combination and whether that matches the distribution of the GAD subjects, and print a message detailing whether those were successful or not
print(if_else(length(unique(control_sample$subjectkey)) == length(control_sample$subjectkey) && all(final_comparitive_count$converged == "Success") , "Control sample distribution is correct & contains only unique subjects", "Error: Controls not sampled correctly"))

#4. Check whether the sampled control distributions at each site-event combination and whether that matches the distribution of the GAD subjects, and store a TRUE/FALSE value that signifies whether those were successful or not
final_comparitive_count$final_cn_sample_unique <- if_else(length(unique(control_sample$subjectkey)) == length(control_sample$subjectkey) && all(final_comparitive_count$converged == "Success") , TRUE, FALSE)
final_cn_sample_unique <- if_else(length(unique(control_sample$subjectkey)) == length(control_sample$subjectkey) && all(final_comparitive_count$converged == "Success") , TRUE, FALSE)

#5. If the controls were sampled successfully, print a message detailing this; if they were not sampled correctly, print the subjects that were duplicated 
if_else(final_cn_sample_unique == TRUE, print("NA; Controls Sampled Successfully", print(control_sample[duplicated(control_sample$subjectkey), "subjectkey"])))


## Output ##

#1. Save a csv file containing the site-visit comparison counts
write.csv(final_comparitive_count, "./data_processed/main_analysis/resampled_gad_hc_comparitive_count.csv", row.names = FALSE)

#2. Save a csv file containing the final control sample dataframe
write.csv(control_sample, "./data_processed/main_analysis/resampled_hc_sample.csv", row.names = FALSE)

#3. Write a csv file containing the GAD group data
write.csv(GAD_grouped_demo, "./data_processed/main_analysis/gad_group.csv", row.names = FALSE)

#4. Write a csv file containing the QC'd + QA'd imaging data
write.csv(clinical_imaging_data, "./data_processed/main_analysis/subset_qcd_imaging_data.csv", row.names = FALSE)
