## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(tidyr)
library(broom)
library(knitr)
library(ggplot2)
library(stringr)
options(digits = 8, scipen = 999) 

# Read in required data 
# ABCD release 5.1 parent report KSADS-COMP data
ABCD_KSADS_release_5.1_parent_report <- read.csv("./data_raw/mh_p_ksads_ss.csv")

# ABCD release 5.1 youth report KSADS-COMP data
ABCD_KSADS_release_5.1_youth_report <- read.csv("./data_raw/mh_y_ksads_ss.csv")

# Cleaned and QC'd imaging data
cleaned_qcd_imaging_data <- read.csv("./data_processed/rsfMRI_data_qcd_unsubset.csv")

# Diagnosis group labels
GAD_HC_group <- read.csv("./data_processed/GAD_HC_resampled_merged_groups.csv") %>% 
  dplyr::select(c(src_subject_id, eventname, group))

# Original HC group (non-resampled)
original_HC_subjects_not_resampled <- read.csv("./data_processed/GAD_HC_subjects_not_resampled.csv") %>% 
  filter(group == "HC") %>% 
  dplyr::select(c(src_subject_id, eventname, group))


## Data Wrangling ## 

#1. Prep the KSADS Data
#1.11 Full join the parent and youth report KSADS dataframes
ABCD_KSADS_release_5.1_data <- full_join(ABCD_KSADS_release_5.1_parent_report, ABCD_KSADS_release_5.1_youth_report)

#1.2 Clean the parent and youth report data
#1.21 Retain only the columns of interest to creating the GAD group(s) and assessing comorbidity
ABCD_KSADS_release_5.1_GAD_data <- ABCD_KSADS_release_5.1_data %>% 
  dplyr::select(c(src_subject_id, eventname, ksads_10_869_p, ksads_10_869_t, ksads_8_863_p, ksads_8_863_t, ksads_7_861_p, ksads_7_861_t, ksads_1_840_p, ksads_1_840_t))

#1.22 Retain only the baseline and 2 year follow up data
ABCD_KSADS_release_5.1_GAD_data_filtered <- ABCD_KSADS_release_5.1_GAD_data %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#1.23 Convert any instances of NA, empty string, 555, or 999 to NA; and any instances of 888 to 0
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_filtered %>%
  mutate_all(~ ifelse(. %in% c(999, 555, ""), NA, ifelse(. == 888, 0, .)))

#1.3 Bring in HC + GAD diagnosis labels
#1.31 Join the HC + GAD labels into the KSADS data
ABCD_KSADS_release_5.1_data_for_merging <- left_join(ABCD_KSADS_release_5.1_data_for_merging, GAD_HC_group)

#1.32 Create a new GAD diangostic label determining the source of the GAD diagnosis for this set of descriptive stats 
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  mutate(
    GAD = if_else(
      replace_na(ksads_10_869_p, 0) == 1 | replace_na(ksads_10_869_t, 0) == 1, 1, 0
    ),
    GAD_Source = case_when(
      replace_na(ksads_10_869_p, 0) == 1 & replace_na(ksads_10_869_t, 0) == 1 ~ "Concordant",
      replace_na(ksads_10_869_p, 0) == 1 ~ "Parent",
      replace_na(ksads_10_869_t, 0) == 1 ~ "Youth",
      TRUE ~ "None"
    ))

#1.4 Create comorbidities of interest
#1.41 Social Anxiety Disorder
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  mutate(
    Social_Anxiety_Disorder = if_else(
      replace_na(ksads_8_863_p, 0) == 1 | replace_na(ksads_8_863_t, 0) == 1, 1, 0
    ),
    Social_Anxiety_Disorder_Source = case_when(
      replace_na(ksads_8_863_p, 0) == 1 & replace_na(ksads_8_863_t, 0) == 1 ~ "Concordant",
      replace_na(ksads_8_863_p, 0) == 1 ~ "Parent",
      replace_na(ksads_8_863_t, 0) == 1 ~ "Youth",
      TRUE ~ "None"
    ))

#1.42 Separation Anxiety Disorder
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  mutate(
    Separation_Anxiety_Disorder = if_else(
      replace_na(ksads_7_861_p, 0) == 1 | replace_na(ksads_7_861_t, 0) == 1, 1, 0
    ),
    Separation_Anxiety_Disorder_Source = case_when(
      replace_na(ksads_7_861_p, 0) == 1 & replace_na(ksads_7_861_t, 0) == 1 ~ "Concordant",
      replace_na(ksads_7_861_p, 0) == 1 ~ "Parent",
      replace_na(ksads_7_861_t, 0) == 1 ~ "Youth",
      TRUE ~ "None"
    ))

#1.43 MDD
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  mutate(
    MDD = if_else(
      replace_na(ksads_1_840_p, 0) == 1 | replace_na(ksads_1_840_t, 0) == 1, 1, 0
    ),
    MDD_Source = case_when(
      replace_na(ksads_1_840_p, 0) == 1 & replace_na(ksads_1_840_t, 0) == 1 ~ "Concordant",
      replace_na(ksads_1_840_p, 0) == 1 ~ "Parent",
      replace_na(ksads_1_840_t, 0) == 1 ~ "Youth",
      TRUE ~ "None"
    ))

#1.5 Retain columns of interest
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>% 
  dplyr::select(c(src_subject_id, eventname, group, GAD_Source, Social_Anxiety_Disorder, Social_Anxiety_Disorder_Source, Separation_Anxiety_Disorder, Separation_Anxiety_Disorder_Source, MDD, MDD_Source))

#1.6 Clean subjects in the KSADS data to merge
#1.61 Remove subjects whos rsfMRI data is incomplete or did not pass QA/QC
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  semi_join(cleaned_qcd_imaging_data, by = c("src_subject_id", "eventname"))

#1.62 Remove rows where the group variable contains "HC", keeping NAs
ABCD_KSADS_release_5.1_data_for_merging <- ABCD_KSADS_release_5.1_data_for_merging %>%
  filter(!grepl("HC", group, ignore.case = TRUE) | is.na(group))

#1.63 For the supplementary analyses, we will be prioritizing the parent report diagnosis data. As such, retain only subjects who have a parent reported diagnosis of each of the respected comorbidities
KSADS_data_for_grouping_filtered <- ABCD_KSADS_release_5.1_data_for_merging %>%
  filter(
    !(
      # Check if any column contains "Youth"
      if_any(all_of(c(
        "GAD_Source",
        "Social_Anxiety_Disorder_Source",
        "Separation_Anxiety_Disorder_Source",
        "MDD_Source"
      )), ~ .x == "Youth") &
        
        # Ensure no column contains "Parent" or "Concordant"
        if_all(all_of(c(
          "GAD_Source",
          "Social_Anxiety_Disorder_Source",
          "Separation_Anxiety_Disorder_Source",
          "MDD_Source"
        )), ~ !(.x %in% c("Parent", "Concordant")))
    )
  )

#1.64 Remove subjects who have no GAD or comorbid diagnoses of interest
KSADS_data_for_grouping_filtered <- KSADS_data_for_grouping_filtered %>% 
  filter(group == "GAD" | Social_Anxiety_Disorder == 1 | Separation_Anxiety_Disorder == 1 | MDD == 1)


#2. Create the comorbidity groups
#2.1 Establish functions to aid in the formation of the GAD and comorbidity/other clinical sam
#2.11 Create a function to randomly assign a subject with multiple diagnoses to a timepoint
random_sample_to_timepoint <- function(df, timepoint_var, eventname_vals = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1")) {
  df %>%
    group_by(src_subject_id) %>%
    mutate(
      # Randomly assign one of the timepoints for subjects with multiple timepoints
      random_timepoint = if_else(
        eventname %in% eventname_vals,
        sample(eventname_vals, 1),
        eventname
      )
    ) %>%
    ungroup()
}

#2.12 Helper function to check if the diagnosis is either 'Parent' or 'Concordant'
is_parent_or_concordant <- function(diagnosis_column) {
  diagnosis_column %in% c("Parent", "Concordant")
}

#2.21 Subset GAD subjects with any comorbidities (parent or concordant reported diagnoses only)
GAD_with_comorbidities <- KSADS_data_for_grouping_filtered %>%
  filter(group == "GAD" &
      # Check if GAD diagnosis is Parent or Concordant
      is_parent_or_concordant(GAD_Source) &
      (
        # Check if Social Anxiety is Parent or Concordant
        is_parent_or_concordant(Social_Anxiety_Disorder_Source) |
          # Check if Separation Anxiety is Parent or Concordant
          is_parent_or_concordant(Separation_Anxiety_Disorder_Source) |
          # Check if MDD diagnosis is Parent or Concordant
          is_parent_or_concordant(MDD_Source)) &
      (Social_Anxiety_Disorder == 1 | Separation_Anxiety_Disorder == 1 | MDD == 1)) %>%
  mutate(comorbidity_group = "GAD_with_comorbidities")

#2.22 Subset GAD subjects only (no comorbidities, parent or concordant reported diagnoses only)
GAD_only <- KSADS_data_for_grouping_filtered %>%
  filter(group == "GAD" & 
           Social_Anxiety_Disorder == 0 & 
           Separation_Anxiety_Disorder == 0 & 
           MDD == 0 &
           is_parent_or_concordant(GAD_Source)) %>% 
  mutate(comorbidity_group = "GAD_Only")

#2.23 Subset Social Anxiety without GAD, Separation Anxiety, or MDD (parent or concordant reported diagnoses only)
Social_Anxiety_only <- KSADS_data_for_grouping_filtered %>%
  filter(Social_Anxiety_Disorder == 1 & 
           group != "GAD" & 
           Separation_Anxiety_Disorder == 0 & 
           MDD == 0 &
           is_parent_or_concordant(Social_Anxiety_Disorder_Source)) %>%
  mutate(comorbidity_group = "Social_Anxiety_Only")

#2.24 Subset Separation Anxiety without GAD (parent or concordant reported diagnoses only)
Separation_Anxiety_only <- KSADS_data_for_grouping_filtered %>%
  filter(Separation_Anxiety_Disorder == 1 & 
           group != "GAD" & 
           Social_Anxiety_Disorder == 0 & 
           MDD == 0 &
           is_parent_or_concordant(Separation_Anxiety_Disorder_Source)) %>%
  mutate(comorbidity_group = "Separation_Anxiety_Only")

#2.25 Subset MDD without GAD (parent or concordant reported diagnoses only)
MDD_only <- KSADS_data_for_grouping_filtered %>%
  filter(MDD == 1 & 
           group != "GAD" & 
           Social_Anxiety_Disorder == 0 & 
           Separation_Anxiety_Disorder == 0 &
           is_parent_or_concordant(MDD_Source)) %>%
  mutate(comorbidity_group = "MDD_Only")

#2.26 Subset subjects with multiple comorbidities (having any combination of comorbid conditions, with parent or concordant diagnoses only)
subjects_w_multiple_comorbidities <- KSADS_data_for_grouping_filtered %>%
  filter((is_parent_or_concordant(Social_Anxiety_Disorder_Source) & 
            is_parent_or_concordant(Separation_Anxiety_Disorder_Source)) |
           (is_parent_or_concordant(Social_Anxiety_Disorder_Source) & 
              is_parent_or_concordant(MDD_Source)) |
           (is_parent_or_concordant(Separation_Anxiety_Disorder_Source) & 
              is_parent_or_concordant(MDD_Source))) %>%
  mutate(comorbidity_group = "Multiple_Comorbidities")

#2.3 Randomly assign subjects_w_multiple_comorbidities to one of the diagnostic groups
subjects_w_multiple_comorbidities <- subjects_w_multiple_comorbidities %>%
  mutate(comorbidity_group = sample(
    c("GAD_with_comorbidities", "Social_Anxiety_Only", "Separation_Anxiety_Only", "MDD_Only"),
    size = n(),
    replace = TRUE))



#2.4 Combine all groups into one dataframe
final_groups <- bind_rows(GAD_with_comorbidities, 
                          GAD_only, 
                          Social_Anxiety_only, 
                          Separation_Anxiety_only, 
                          MDD_only, 
                          subjects_w_multiple_comorbidities)

#2.5 Now we will apply random sampling for continuous GAD subjects to either timepoint
final_groups <- final_groups %>%
  random_sample_to_timepoint(timepoint_var = "eventname")


#3. Replace the existing HC with the non resampled HC subjects so that they can be resampled again according to the comorbidity groups
#3.1 Prepare the original HC group dataframe for merging
original_HC_subjects_for_merging <- original_HC_subjects_not_resampled %>% 
  mutate(GAD_Source = rep("None"),
         Social_Anxiety_disorder = 0,
         Social_Anxiety_Disorder_Source = rep("None"),
         Separation_Anxiety_Disorder = 0,
         Separation_Anxiety_Disorder_Source = rep("None"),
         MDD = 0,
         MDD_Source = rep("None"))

#3.2 Merge the non-resampled HC subjects into the updated supplemental (comorbidity) sample for resampling
supplemental_comorbidity_analysis_sample_not_resampled <- full_join( , original_HC_subjects_for_merging)


## Output ## 

#1. Write the comorbidity sample + non-resampled HC subject dataframe as a csv
write.csv(supplemental_comorbidity_analysis_sample_not_resampled, row.names = FALSE)