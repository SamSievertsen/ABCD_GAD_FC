# Load in Necessary Packages
library(knitr)
library(kableExtra)
library(dplyr)
library(tidyr)


## Set Up ##

# Read in the ksads-comp data without the desc row data 
KSADS_Data_no_desc <- read.csv("./data_raw/ksads_raw_data.csv") 
KSADS_Data_no_desc <- KSADS_Data_no_desc[-1,]

# Read in the GAD vs HC grouped connectivity analysis data
site_visit_analysis_data <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/site_visit_analysis_data.csv")

# Read in the release 5.1 data
ABCD_KSADS_5.1_release_parent <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/v5.1_mh_p_ksads_ss.csv")
ABCD_KSADS_5.1_release_youth <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/v5.1_mh_y_ksads_ss.csv")
ABCD_4.0_Youth_Report_KSADS <- read.table("/Users/samsievertsen/Desktop/ABCC_Package_1203705_Tabulated-Data/abcd_ksad501.txt")


## Data Wrangling ##

#1. Youth report
#1.1 Version 4.0
#1.11 Set the first row as column names
colnames(ABCD_4.0_Youth_Report_KSADS) <- ABCD_4.0_Youth_Report_KSADS[1, ]

#1.12 Remove the first two rows
ABCD_4.0_Youth_Report_KSADS <- ABCD_4.0_Youth_Report_KSADS[-c(1, 2), ]

#1.13 Retain only baseline and 2 year follow up
ABCD_4.0_Youth_Report_KSADS <- ABCD_4.0_Youth_Report_KSADS %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#1.14 Retain only the columns of interest
ABCD_4.0_Youth_Report_KSADS <- ABCD_4.0_Youth_Report_KSADS %>% 
  dplyr::select(c(src_subject_id, subjectkey, interview_age, sex, eventname, ksads_10_869_t, ksads_8_863_t, ksads_7_861_t))
  
#1.2 Version 5.1
#1.21 Retain only baseline and 2 year follow up
ABCD_KSADS_5.1_release_youth <- ABCD_KSADS_5.1_release_youth %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#1.22 Only retain columns of interest, AKA required demo data and MDD diagnoses
ABCD_KSADS_5.1_release_youth <- ABCD_KSADS_5.1_release_youth %>% 
  dplyr::select(c(src_subject_id, eventname, ksads_1_840_t))

#1.3 Merged youth report data
#1.31 Join the release 4.0 and 5.1 data
Merged_v4_v5_KSADS_Youth_Data <- full_join(ABCD_4.0_Youth_Report_KSADS, ABCD_KSADS_5.1_release_youth)

#1.32 For merging purposes, make the interview age column numeric
Merged_v4_v5_KSADS_Youth_Data$interview_age <- as.numeric(Merged_v4_v5_KSADS_Youth_Data$interview_age)

#2. Parent report
#2.1 Release 5.1 
#2.11 Subset the release 5.1 data for just current MDD 
ABCD_KSADS_5.1_release_parent_subset <- ABCD_KSADS_5.1_release_parent %>% 
  dplyr::select(c(src_subject_id, eventname, ksads_1_840_p))

#2.2 Release 4.0 
#2.21 Subset columns of interest in the V4 data, ensuring the removal of the current MDD column too (since V5 data contains current MDD)
KSADS_Data_no_desc_subset <- KSADS_Data_no_desc %>% 
  dplyr::select(c(src_subject_id, subjectkey, interview_age, sex, eventname, ksads_10_869_p, ksads_8_863_p, ksads_7_861_p))

#2.3 Merged data 
#2.31 Merge the current MDD data into the KSADS dataframe
Merged_v4_v5_KSADS_Parent_Data <- left_join(KSADS_Data_no_desc_subset, ABCD_KSADS_5.1_release_parent_subset)

#2.32 Subset the KSADS dataframe to only include diagnoses + data of interest
Merged_v4_v5_KSADS_Parent_Data <- Merged_v4_v5_KSADS_Parent_Data %>% 
  dplyr::select(c(subjectkey, interview_age, sex, eventname, ksads_10_869_p, ksads_8_863_p, ksads_7_861_p, ksads_1_840_p))

#2.33 Change the interview age data to numeric type
Merged_v4_v5_KSADS_Parent_Data$interview_age <- as.numeric(KSADS_Data_no_desc$interview_age)

#2.34 Subset the assessment timepoints of interest for the merged KSADS data
Merged_v4_v5_KSADS_Parent_Data <- Merged_v4_v5_KSADS_Parent_Data %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")


#3. Merged Youth and Parent Report
#3.11 Join the release 4.0 + 5.1 Youth and Parent Data
Merged_v4_v5_KSADS_Data <- full_join(Merged_v4_v5_KSADS_Parent_Data, Merged_v4_v5_KSADS_Youth_Data)

#3.12 Ensure all NA, missing, 999, or 555 are converted to NA, and all 888 are converted to 0
Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>%
  mutate_all(~ ifelse(. %in% c(999, 555, ""), NA, ifelse(. == 888, 0, .)))

#3.131 Remove duplicated ID column
Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>% 
  dplyr::select(-subjectkey)

#3.132 Reorder columns
Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>% 
  dplyr::select(c(src_subject_id, everything(.)))

#3.14 Make all KSADS diagnosis columns numeric type
Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>%
  mutate_at(vars(5:12), as.numeric)

#3.2 Create combined parent OR youth reported current diagnosis columns
#3.21 Social Anxiety Disorder
Merged_v4_v5_KSADS_Data$Social_Anxiety_Disorder <- if_else(
  replace_na(Merged_v4_v5_KSADS_Data$ksads_8_863_p, 0) == 1 | 
    replace_na(Merged_v4_v5_KSADS_Data$ksads_8_863_t, 0) == 1, 1, 0)

#3.22 Separation Anxiety Disorder
Merged_v4_v5_KSADS_Data$Separation_Anxiety_Disorder <- if_else(
  replace_na(Merged_v4_v5_KSADS_Data$ksads_7_861_p, 0) == 1 | 
    replace_na(Merged_v4_v5_KSADS_Data$ksads_7_861_t, 0) == 1, 1, 0)

#3.23 MDD
Merged_v4_v5_KSADS_Data$MDD <- if_else(
  replace_na(Merged_v4_v5_KSADS_Data$ksads_1_840_p, 0) == 1 | 
    replace_na(Merged_v4_v5_KSADS_Data$ksads_1_840_t, 0) == 1, 1, 0)

#3.24 GAD 
Merged_v4_v5_KSADS_Data$GAD <- if_else(
  replace_na(Merged_v4_v5_KSADS_Data$ksads_10_869_p, 0) == 1 | 
    replace_na(Merged_v4_v5_KSADS_Data$ksads_10_869_t, 0) == 1, 1, 0)

#3.3 Retain only columns of interest for stats purposes
Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>% 
  dplyr::select(c(src_subject_id, eventname, Social_Anxiety_Disorder, Separation_Anxiety_Disorder, MDD, GAD))

#3.4 Determine who, according to the criteria in the original GAD group analysis (complete rsfMRI data and passes rsfMRI QC/QA), can be kept in the sample for comparative analysis
#3.41 Identify the common columns between the current iteration of the analysis data and the rsfMRI data
shared_rsfMRI_analysis_data_columns <- intersect(names(ABCD_rsfMRI_Data_no_desc), names(site_visit_analysis_data))

#3.42 Subset the ABCD_rsfMRI_Data_no_desc data to retain only the shared columns
ABCD_rsfMRI_Data_matched <- ABCD_rsfMRI_Data_no_desc[, shared_rsfMRI_analysis_data_columns]

#3.43 Merge ABCD_rsfMRI_Data_matched with only the necessary columns from recommended_imaging
ABCD_rsfMRI_Data_merged <- merge(
  ABCD_rsfMRI_Data_matched,
  recommended_imaging[, c("subjectkey", "eventname", "imgincl_rsfmri_include")],
  by = c("subjectkey", "eventname"), all.x = TRUE)

#3.44 Filter the data to keep only rows where 'imgincl_rsfmri_include' is 1 (passed quality control)
ABCD_rsfMRI_Data_filtered <- ABCD_rsfMRI_Data_merged[ABCD_rsfMRI_Data_merged$imgincl_rsfmri_include == 1, ]

#3.45 Drop the 'imgincl_rsfmri_include' column
ABCD_rsfMRI_Data_filtered <- ABCD_rsfMRI_Data_filtered[, !names(ABCD_rsfMRI_Data_filtered) %in% "imgincl_rsfmri_include"]

#3.46 Keep only subjects who have complete rsfMRI data
ABCD_rsfMRI_Data_cleaned <- ABCD_rsfMRI_Data_filtered[
  complete.cases(ABCD_rsfMRI_Data_filtered[, 6:99]) & 
    !apply(ABCD_rsfMRI_Data_filtered[, 6:99] == "", 1, any),]

#3.47 Create a df containing the subject IDs and eventnames of participants who meet imaging criteria 
ABCD_subjects_meeting_imaging_criteria <- ABCD_rsfMRI_Data_cleaned %>% 
  dplyr::select(c(subjectkey, eventname))

#3.5 In the data to be analyzed for summary stats, retain only the subjects who meet criteria to be included in future imaging analyses
Filtered_Merged_v4_v5_KSADS_Data <- Merged_v4_v5_KSADS_Data %>%
  inner_join(ABCD_subjects_meeting_imaging_criteria, by = c("src_subject_id" = "subjectkey", "eventname"))


## Descriptive Stats ##
#1. Create the whole sample summary table
Whole_Sample_Dx_Summary <- Filtered_Merged_v4_v5_KSADS_Data %>%
  group_by(eventname) %>%
  summarize(
    # Total number of subjects
    N_Subjects = n(),
    # GAD diagnosis
    N_GAD_Dx = sum(as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_GAD_Dx = (N_GAD_Dx / N_Subjects * 100),
    # Social Anxiety diagnosis
    N_Social_Anxiety_Dx = sum(as.numeric(Social_Anxiety_Disorder) == 1, na.rm = TRUE),
    Percent_Social_Anxiety_Dx = (N_Social_Anxiety_Dx / N_Subjects * 100),
    # Separation Anxiety diagnosis
    N_Separation_Anxiety_Dx = sum(as.numeric(Separation_Anxiety_Disorder) == 1, na.rm = TRUE),
    Percent_Separation_Anxiety_Dx = (N_Separation_Anxiety_Dx / N_Subjects * 100),
    # MDD diagnosis
    N_MDD = sum(as.numeric(MDD) == 1, na.rm = TRUE),
    Percent_MDD = (N_MDD / N_Subjects * 100),
    # Comorbidity: Social Anxiety + GAD (using GAD as reference point)
    N_GAD_with_Social_Anxiety_Comorbidity = sum(as.numeric(Social_Anxiety_Disorder) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_GAD_with_Social_Anxiety_Comorbidity = (N_GAD_with_Social_Anxiety_Comorbidity / N_GAD_Dx * 100),
    # Comorbidity: Separation Anxiety + GAD (using GAD as reference point)
    N_GAD_with_Separation_Anxiety_Comorbidity = sum(as.numeric(Separation_Anxiety_Disorder) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_GAD_with_Separation_Anxiety_Comorbidity = (N_GAD_with_Separation_Anxiety_Comorbidity / N_GAD_Dx * 100),
    # Comorbidity: MDD + GAD (using GAD as reference point)
    N_GAD_with_MDD_Comorbidity = sum(as.numeric(MDD) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_GAD_with_MDD_Comorbidity = (N_GAD_with_MDD_Comorbidity / N_GAD_Dx * 100),
    # Comorbidity: Social Anxiety + GAD (using Social Anxiety as reference point)
    N_Social_Anxiety_with_GAD_Comorbidity = sum(as.numeric(Social_Anxiety_Disorder) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_Social_Anxiety_with_GAD_Comorbidity = (N_Social_Anxiety_with_GAD_Comorbidity / N_Social_Anxiety_Dx * 100),
    # Comorbidity: Separation Anxiety + GAD (using Separation Anxiety as reference point)
    N_Separation_Anxiety_with_GAD_Comorbidity = sum(as.numeric(Separation_Anxiety_Disorder) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_Separation_Anxiety_with_GAD_Comorbidity = (N_Separation_Anxiety_with_GAD_Comorbidity/ N_Separation_Anxiety_Dx * 100),
    # Comorbidity: MDD + GAD (using MDD as reference point)
    N_MDD_with_GAD_Comorbidity = sum(as.numeric(MDD) == 1 & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_MDD_with_GAD_Comorbidity = (N_MDD_with_GAD_Comorbidity / N_MDD * 100),
    # Any comorbidity (with GAD): At least one non-GAD diagnosis with GAD (unduplicated)
    N_Any_Comorbidity_with_GAD = sum((as.numeric(Social_Anxiety_Disorder) == 1 | as.numeric(Separation_Anxiety_Disorder) == 1 | as.numeric(MDD) == 1) & as.numeric(GAD) == 1, na.rm = TRUE),
    Percent_Any_Comorbidity_with_GAD = (N_Any_Comorbidity_with_GAD / N_GAD_Dx * 100),
    # Subjects with only GAD (no comorbidity)
    N_Only_GAD = N_GAD_Dx - N_Any_Comorbidity_with_GAD,
    Percent_Only_GAD = (N_Only_GAD / N_GAD_Dx * 100))

#2. Create a pivoted version of the whole sample summary table
Whole_Sample_Long_Dx_Summary <- Whole_Sample_Dx_Summary %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Percent_"), 
    names_to = c(".value", "Diagnosis"), 
    names_pattern = "(N|Percent)_(.*)") %>%
  mutate(Diagnosis = gsub("_", " ", Diagnosis))

#3. Create a cleaned version of the summary table with event names as columns and diagnoses as rows
Whole_Sample_Formatted_Dx_Summary <- Whole_Sample_Long_Dx_Summary %>%
  select(eventname, Diagnosis, N, Percent) %>%
  pivot_wider(names_from = eventname, values_from = c(N, Percent), names_sep = "_") %>%
  rename_with(~ gsub("year_1_arm_1|y_arm_1", "", .)) %>%
  rename_with(~ gsub("_", " ", .)) %>%
  dplyr::select(c("Diagnosis", `N baseline `, `Percent baseline `, `N 2 year follow up `, `Percent 2 year follow up `)) %>% 
  kable("html", escape = FALSE) %>%
  kable_styling(full_width = FALSE, position = "left") %>%
  column_spec(1, bold = TRUE)

#3.1 Print the summary table for extraction
Whole_Sample_Formatted_Dx_Summary
