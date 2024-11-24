## Set up ##

# Load packages for loading, wrangling, mutating, and visualizing data
library(dplyr)
library(DT)
library(lubridate)
library(writexl)
library(readxl)
library(stringr)

# Read in relevant data
abcd_ksad_data_amended_dx <- read.csv("./data_raw/abcd_ksad_data_amended_dx.csv")


## Create control data set for no diagnoses at baseline (excluding diagnoses highlighted in green on SP). See how many time points are available for clinical data, make additional control group containing all time points (maintain visit label variable) ##

#1. Control data set for no diagnoses at baseline
#1.1 Filter the KSADS data frame for just baseline visits
baseline_only_KSADS_data <- subset(abcd_ksad_data_amended_dx, eventname == "baseline_year_1_arm_1")

#1.2 Create a control data frame containing subjects not positive for any specified diagnoses @ baseline (highlighted dx already excluded)
control_data_0_dx_at_baseline <- data.frame(subset(baseline_only_KSADS_data, ksads_6_860_p != '1' & ksads_6_859_p != '1' & ksads_19_892_p != '1' & ksads_19_891_p != '1' & ksads_13_932_p != '1' & ksads_13_933_p != '1' & ksads_13_934_p != '1' & ksads_13_931_p != '1' & ksads_13_929_p != '1' & ksads_13_930_p != '1' & ksads_14_855_p != '1' & ksads_14_854_p != '1' & ksads_14_853_p != '1' & ksads_13_938_p != '1' & ksads_13_939_p != '1' & ksads_13_940_p != '1' & ksads_2_831_p != '1' & ksads_2_830_p != '1' & ksads_2_832_p != '1' & ksads_2_834_p != '1' & ksads_2_833_p != '1' & ksads_2_836_p != '1' & ksads_2_835_p != '1' & ksads_2_837_p != '1' & ksads_13_936_p != '1' & ksads_13_935_p != '1' & ksads_13_937_p != '1' & ksads_20_880_p != '1' & ksads_20_871_p != '1' & ksads_16_900_p != '1' & ksads_16_899_p != '1' & ksads_16_898_p != '1' & ksads_16_897_p != '1' & ksads_3_848_p != '1' & ksads_10_870_p != '1' & ksads_10_869_p != '1' & ksads_20_887_p != '1' & ksads_20_878_p != '1' & ksads_1_841_p != '1' & ksads_1_842_p != '1' & ksads_1_840_p != '1' & ksads_11_918_p != '1' & ksads_11_917_p != '1' & ksads_20_884_p != '1' & ksads_20_875_p != '1' & ksads_15_902_p != '1' & ksads_15_901_p != '1' & ksads_20_888_p != '1' & ksads_20_879_p != '1' & ksads_20_885_p != '1' & ksads_20_876_p != '1' & ksads_6_908_p != '1' & ksads_10_913_p != '1' & ksads_10_914_p != '1' & ksads_5_906_p != '1' & ksads_25_915_p != '1' & ksads_25_916_p != '1' & ksads_7_909_p != '1' & ksads_7_910_p != '1' & ksads_8_911_p != '1' & ksads_8_912_p != '1' & ksads_5_907_p != '1' & ksads_13_941_p != '1' & ksads_13_944_p != '1' & ksads_13_942_p != '1' & ksads_13_943_p != '1' & ksads_18_903_p != '1' & ksads_11_920_p != '1' & ksads_11_919_p != '1' & ksads_21_924_p != '1' & ksads_21_923_p != '1' & ksads_5_858_p != '1' & ksads_5_857_p != '1' & ksads_1_844_p != '1' & ksads_1_845_p != '1' & ksads_1_843_p != '1' & ksads_20_886_p != '1' & ksads_20_877_p != '1' & ksads_21_922_p != '1' & ksads_21_921_p != '1' & ksads_20_882_p != '1' & ksads_20_873_p != '1' & ksads_25_866_p != '1' & ksads_25_865_p != '1' & ksads_7_862_p != '1' & ksads_7_861_p != '1' & ksads_8_864_p != '1' & ksads_8_863_p != '1' & ksads_9_868_p != '1' & ksads_9_867_p != '1' & ksads_20_883_p != '1' & ksads_20_881_p != '1' & ksads_20_874_p != '1' & ksads_20_872_p != '1' & ksads_20_889_p != '1' & ksads_20_890_p != '1' & ksads_19_896_p != '1' & ksads_19_895_p != '1' & ksads_14_856_p != '1' & ksads_2_838_p != '1' & ksads_2_839_p != '1' & ksads_1_846_p != '1' & ksads_1_847_p != '1' & ksads_4_851_p != '1' & ksads_4_852_p != '1' & ksads_20_894_p != '1' & ksads_20_893_p != '1' & ksads_17_905_p != '1' & ksads_17_904_p != '1'))

#1.3 Create subject list for baseline control data
baseline_control_sub_list <- data.frame(control_data_0_dx_at_baseline$subjectkey)

#1.4 Create excel files for the data above
write_xlsx(control_data_0_dx_at_baseline, "ABCD_control_data_baseline.xlsx")
write_xlsx(baseline_control_sub_list, "ABCD_control_sublist_baseline.xlsx")

#2. Number of assessment time points in KSADS clinical data
unique_timepoints <- unique(abcd_ksad_data_amended_dx[c("eventname")])
print(unique_timepoints)

#3. Control data set for no diagnoses at baseline, 1 year and 2 year followups
#3.1 Read in the data
ksads_raw_data <- read_xlsx("./data_raw/ksads_raw_data.csv")

#4.1 Remove the first column containing variable descriptions 
ksads_diagnoses <- ksads_raw_data[,1]

#4.2 Remove the first row of the clinical to account for the descriptions of the variables in the pertinent diagnoses and raw data
abcd_ksads_no_desc <- abcd_ksad_data_amended_dx[-1,]

#4.31 Create a vector containing the relevant demographic and disorder columns to retain in the raw clinical data
cols_2_keep <- c("subjectkey", "eventname", "interview_age", "sex", "interview_date", ksads_diagnoses$Variable)

#4.32 Create a new df containing only the relevant demographic and disorder columns
abcd_ksads_for_analysis <- abcd_ksads_no_desc[, names(abcd_ksads_no_desc) %in% cols_2_keep]

#4.4 re code the "555" and "888" values in the dx relevant column range to allow keeping only "0" values in later sub-setting steps 
abcd_ksads_for_analysis[ , 6:144 ][ abcd_ksads_for_analysis[ , 6:144 ] == "555" ] <- NA
abcd_ksads_for_analysis[ , 6:144 ][ abcd_ksads_for_analysis[ , 6:144 ] == "888" ] <- 0

#4.51 Create a df containing only the demographic information/variables
abcd_ksads_demographic <- abcd_ksads_for_analysis[ , 1:5 ]

#4.52 Create a df containing only the disorder information/variables
abcd_ksads_diagnoses_data <- abcd_ksads_for_analysis[ , 6:144 ]

#4.6 Alter the values in the disorder information/variables df columns to numeric from character
abcd_ksads_numeric_diagnoses <- apply(abcd_ksads_diagnoses_data, 2, as.numeric) %>% as.data.frame()

#4.7 Change the row names in the disorder information/variables df to unique subject ID's containing event label
rownames(abcd_ksads_numeric_diagnoses) <- paste(abcd_ksads_demographic$subjectkey, abcd_ksads_demographic$eventname, sep = "_")

#4.8 Create a df containing the subjects/subject data who received a positive ("1") diagnosis
abcd_ksads_positive_dx <- abcd_ksads_numeric_diagnoses[rowSums(abcd_ksads_numeric_diagnoses>=1, na.rm = TRUE),] %>% as.data.frame()

#4.91 Create a replicated positive dx df to manipulate before comparing against original list
abcd_ksads_positive_dx_df <- abcd_ksads_positive_dx

#4.92 Create a column in the replicated dx df containing full length subject ID's (event label , etc.)
abcd_ksads_positive_dx_df$subject_ID <- rownames(abcd_ksads_positive_dx_df)

#3.10 Remove the string after the Subject ID so that they can be controlled for unique values positive for a dx at any time point
abcd_ksads_positive_dx_df$subject_ID <- str_extract(abcd_ksads_positive_dx_df$subject_ID, "[^_]*_[^_]*")

#4.11 Create a new column in the numeric dx df adding up any positive dx by/for each subject row
abcd_ksads_numeric_diagnoses$positive_dx_count <- rowSums(abcd_ksads_numeric_diagnoses>=1, na.rm = T)

#4.12 Create a subject ID column in the numeric dx df 
abcd_ksads_numeric_diagnoses$subject_ID <- abcd_ksads_no_desc$subjectkey

#4.13 Create a df of subjects who were positive for a dx at any timepoint
abcd_ksads_patients_any_timepoint <- abcd_ksads_numeric_diagnoses[abcd_ksads_numeric_diagnoses$positive_dx_count>0 ,]

#5.1 Create a df of subjects who were controls at any timepoint by sub-setting the raw data against the positive patients df 
abcd_ksads_controls_all_timepoints <- abcd_ksads_no_desc[!abcd_ksads_no_desc$subjectkey %in% abcd_ksads_patients_any_timepoint$subject_ID, ]

#5.2 Mutate the controls df to include only unique subject ID's (deleting repeat ID's to get true control list)
abcd_control_subjects <-
  abcd_ksads_controls_all_timepoints$subjectkey[!duplicated(abcd_ksads_controls_all_timepoints$subjectkey)] %>% as.data.frame()

#5.3 Rename the column in the unique control ID df to a subject ID appropriate name (i.e., for merging df with smri data)
abcd_control_subjects <- rename(abcd_control_subjects, participant = .)


## Output ##

#1. Write excel files for the control (all timepoints) data
write_xlsx(abcd_control_subjects, "./data_processed/main_analysis/ABCD_control_data_all_timepoints.xlsx")
