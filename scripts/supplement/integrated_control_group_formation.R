## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
options(digits = 8, scipen = 999) 

# Read in required data 
# ABCD release 5.1 parent report KSADS-COMP data
ABCD_KSADS_release_5.1_parent_report <- read.csv("data_raw/mh_p_ksads_ss.csv")

# ABCD release 5.1 youth report KSADS-COMP data
ABCD_KSADS_release_5.1_youth_report <- read.csv("data_raw/mh_y_ksads_ss.csv")

# KSADS-COMP variable names of interest to the healthy control group formation
KSADS_variable_names <- read.csv("data_raw/ksads_variable_desc.csv")

## Data Wrangling ## 

#1. Merge the parent and youth report data
#1.1 Full join the parent and youth report KSADS dataframes
ABCD_KSADS_release_5.1_data <- full_join(ABCD_KSADS_release_5.1_parent_report, ABCD_KSADS_release_5.1_youth_report)

#2 Clean the parent and youth report data
#2.1 Retain only the columns of interest to creating the healthy control group 
ABCD_KSADS_release_5.1_data <- ABCD_KSADS_release_5.1_data %>% 
  dplyr::select(c(src_subject_id, eventname, all_of(KSADS_variable_names$variable_name)))

#2.2 Retain only the baseline and 2 year follow up data
ABCD_KSADS_release_5.1_data <- ABCD_KSADS_release_5.1_data %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.3 Convert any instances of NA, empty string, 555, or 999 to NA; and any instances of 888 to 0
ABCD_KSADS_release_5.1_data <- ABCD_KSADS_release_5.1_data %>%
  mutate_all(~ ifelse(. %in% c(999, 555, ""), NA, ifelse(. == 888, 0, .)))

#2.4 Ensure all relevant columns are numeric type
ABCD_KSADS_release_5.1_data <- ABCD_KSADS_release_5.1_data %>%
  mutate_at(vars(3:220), as.numeric)


## Healthy Control Group Creation ##

#1. Identify all diagnosis-related columns, which start with "ksads_"
KSADS_diagnosis_columns <- grep("^ksads_", colnames(ABCD_KSADS_release_5.1_data), value = TRUE)

#2. Create a new variable indicating if a subject has any diagnoses across all rows
ABCD_KSADS_release_5.1_data <- ABCD_KSADS_release_5.1_data %>%
  group_by(src_subject_id) %>%
  mutate(has_diagnosis = any(across(all_of(KSADS_diagnosis_columns), ~ .x == 1))) %>%
  mutate(any_non_na = rowSums(!is.na(across(all_of(KSADS_diagnosis_columns)))) > 0) %>%
  ungroup()

#3. Store subjects who have no diagnoses and at least one non-NA value as healthy controls
healthy_control_group <- ABCD_KSADS_release_5.1_data %>%
  filter(is.na(has_diagnosis) & any_non_na) %>%
  pull(src_subject_id) %>% 
  unique()

## Output ##

#1. Create a dataframe containing the subject IDs of healthy controls
#1.1 Convert the current vector to a datframe 
healthy_control_group <- data.frame(src_subject_id = healthy_control_group)

#1.2 Add the group and subgroup variables
healthy_control_group <- healthy_control_group %>% 
  mutate(group = rep("HC"),
         subgroup = rep("HC"))

#2. Write the dataframe as a csv 
write.csv(healthy_control_group, "./data_processed/healty_control_group.csv", row.names = FALSE)


