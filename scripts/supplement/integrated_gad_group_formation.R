## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(tidyr)
options(digits = 8, scipen = 999) 

# Read in required data 
# ABCD release 5.1 parent report KSADS-COMP data
ABCD_KSADS_release_5.1_parent_report <- read.csv("data_raw/mh_p_ksads_ss.csv")

# ABCD release 5.1 youth report KSADS-COMP data
ABCD_KSADS_release_5.1_youth_report <- read.csv("data_raw/mh_y_ksads_ss.csv")


## Data Wrangling ##

#1. Merge the parent and youth report data
#1.1 Full join the parent and youth report KSADS dataframes
ABCD_KSADS_release_5.1_data <- full_join(ABCD_KSADS_release_5.1_parent_report, ABCD_KSADS_release_5.1_youth_report)

#2 Clean the parent and youth report data
#2.1 Retain only the columns of interest to creating the GAD group(s)
ABCD_KSADS_release_5.1_GAD_data <- ABCD_KSADS_release_5.1_data %>% 
  dplyr::select(c(src_subject_id, eventname, ksads_10_869_p, ksads_10_869_t))

#2.2 Retain only the baseline and 2 year follow up data
ABCD_KSADS_release_5.1_GAD_data <- ABCD_KSADS_release_5.1_GAD_data %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.3 Convert any instances of NA, empty string, 555, or 999 to NA; and any instances of 888 to 0
ABCD_KSADS_release_5.1_GAD_data <- ABCD_KSADS_release_5.1_GAD_data %>%
  mutate_all(~ ifelse(. %in% c(999, 555, ""), NA, ifelse(. == 888, 0, .)))


## GAD Group(s) Formation ##

#1. Determine subjects with current GAD at either timepoint
#1.1 Create the primary current GAD group
GAD_group <- ABCD_KSADS_release_5.1_GAD_data %>% 
  mutate(group = if_else(ksads_10_869_p == 1 | ksads_10_869_t == 1, "GAD", NA))

#1.2 Create the first analysis (current GAD vs HC) groups
#1.21 Determine at which time point (baseline, follow up, both) participants had GAD
GAD_group <- GAD_group %>%
  group_by(src_subject_id) %>%
  mutate(GAD_timepoint = case_when(
      any(group[eventname == "baseline_year_1_arm_1"] == "GAD", na.rm = TRUE) &
        any(group[eventname == "2_year_follow_up_y_arm_1"] == "GAD", na.rm = TRUE) ~ "both",
      any(group[eventname == "baseline_year_1_arm_1"] == "GAD", na.rm = TRUE) ~ "baseline",
      any(group[eventname == "2_year_follow_up_y_arm_1"] == "GAD", na.rm = TRUE) ~ "followup",
      TRUE ~ NA_character_)) %>%
  ungroup()

#1.22 Craft a variable to be used in the main GAD vs HC group analysis wherein subjects with GAD at both timepoints are randomly assigned to either the baseline or follow up group
GAD_group <- GAD_group %>%
  group_by(src_subject_id) %>%
  mutate(random_assignment = if_else(
    all(GAD_timepoint == "both"),
    sample(c("baseline", "followup"), size = 1, replace = TRUE), 
    NA_character_)) %>%
  ungroup() %>% 
  mutate(analysis_group = case_when(
    GAD_timepoint %in% c("baseline", "followup") ~ GAD_timepoint,
    GAD_timepoint == "both" ~ random_assignment,
    TRUE ~ NA_character_))
#1.221 Remove the random_assignment variable
GAD_group <- GAD_group %>% dplyr::select(-random_assignment)

#1.3 Create the GAD subgroups for the repeated measures analysis
GAD_group <- GAD_group %>%
  group_by(src_subject_id) %>%
  mutate(subgroup = case_when(
    any(eventname == "baseline_year_1_arm_1" & group == "GAD") & 
      any(eventname == "2_year_follow_up_y_arm_1" & group == "GAD") ~ "Continuous GAD",
    eventname == "baseline_year_1_arm_1" & group == "GAD" & 
      any(eventname == "2_year_follow_up_y_arm_1" & is.na(group)) ~ "GAD Remitter",
    eventname == "baseline_year_1_arm_1" & is.na(group) & 
      any(eventname == "2_year_follow_up_y_arm_1" & group == "GAD") ~ "GAD Converter",
    TRUE ~ NA_character_)) %>%
  tidyr::fill(subgroup) %>%
  ungroup()


## Output ##

#1. Clean the GAD group data
#1.1 Retain only GAD sample subjects
GAD_group_data <- GAD_group %>% 
  filter(group == "GAD")

#1.2 Retain only columns of interest 
GAD_group_data <- GAD_group_data %>% 
  dplyr::select(c(src_subject_id, eventname, group, GAD_timepoint, analysis_group, subgroup))

#2. Write the GAD group data as a csv
write.csv(GAD_group_data, "./data_processed/supplement/GAD_group.csv", row.names = FALSE)
