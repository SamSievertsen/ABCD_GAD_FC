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
ABCD_rsfMRI_Data <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/ABCD_rsfMRI_Data.csv")


## Data Wrangling ## 

# Remove the first row of the df to account for the descriptions of the variables in the pertinent diagnoses and raw data
ABCD_rsfMRI_Data_no_desc <- ABCD_rsfMRI_Data[-1,]
# Create the site_name variable containing the scanner and site info
ABCD_rsfMRI_Data_no_desc$site_name <- substr(ABCD_rsfMRI_Data_no_desc$rsfmri_c_ngd_visitid, 1, 4)
# Read in the Family ID data 
abcd_family_id_data <- read.delim("C:/Users/Sam Sievertsen/Desktop/ABCD/ABCC_Package_1203705_Tabulated-Data/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, rel_family_id))
# Change the name of the family ID column
abcd_family_id_data <- abcd_family_id_data %>% rename(family_id = rel_family_id)
# Merge the family ID data with the resting state data
ABCD_rsfMRI_Data_no_desc <- left_join(ABCD_rsfMRI_Data_no_desc, abcd_family_id_data)
# Read in analysis four (analysis 2 group connectivity difference) subjectkey's and merge with resting state data
analysis_five_GAD_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_GAD_subjectkeys_groups.csv")
analysis_five_CN_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_CN_subjectkeys_groups.csv")
analysis_five_sample_groups <- full_join(analysis_five_CN_sub_groups, analysis_five_GAD_sub_groups)
# Join the grouped subjectkeys with the rsfMRI imaging data. Keep only the vars of interest 
analysis_five_grouped_imaging_data <- merge(analysis_five_sample_groups, ABCD_rsfMRI_Data_no_desc)
analysis_five_grouped_imaging_data$age_in_years <- floor(as.numeric(as.character(analysis_five_grouped_imaging_data$interview_age)) / 12)
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data %>% 
  dplyr::select(c(subjectkey, eventname, group, interview_age, age_in_years, sex, site_name, family_id, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))
# Change the values of the "group" variable to represent their experimental identity
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "control"] <- "Control"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "baseline_GAD"] <- "GAD Remitter"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "followup_GAD"] <- "GAD Converter"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "GAD_Both"] <- "Continuous GAD"
# make numerical values numeric
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data %>%
  mutate_at(vars(9:14), as.numeric)
# make factor values factor type
analysis_five_grouped_imaging_data$eventname <- as.factor(analysis_five_grouped_imaging_data$eventname)
analysis_five_grouped_imaging_data$group <- as.factor(analysis_five_grouped_imaging_data$group)
analysis_five_grouped_imaging_data$sex <- as.factor(analysis_five_grouped_imaging_data$sex)
analysis_five_grouped_imaging_data$site_name <- as.factor(analysis_five_grouped_imaging_data$site_name)
analysis_five_grouped_imaging_data$subjectkey <- as.factor(analysis_five_grouped_imaging_data$subjectkey)
analysis_five_grouped_imaging_data$family_id <- as.factor(analysis_five_grouped_imaging_data$family_id)
# Remove any rows with NA or empty values and retain subjects who have both baseline and followup scans only (to perform subtraction)
analysis_five_grouped_imaging_data <- na.omit(analysis_five_grouped_imaging_data)
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data[rowSums(analysis_five_grouped_imaging_data == "") == 0, ]
analysis_five_grouped_imaging_data <-  analysis_five_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") &
           any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()
# Set the reference level of the diagnostic group to "controls" and the eventname to "baseline" 
analysis_five_grouped_imaging_data$group <- relevel(analysis_five_grouped_imaging_data$group, ref = "Control")
analysis_five_grouped_imaging_data$eventname <- relevel(analysis_five_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")