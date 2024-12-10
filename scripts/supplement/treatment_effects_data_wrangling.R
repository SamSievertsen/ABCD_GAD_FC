## Set Up ##

# Load in necessary packages and configure environmental variables
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)
options(digits = 8, scipen = 999) 

# Read in the treatment data
mental_health_treatment_data <- read.csv("./data_raw/mh_p_ksads_bg.csv")

# Read in the medication data
medication_data <- read_xlsx("./data_raw/ph_p_meds.xlsx")


## Data Wrangling ##

#1. Clean the mental health treatment data and create variables for analysis
#1.1 Retain only columns of interest
mental_health_treatment_data_filtered <- mental_health_treatment_data %>% dplyr::select(c(src_subject_id, eventname, kbi_p_c_age_services, kbi_p_c_scheck1, kbi_p_c_scheck2, kbi_p_c_scheck3, kbi_p_c_scheck7, kbi_p_c_scheck8, kbipcserviceschecklistl1, kbipcserviceschecklistl2, kbipcserviceschecklistl3, kbipcserviceschecklistl7, kbipcserviceschecklistl8))

#1.2 Retain only time points of interest
mental_health_treatment_data_filtered <- mental_health_treatment_data_filtered %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "1_year_follow_up_y_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#1.3 Create the binary treatment variables of interest
#1.31 Create binary variables for treatment types at each time point
mental_health_treatment_data_filtered <- mental_health_treatment_data_filtered %>%
  mutate(
    
    #1.311 Baseline treatment variables
    outpatient_baseline = if_else(kbi_p_c_scheck1 == 1, 1, 0),
    partial_hospitalization_baseline = if_else(kbi_p_c_scheck2 == 1, 1, 0),
    inpatient_baseline = if_else(kbi_p_c_scheck3 == 1, 1, 0),
    psychotherapy_baseline = if_else(kbi_p_c_scheck7 == 1, 1, 0),
    medication_management_baseline = if_else(kbi_p_c_scheck8 == 1, 1, 0),
    
    #1.312 Follow-up treatment variables
    outpatient_followup = if_else(
      eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") &
        kbipcserviceschecklistl1 == 1, 1, 0),
    partial_hospitalization_followup = if_else(
      eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") &
        kbipcserviceschecklistl2 == 1, 1, 0),
    inpatient_followup = if_else(
      eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") &
        kbipcserviceschecklistl3 == 1, 1, 0),
    psychotherapy_followup = if_else(
      eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") &
        kbipcserviceschecklistl7 == 1, 1, 0),
    medication_management_followup = if_else(
      eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") &
        kbipcserviceschecklistl8 == 1, 1, 0))

#1.32 Create non-timepoint-specific binary variables for treatment types
mental_health_treatment_data_filtered <-
  mental_health_treatment_data_filtered %>%
  mutate(
    
    #1.321 Collapsing timepoint-specific variables into overall binary variables
    outpatient = pmax(outpatient_baseline, outpatient_followup, na.rm = TRUE),
    partial_hospitalization = pmax(partial_hospitalization_baseline, partial_hospitalization_followup, na.rm = TRUE),
    inpatient = pmax(inpatient_baseline, inpatient_followup, na.rm = TRUE),
    psychotherapy = pmax(psychotherapy_baseline, psychotherapy_followup, na.rm = TRUE),
    medication_management = pmax(medication_management_baseline, medication_management_followup, na.rm = TRUE),
    
    #1.322 Create a variable denoting whether subjects received any type of mental health treatment
    any_mental_health_tx = if_else(outpatient == 1 | partial_hospitalization == 1 | inpatient == 1 | psychotherapy == 1 | medication_management == 1, 1, 0))

#1.4 Retain columns of interest for merging with the medication data + creating additional variables
mental_health_treatment_data_for_merging <- mental_health_treatment_data_filtered %>% 
  dplyr::select(c(src_subject_id, eventname, outpatient, partial_hospitalization, inpatient, psychotherapy, medication_management, any_mental_health_tx))


#2. Clean the medication data
#2.1 Retain only the columns of interest
medication_data_filtered <- medication_data %>% dplyr::select(c(src_subject_id, eventname, med1_rxnorm_p, med2_rxnorm_p, med3_rxnorm_p, med4_rxnorm_p, med5_rxnorm_p, med6_rxnorm_p, med7_rxnorm_p))

#2.2 Retain only time points of interest
medication_data_filtered <- medication_data_filtered %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "1_year_follow_up_y_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.3 Create a variable for each subject at each timepoint denoting whether they are taking an anti-depressant or anxiolytic medication
#2.31 Define the lists of anit-depressant and anxiolytic medication strings
#2.311 Anti-depressants
antidepressants <- c("escitalopram", "fluoxetine", "sertraline", "mirtazapine", 
                     "citalopram", "fluvoxamine", "imipramine", "bupropion", 
                     "amitriptyline", "trazodone")

#2.312 Anxiolytics
anxiolytics <- c("hydroxyzine", "buspirone", "diazepam", "alprazolam")

#2.32 Check whether each subject at each timepoint was taking the medication(s) of interest and create related variables
medication_data_filtered <- medication_data_filtered %>%
  mutate(
    
    # 2.321 Check for any anti-depressants taken in the medication columns
    taking_antidepressant = as.numeric(if_any(starts_with("med"), 
                                              ~ grepl(paste(antidepressants, collapse = "|"), 
                                                      ., ignore.case = TRUE))),
    
    # 2.322 Check for any anxiolytics taken in the medication columns
    taking_anxiolytic = as.numeric(if_any(starts_with("med"), 
                                          ~ grepl(paste(anxiolytics, collapse = "|"), 
                                                  ., ignore.case = TRUE))),
    
    # 2.323 Denote whether subjects were taking either anti-depressants or anxiolytics
    taking_medication = if_else(taking_antidepressant == 1 | taking_anxiolytic == 1, 1, 0))

#2.4 Retain columns of interest for merging with the treatment data + creating additional variables
medication_data_for_merging <- medication_data_filtered %>% 
  dplyr::select(c(src_subject_id, eventname, taking_antidepressant, taking_anxiolytic, taking_medication))


#3. Merge the mental health treatment and medication data + create additional variables for analysis
#3.1 Join the mental health + medication treatment data
treatment_data_raw <- full_join(mental_health_treatment_data_for_merging, medication_data_for_merging)

#3.2 Create a variable representing whether subjects at each timepoint received any type of treatment 
treatment_data_raw <- treatment_data_raw %>% 
  mutate(received_any_treatment = if_else(outpatient == 1 | partial_hospitalization == 1 | inpatient == 1 | psychotherapy == 1 | medication_management == 1 | taking_medication == 1, 1, 0))

#3.3 Create a variable denoting when subjects received treatment
treatment_data_raw <- treatment_data_raw %>%
  group_by(src_subject_id) %>%
  mutate(
    
    #3.31 Check if the subject has either baseline + 2-year follow-up or baseline + 1-year follow-up data
    has_baseline_and_followup = any(eventname == "baseline_year_1_arm_1") &
      any(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1")),
    
    #3.32 Check if treatment was received by baseline
    baseline_treatment = any(eventname == "baseline_year_1_arm_1" & received_any_treatment == 1, na.rm = TRUE),
    
    #3.33 Check if treatment was initiated between baseline and follow-up (1 year or 2 year)
    followup_treatment = any(eventname %in% c("1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1") & 
                               received_any_treatment == 1, na.rm = TRUE),
    
    #3.34 Assign treatment status based on conditions, only for subjects with baseline and follow-up data
    treatment_status = if_else(
      has_baseline_and_followup,
      case_when(
        baseline_treatment & followup_treatment ~ "baseline_and_followup_tx",
        baseline_treatment & !followup_treatment ~ "tx_before_baseline",
        !baseline_treatment & followup_treatment ~ "tx_after_baseline",
        !baseline_treatment & !followup_treatment ~ "no_tx_history",
        TRUE ~ NA_character_), NA_character_ # Assign NA if the subject doesn't have required time points
    )
  ) %>%
  ungroup()


#3.4 Retain columns of interest at timepoints of interest for analysis purposes 
treatment_data <- treatment_data_raw %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1") %>% 
  dplyr::select(-c(baseline_treatment, followup_treatment, has_baseline_and_followup))


## Output ## 

#1. Write the treatment data as a csv file 
write.csv(treatment_data, "./data_processed/supplement/treatment_data.csv", row.names = FALSE)
