library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(stringr)

mh_tx_test_5.1 <- read.csv("./data_raw/mh_p_ksads_bg.csv")

mh_tx_test_5.1_filtered <- mh_tx_test_5.1 %>% dplyr::select(c(src_subject_id, eventname, kbi_p_c_mh_sa, kbi_p_c_age_services, kbi_p_c_scheck1, kbi_p_c_scheck2, kbi_p_c_scheck3, kbi_p_c_scheck7, kbi_p_c_scheck8, kbi_p_c_mental_health, kbi_ss_c_mental_health_p, kbi_p_c_mh_sa_l, kbipcserviceschecklistl1, kbipcserviceschecklistl2, kbipcserviceschecklistl3, kbipcserviceschecklistl7, kbipcserviceschecklistl8, kbi_p_c_mental_health_l, kbi_ss_c_mental_health_p_l, kbi_p_c_serviceschecklist3l8, kbi_p_c_serviceschecklist3l1, kbi_p_c_serviceschecklist3l2, kbi_p_c_serviceschecklist3l3, kbi_p_c_serviceschecklist3l7))

meds_test_5.1 <- read_xlsx("./data_raw/ph_p_meds.xlsx")

meds_test_5.1_filtered <- meds_test_5.1 %>% dplyr::select(c(src_subject_id, eventname, med1_rxnorm_p, med2_rxnorm_p, med3_rxnorm_p, med4_rxnorm_p, med5_rxnorm_p, med6_rxnorm_p, med7_rxnorm_p))

main_data <- read.csv("/Users/samsievertsen/Desktop/SamResearch/ABCD_Data/site_visit_analysis_data.csv")

main_data_filtered <- main_data %>% dplyr::select(c(subjectkey, eventname, group, analysis_group))


######
main_data_filtered_w_meds <- left_join(main_data_filtered, meds_test_5.1_filtered, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))

# Define the lists of medication strings
antidepressants <- c("escitalopram", "fluoxetine", "sertraline", "mirtazapine", 
                     "citalopram", "fluvoxamine", "imipramine", "bupropion", 
                     "amitriptyline", "trazodone")
anxiolytics <- c("hydroxyzine", "buspirone", "diazepam", "alprazolam")

# Check medications for each subject and create flags
meds_table <- main_data_filtered_w_meds %>%
  mutate(
    # Check for any anti-depressant in the medication columns
    taking_antidepressant = if_any(starts_with("med"), 
                                   ~ grepl(paste(antidepressants, collapse = "|"), 
                                           ., ignore.case = TRUE)),
    # Check for any anxiolytic in the medication columns
    taking_anxiolytic = if_any(starts_with("med"), 
                               ~ grepl(paste(anxiolytics, collapse = "|"), 
                                       ., ignore.case = TRUE))
  ) %>%
  # Group by analysis_group and eventname
  group_by(analysis_group, eventname) %>%
  # Summarize counts of subjects taking each medication type
  summarise(
    subjects_on_antidepressants = sum(taking_antidepressant, na.rm = TRUE),
    subjects_on_anxiolytics = sum(taking_anxiolytic, na.rm = TRUE),
    total_subjects = n_distinct(subjectkey),
    .groups = "drop"
  )

# View the result
print(meds_table)



##############################
main_data_filtered_w_tx <- left_join(main_data_filtered, mh_tx_test_5.1_filtered, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))

tabled_mh_tx_values <- main_data_filtered_w_tx %>%
  # Pivot columns 5:26 into a long format
  pivot_longer(cols = 5:26, names_to = "variable", values_to = "value") %>%
  # Group by eventname, group, and variable
  group_by(eventname, analysis_group, variable, value) %>%
  # Count occurrences of each value
  summarise(count = n(), .groups = "drop") %>%
  # Optionally pivot back to wide format (if desired)
  pivot_wider(names_from = value, values_from = count, values_fill = 0)

