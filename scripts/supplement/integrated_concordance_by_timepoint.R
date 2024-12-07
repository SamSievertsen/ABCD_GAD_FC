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

# ABCD release 4.0 CBCL data
ABCD_CBCL_data <- read.csv("./data_raw/abcd_cbcls01.csv")
ABCD_CBCL_data <- ABCD_CBCL_data[-1,]

# ABCD release 4.0 BPM data
ABCD_BPM_data <- read.csv("./data_raw/abcd_yssbpm01.csv")
ABCD_BPM_data <- ABCD_BPM_data[-1,]

# Cleaned and QC'd imaging data
cleaned_qcd_imaging_data <- read.csv("./data_processed/rsfMRI_data_qcd_unsubset.csv")

# Diagnosis group labels
GAD_HC_group <- read.csv("./data_processed/GAD_HC_resampled_merged_groups.csv") %>% 
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
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_filtered %>%
  mutate_all(~ ifelse(. %in% c(999, 555, ""), NA, ifelse(. == 888, 0, .)))

#1.3 Bring in HC + GAD diagnosis labels
#1.31 Join the HC + GAD labels into the KSADS data
ABCD_KSADS_release_5.1_GAD_data_for_merging <- left_join(ABCD_KSADS_release_5.1_GAD_data_for_merging, GAD_HC_group)

#1.32 Create a new GAD diangostic label determining the source of the GAD diagnosis for this set of descriptive stats 
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
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
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
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
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
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
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
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
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>% 
  dplyr::select(c(src_subject_id, eventname, group, GAD_Source, Social_Anxiety_Disorder, Social_Anxiety_Disorder_Source, Separation_Anxiety_Disorder, Separation_Anxiety_Disorder_Source, MDD, MDD_Source))

#1.6 Clean subjects in the KSADS data to merge
#1.61 Remove subjects whos rsfMRI data is incomplete or did not pass QA/QC
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
  semi_join(cleaned_qcd_imaging_data, by = c("src_subject_id", "eventname"))

#1.62 Remove subjects that get tagged as GAD but aren't included in the final sample
ABCD_KSADS_release_5.1_GAD_data_for_merging <- ABCD_KSADS_release_5.1_GAD_data_for_merging %>%
filter(!(is.na(group) & GAD_Source != "None"))


#2. Prep the CBCL/BPM data
#2.1 CBCL
#2.11 Retain only CBCL columns of interest
CBCL_data_filtered <- ABCD_CBCL_data %>%
  dplyr::select(c(src_subject_id, eventname, sex, cbcl_scr_syn_anxdep_t, cbcl_scr_syn_anxdep_nm, cbcl_scr_syn_internal_t, cbcl_scr_syn_internal_nm, cbcl_scr_syn_external_t, cbcl_scr_syn_external_nm))

#2.12 Subset the baseline and 2 year follow up data
CBCL_data_filtered <- CBCL_data_filtered %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.13 Keep rows where the CBCL scores are not NA or empty, and where the number of questions missing that went into the calculation of each CBCL score is 0
CBCL_data_filtered_complete <- subset(CBCL_data_filtered, !is.na(cbcl_scr_syn_internal_t) & cbcl_scr_syn_internal_t != "" & !is.na(cbcl_scr_syn_external_t) & cbcl_scr_syn_external_t != "" & !is.na(cbcl_scr_syn_anxdep_t) & cbcl_scr_syn_anxdep_t != "" & cbcl_scr_syn_internal_nm == 0 & cbcl_scr_syn_external_nm == 0 & cbcl_scr_syn_anxdep_nm == 0)

#2.14 Subset only columns of interest
CBCL_data_for_merging <- CBCL_data_filtered_complete %>%
  dplyr::select(c(src_subject_id, eventname, sex, cbcl_scr_syn_anxdep_t, cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t))


#2.2 BPM
#2.21 Retain only BPM columns of interest
BPM_data_filtered <- ABCD_BPM_data %>% 
  dplyr::select(c(src_subject_id, eventname, bpm_y_scr_internal_t, bpm_y_scr_internal_nm, bpm_y_scr_external_t, bpm_y_scr_external_nm))

#2.22 Subset the baseline and 2 year follow up data
BPM_data_filtered <- BPM_data_filtered %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.23 Keep rows where the BPM scores are not NA or empty, and where the number of questions missing that went into the calculation of each BPM score is 0
BPM_data_filtered_complete <- subset(BPM_data_filtered, !is.na(bpm_y_scr_internal_t) & bpm_y_scr_internal_t != "" & !is.na(bpm_y_scr_external_t) & bpm_y_scr_external_t != "" & bpm_y_scr_internal_nm == 0 & bpm_y_scr_external_nm == 0)

#2.23 Again retain only the columns of interest
BPM_data_for_merging <- BPM_data_filtered_complete %>% 
  dplyr::select(c(src_subject_id, eventname, bpm_y_scr_internal_t, bpm_y_scr_external_t))

#2.3 Merge the CBCL + BPM data
CBCL_BPM_data_for_merging <- full_join(CBCL_data_for_merging, BPM_data_for_merging)

#3. Merge + clean the KSADS, CBCL & BPM data
#3.1 Join the KSADS, CBCL & BPM data
KSADS_CBCL_BPM_merged_data <- full_join(ABCD_KSADS_release_5.1_GAD_data_for_merging, CBCL_BPM_data_for_merging)

#3.2 Ensure all relevant columns are numeric type
KSADS_CBCL_BPM_merged_data <- KSADS_CBCL_BPM_merged_data %>%
  mutate_at(vars(12:16), as.numeric)


## Data Analysis ##

#1. Check whether these self-reported-only GAD youth look different than parent report only GAD on the parent-report + youth-report CBCL and BPM
#1.11 Summarize the mean CBCL and BPM scores for each group
parent_vs_youth_report_GAD_cbcl_comparison <- KSADS_CBCL_BPM_merged_data %>%
  filter(!is.na(GAD_Source)) %>%
  group_by(GAD_Source, eventname) %>%
  summarise(
    Group_N = n(), 
    CBCL_anxdep_mean = round(mean(cbcl_scr_syn_anxdep_t, na.rm = TRUE), digits = 2),
    CBCL_anxdep_sd = round(sd(cbcl_scr_syn_anxdep_t, na.rm = TRUE), digits = 2),
    CBCL_internalizing_mean = round(mean(cbcl_scr_syn_internal_t, na.rm = TRUE), digits = 2),
    CBCL_internalizing_sd = round(sd(cbcl_scr_syn_internal_t, na.rm = TRUE), digits = 2), 
    CBCL_externalizing_mean = round(mean(cbcl_scr_syn_external_t, na.rm = TRUE), digits = 2),
    CBCL_externalizing_sd = round(sd(cbcl_scr_syn_external_t, na.rm = TRUE), digits = 2),
    BPM_internalizing_mean = round(mean(bpm_y_scr_internal_t, na.rm = TRUE), digits = 2),
    BPM_internalizing_sd = round(sd(bpm_y_scr_internal_t, na.rm = TRUE), digits = 2), 
    BPM_externalizing_mean = round(mean(bpm_y_scr_external_t, na.rm = TRUE), digits = 2),
    BPM_externalizing_sd = round(sd(bpm_y_scr_external_t, na.rm = TRUE), digits = 2))

#1.121 Reshape data for easier plotting by selecting relevant columns and metrics
parent_vs_youth_GAD_CBCL_long <- KSADS_CBCL_BPM_merged_data %>%
  filter(!is.na(GAD_Source)) %>%
  dplyr::select(c(GAD_Source, eventname, cbcl_scr_syn_anxdep_t, cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t, bpm_y_scr_internal_t, bpm_y_scr_external_t)) %>%
  pivot_longer(
    cols = starts_with("cbcl") | starts_with("bpm"),
    names_to = "Metric",
    values_to = "Score") %>%
  filter(!is.na(Score))

#1.122 Recode the eventname values
parent_vs_youth_GAD_CBCL_long <- parent_vs_youth_GAD_CBCL_long %>%
  mutate(eventname = recode(eventname, "baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 yr Follow Up"))

#1.123 Rename the `Metric` column values
parent_vs_youth_GAD_CBCL_long <- parent_vs_youth_GAD_CBCL_long %>%
  mutate(Metric = recode(Metric,
                         "cbcl_scr_syn_anxdep_t" = "CBCL AnxDep T",
                         "cbcl_scr_syn_internal_t" = "CBCL Internalizing T",
                         "cbcl_scr_syn_external_t" = "CBCL Externalizing T",
                         "bpm_y_scr_internal_t" = "BPM Internalizing T",
                         "bpm_y_scr_external_t" = "BPM Externalizing T"))

#1.124 Plot the CBCL and BPM data by group
ggplot(parent_vs_youth_GAD_CBCL_long, aes(x = eventname, y = Score, fill = GAD_Source)) + 
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.3, width = 0.7) + 
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA) + 
  facet_wrap(~ Metric, scales = "free_y", nrow = 3) + 
  labs(title = "Parent vs. Youth Reported GAD Scores Across Metrics and Events", x = "Event Name", y = "Score", fill = "Source") + 
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top", strip.text = element_text(size = 12), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))

#1.2 Perform statistical testing to check whether the groups differ from each other or not in 
parent_vs_youth_report_GAD_cbcl_ANCOVA <- list(
  anova_cbcl_anxdep = car::Anova(
    lm(cbcl_scr_syn_anxdep_t ~ GAD_Source, 
       data = KSADS_CBCL_BPM_merged_data, 
       subset = GAD_Source %in% c("Youth", "Parent")),
    type = "II"),
  anova_cbcl_internal = car::Anova(
    lm(cbcl_scr_syn_internal_t ~ GAD_Source, 
       data = KSADS_CBCL_BPM_merged_data, 
       subset = GAD_Source %in% c("Youth", "Parent")),
    type = "II"),
  anova_cbcl_external = car::Anova(
    lm(cbcl_scr_syn_external_t ~ GAD_Source, 
       data = KSADS_CBCL_BPM_merged_data, 
       subset = GAD_Source %in% c("Youth", "Parent")),
    type = "II"),
  anova_bpm_internal = car::Anova(
    lm(bpm_y_scr_internal_t ~ GAD_Source, 
       data = KSADS_CBCL_BPM_merged_data, 
       subset = GAD_Source %in% c("Youth", "Parent")),
    type = "II"),
  anova_bpm_external = car::Anova(
    lm(bpm_y_scr_external_t ~ GAD_Source, 
       data = KSADS_CBCL_BPM_merged_data, 
       subset = GAD_Source %in% c("Youth", "Parent")),
    type = "II")
)


#2. Check whether the self-reported-only GAD youth have a different KSADS comorbidity profile (social anxiety disorder, separation anxiety disorder, MDD) than the parent reported GAD youth
#2.1 Create a whole sample summary dataframe containing N and % of GAD and comorbid disorders of interest
Whole_Sample_Dx_Summary <- KSADS_CBCL_BPM_merged_data %>%
  filter(group == "GAD") %>%
  group_by(GAD_Source, eventname) %>%  
  summarize(
    
    #2.11 Total number of subjects
    N_Subjects = n(),
    
    #2.12 GAD diagnosis
    N_GAD_Dx = sum(group == "GAD", na.rm = TRUE),
    Percent_GAD_Dx = (N_GAD_Dx / N_Subjects * 100),
    
    #2.131 Comorbidity: Social Anxiety + GAD (using GAD as reference point)
    N_GAD_with_Social_Anxiety_Comorbidity = sum(Social_Anxiety_Disorder == 1 & group == "GAD", na.rm = TRUE),
    Percent_GAD_with_Social_Anxiety_Comorbidity = (N_GAD_with_Social_Anxiety_Comorbidity / N_GAD_Dx * 100),
    
    #2.132 Comorbidity: Separation Anxiety + GAD (using GAD as reference point)
    N_GAD_with_Separation_Anxiety_Comorbidity = sum(Separation_Anxiety_Disorder == 1 & group == "GAD", na.rm = TRUE),
    Percent_GAD_with_Separation_Anxiety_Comorbidity = (N_GAD_with_Separation_Anxiety_Comorbidity / N_GAD_Dx * 100),
    
    #2.133 Comorbidity: MDD + GAD (using GAD as reference point)
    N_GAD_with_MDD_Comorbidity = sum(MDD == 1 & group == "GAD", na.rm = TRUE),
    Percent_GAD_with_MDD_Comorbidity = (N_GAD_with_MDD_Comorbidity / N_GAD_Dx * 100),
    
    #2.134 Any comorbidity (with GAD): At least one non-GAD diagnosis with GAD (unduplicated)
    N_Any_Comorbidity_with_GAD = sum((Social_Anxiety_Disorder == 1 | Separation_Anxiety_Disorder == 1 | MDD == 1) & group == "GAD", na.rm = TRUE),
    Percent_Any_Comorbidity_with_GAD = (N_Any_Comorbidity_with_GAD / N_GAD_Dx * 100),
    
    #2.14 Subjects with only GAD (no comorbidity)
    N_Only_GAD = N_GAD_Dx - N_Any_Comorbidity_with_GAD,
    Percent_Only_GAD = (N_Only_GAD / N_GAD_Dx * 100)) %>%
  
  #2.15 Rounding to 2 decimal places for all numeric columns
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  ungroup() 

# 2.2 Pivot the comorbidity summary dataframe to a format suitable for Notion
Whole_Sample_Long_Dx_Summary <- Whole_Sample_Dx_Summary %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Percent_"), 
    names_to = c(".value", "Diagnosis"), 
    names_pattern = "(N|Percent)_(.*)") %>%
  mutate(Diagnosis = gsub("_", " ", Diagnosis))


#3. Check for youth and parent report convergence
#3.1 Count how many youth and parents both reported each disorder
Divergence_Concordance_Summary <- KSADS_CBCL_BPM_merged_data %>%
  
  #3.11 Filter for subjects in the GAD group only
  filter(group == "GAD") %>%
  
  #3.12 Group data by timepoint (eventname) to calculate divergence and concordance at each timepoint
  group_by(eventname) %>%
  
  #3.13 Summarize the data to count the different reporting patterns for each disorder
  summarize(
    
    #3.131 GAD: Count cases where only the parent reported GAD
    Parent_GAD_Only = sum(GAD_Source == "Parent", na.rm = TRUE),
    
    #3.132 Count total GAD cases to use in percentage calculations
    Total_GAD = sum(group == "GAD", na.rm = TRUE),
    
    #3.133 Calculate percentage of cases where only the parent reported GAD
    Percent_Parent_GAD_Only = (Parent_GAD_Only / Total_GAD * 100),
    
    #3.134 GAD: Count cases where only the youth reported GAD
    Youth_GAD_Only = sum(GAD_Source == "Youth", na.rm = TRUE),
    
    #3.135 Calculate percentage of cases where only the youth reported GAD
    Percent_Youth_GAD_Only = (Youth_GAD_Only / Total_GAD * 100),
    
    #3.136 GAD: Count cases where both parent and youth reported GAD
    Concordant_GAD = sum(GAD_Source == "Concordant", na.rm = TRUE),
    
    #3.137 Calculate percentage of cases where both parent and youth reported GAD
    Percent_Concordant_GAD = (Concordant_GAD / Total_GAD * 100),
    
    #3.141 MDD: Count cases where only the parent reported MDD
    Parent_MDD_Only = sum(MDD_Source == "Parent", na.rm = TRUE),
    
    #3.142 Count total MDD cases to use in percentage calculations
    Total_MDD = sum(MDD == 1, na.rm = TRUE),
    
    #3.143 Calculate percentage of cases where only the parent reported MDD
    Percent_Parent_MDD_Only = (Parent_MDD_Only / Total_MDD * 100),
    
    #3.144 MDD: Count cases where only the youth reported MDD
    Youth_MDD_Only = sum(MDD_Source == "Youth", na.rm = TRUE),
    
    #3.145 Calculate percentage of cases where only the youth reported MDD
    Percent_Youth_MDD_Only = (Youth_MDD_Only / Total_MDD * 100),
    
    #3.146 MDD: Count cases where both parent and youth reported MDD
    Concordant_MDD = sum(MDD_Source == "Concordant", na.rm = TRUE),
    
    #3.146 Calculate percentage of cases where both parent and youth reported MDD
    Percent_Concordant_MDD = (Concordant_MDD / Total_MDD * 100),
    
    #3.151 Social Anxiety Disorder: Count cases where only the parent reported Social Anxiety Disorder
    Parent_Social_Anxiety_Only = sum(Social_Anxiety_Disorder_Source == "Parent", na.rm = TRUE),
    
    #3.152 Count total Social Anxiety Disorder cases
    Total_Social_Anxiety = sum(Social_Anxiety_Disorder == 1, na.rm = TRUE),
    
    #3.153 Calculate percentage of cases where only the parent reported Social Anxiety Disorder
    Percent_Parent_Social_Anxiety_Only = (Parent_Social_Anxiety_Only / Total_Social_Anxiety * 100),
    
    #3.154 Social Anxiety Disorder: Count cases where only the youth reported Social Anxiety Disorder
    Youth_Social_Anxiety_Only = sum(Social_Anxiety_Disorder_Source == "Youth", na.rm = TRUE),
    
    #3.155 Calculate percentage of cases where only the youth reported Social Anxiety Disorder
    Percent_Youth_Social_Anxiety_Only = (Youth_Social_Anxiety_Only / Total_Social_Anxiety * 100),
    
    #3.156 Social Anxiety Disorder: Count cases where both parent and youth reported Social Anxiety Disorder
    Concordant_Social_Anxiety = sum(Social_Anxiety_Disorder_Source == "Concordant", na.rm = TRUE),
    
    #3.157 Calculate percentage of cases where both parent and youth reported Social Anxiety Disorder
    Percent_Concordant_Social_Anxiety = (Concordant_Social_Anxiety / Total_Social_Anxiety * 100),
    
    #3.161 Separation Anxiety Disorder: Count cases where only the parent reported Separation Anxiety Disorder
    Parent_Separation_Anxiety_Only = sum(Separation_Anxiety_Disorder_Source == "Parent", na.rm = TRUE),
    
    #3.162 Count total Separation Anxiety Disorder cases
    Total_Separation_Anxiety = sum(Separation_Anxiety_Disorder == 1, na.rm = TRUE),
    
    #3.163 Calculate percentage of cases where only the parent reported Separation Anxiety Disorder
    Percent_Parent_Separation_Anxiety_Only = (Parent_Separation_Anxiety_Only / Total_Separation_Anxiety * 100),
    
    #3.164 Separation Anxiety Disorder: Count cases where only the youth reported Separation Anxiety Disorder
    Youth_Separation_Anxiety_Only = sum(Separation_Anxiety_Disorder_Source == "Youth", na.rm = TRUE),
    
    #3.165 Calculate percentage of cases where only the youth reported Separation Anxiety Disorder
    Percent_Youth_Separation_Anxiety_Only = (Youth_Separation_Anxiety_Only / Total_Separation_Anxiety * 100),
    
    #3.166 Separation Anxiety Disorder: Count cases where both parent and youth reported Separation Anxiety Disorder
    Concordant_Separation_Anxiety = sum(Separation_Anxiety_Disorder_Source == "Concordant", na.rm = TRUE),
    
    #3.167 Calculate percentage of cases where both parent and youth reported Separation Anxiety Disorder
    Percent_Concordant_Separation_Anxiety = (Concordant_Separation_Anxiety / Total_Separation_Anxiety * 100)
  ) %>%
  
  #3.17 Round all numeric columns to 2 decimal places for easier interpretation
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# 3.2 Reshape the data to be in a format suitable for Notion
Divergence_Concordance_Long <- Divergence_Concordance_Summary %>%
  pivot_longer(
    cols = -eventname,
    names_to = "temp",  # Temporarily store all name components in one column
    values_to = "value"
  ) %>%
  # Use separate to properly handle the split based on different patterns
  mutate(
    # Extract reporter (e.g., 'Parent', 'Youth', 'Concordant')
    reporter = str_extract(temp, "(?<=^)([A-Za-z]+)(?=_(GAD|MDD|Social_Anxiety|Separation_Anxiety|Concordant))"), 
    # Extract disorder (e.g., 'GAD', 'MDD', etc.)
    disorder = str_extract(temp, "(?<=_)(GAD|MDD|Social_Anxiety|Separation_Anxiety)"),
    # Extract the type (e.g., 'Only', 'Total', 'Percent')
    type = str_extract(temp, "(Only|Total|Percent)"),
    # Clean up the disorder names
    disorder = case_when(
      grepl("GAD", disorder) ~ "GAD",
      grepl("MDD", disorder) ~ "MDD",
      grepl("Social_Anxiety", disorder) ~ "Social_Anxiety",
      grepl("Separation_Anxiety", disorder) ~ "Separation_Anxiety",
      TRUE ~ disorder
    ),
    # Adjust 'reporter' handling
    reporter = case_when(
      grepl("Concordant", reporter) ~ "Concordant",  # Make sure 'Concordant' is treated as a reporter
      TRUE ~ reporter
    ),
    # Adjust 'type' handling
    type = case_when(
      grepl("Only", type) ~ "N_Only",               # Treat "Only" as N_Only
      grepl("Total", type) ~ "Total_N",              # Treat "Total" as Total_N
      grepl("Percent", type) ~ "percent",            # Treat "Percent" as percent
      TRUE ~ NA_character_
    )
  ) %>%
  # Pivot the data back into a wide format
  pivot_wider(
    names_from = type,          # Pivot wide based on 'type'
    values_from = value,        # Use 'value' for the wide columns
    values_fn = list,           # Store multiple values in list-columns
    values_fill = list(value = NA)  # Fill NAs where missing
  ) %>%
  # If needed, simplify list-columns by extracting the first value
  mutate(
    N_Only = sapply(N_Only, function(x) if (length(x) > 1) first(x) else x),
    Total_N = sapply(Total_N, function(x) if (length(x) > 1) first(x) else x),
    percent = sapply(percent, function(x) if (length(x) > 1) first(x) else x)
  ) 
  
# Consolidate rows by grouping and coalescing non-NA values within each group
Divergence_Concordance_Long <- Divergence_Concordance_Long %>% 
  group_by(eventname, reporter, disorder) %>%
  summarize(
    N_Only = coalesce(!!!N_Only),
    Total_N = coalesce(!!!Total_N),
    percent = coalesce(!!!percent),
    .groups = "drop"
  )



