## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggseg)
library(ggseg3d)
library(ggsegGordon)
library(gridGraphics) 
library(patchwork)
options(digits = 8, scipen = 999) 

# Read in required data 
# Read in Merged Clinical Comorbidity + HC Sample
comorbid_clinical_HC_sample <- read.csv("./data_processed/supplement/comorbidity_hc_resampled_merged_groups.csv")

# Read in Family ID Data
family_ID_data <-  read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(src_subject_id, rel_family_id))

# Read in Resting state fMRI Data
cleaned_qcd_rsfMRI_data <- read.csv("./data_processed/rsfMRI_data_qcd_unsubset.csv")

# Read in Demographic Data
demographic_data_raw <- read.csv("./data_raw/ABCD_parent_demographic_data.csv")

# Read in race & ethnicity data
ethnicity_data <- read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, race_ethnicity))

# Read in CBCL data
cbcl_data_raw <- read.csv("./data_raw/abcd_cbcls01.csv")


## Data Wrangling ##

#1. Merge the family ID data into the Clinical Comorbidity + HC sample in prep for group connectivity difference analysis
#1.1 Join the family ID data with the Clinical Comorbidity + HC sample
comorbid_clinical_HC_sample <- left_join(comorbid_clinical_HC_sample, family_ID_data)

#1.2 Retain columns of interest in the desired order
comorbid_clinical_HC_sample_for_merging <- comorbid_clinical_HC_sample %>% 
  dplyr::select(c(src_subject_id, rel_family_id, interview_age, sex, site_name, eventname, broad_clinical_group, comorbidity_group))

#1.3 Join the rsfMRI data with the HC + Clinical Comorbid sample
comorbid_clinical_HC_group_connectivity_analysis_data <- left_join(comorbid_clinical_HC_sample_for_merging, cleaned_qcd_rsfMRI_data)

#1.4 Prep data for analyses
#1.41 Convert relevant columns to numeric type
comorbid_clinical_HC_group_connectivity_analysis_data <- mutate_at(comorbid_clinical_HC_group_connectivity_analysis_data, vars(9:103), as.numeric)

#1.42 Convert relevant columns to factor type
comorbid_clinical_HC_group_connectivity_analysis_data <- mutate_at(comorbid_clinical_HC_group_connectivity_analysis_data, vars(c(rel_family_id, sex, site_name, eventname, broad_clinical_group, comorbidity_group)), as.factor)

#1.43 Set the reference level of the group variable for analysis purposes
comorbid_clinical_HC_group_connectivity_analysis_data$comorbidity_group <- relevel(comorbid_clinical_HC_group_connectivity_analysis_data$comorbidity_group, ref = "HC")


#2. Prep data for the creation of demographic tables
#2.11 Remove the first row the  to account for the descriptions of the variables in the data
demographic_data <- demographic_data_raw[-1,]

#2.12 Create a modified "race" and "ethnicity" column based on the raw data
demographic_data$race <-
  ifelse(
    rowSums(demographic_data[, grepl("demo_race_a_p___", names(demographic_data))] == 1) > 1,
    "multi_racial",
    ifelse(
      demographic_data$demo_race_a_p___10 == 1,
      "White",
      ifelse(
        demographic_data$demo_race_a_p___11 == 1,
        "Black_African_American",
        ifelse(
          demographic_data$demo_race_a_p___12 == 1,
          "American_Indian_Native_American",
          ifelse(
            demographic_data$demo_race_a_p___13 == 1,
            "Alaska_Native",
            ifelse(
              demographic_data$demo_race_a_p___14 == 1,
              "Native_Hawaiian",
              ifelse(
                demographic_data$demo_race_a_p___15 == 1,
                "Guamanian",
                ifelse(
                  demographic_data$demo_race_a_p___16 == 1,
                  "Samoan",
                  ifelse(
                    demographic_data$demo_race_a_p___17 == 1,
                    "Other_Pacific_Islander",
                    ifelse(
                      demographic_data$demo_race_a_p___18 == 1,
                      "Asian_Indian",
                      ifelse(
                        demographic_data$demo_race_a_p___19 == 1,
                        "Chinese",
                        ifelse(
                          demographic_data$demo_race_a_p___20 == 1,
                          "Filipino",
                          ifelse(
                            demographic_data$demo_race_a_p___21 == 1,
                            "Japanese",
                            ifelse(
                              demographic_data$demo_race_a_p___22 == 1,
                              "Korean",
                              ifelse(
                                demographic_data$demo_race_a_p___23 == 1,
                                "Vietnamese",
                                ifelse(
                                  demographic_data$demo_race_a_p___24 == 1,
                                  "Other_Asian",
                                  ifelse(
                                    demographic_data$demo_race_a_p___25 == 1,
                                    "Other_Race",
                                    ifelse(
                                      demographic_data$demo_race_a_p___77 == 1,
                                      "Refuse_To_Answer",
                                      ifelse(demographic_data$demo_race_a_p___99 == 1, "Unsure", NA)
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

#2.13 Change all coded race/ethnicity values to their names according to NDA
#2.131 White
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == 1] <- "White"

#2.132 Black
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == 2] <- "Black"

#2.133 Hispanic
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == 3] <- "Hispanic"

#2.134 Asian
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == 4] <- "Asian"

#2.135 Other
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == 5] <- "Other"

#2.136 Unsure
ethnicity_data$race_ethnicity[ethnicity_data$race_ethnicity == ""] <- "Unsure"

#2.141 Create a new dataframe containing only relevant demographic data 
demographic_vars <- demographic_data %>% dplyr::select(subjectkey, race)

#2.142 Join the demographic and race/ethnicity variables
demographic_ethnicity_vars <- left_join(demographic_vars, ethnicity_data)

#2.15 Merge the demographic race variables with the clinical data for use in the creation of demographic summary stat table(s)
demographic_race_ethnicity_merged_data <- left_join(comorbid_clinical_HC_sample, demographic_ethnicity_vars, by = c("src_subject_id" = "subjectkey"))

#2.16 Create a variable for scanner model 
demographic_race_ethnicity_merged_data$scanner_model <- dplyr::case_when(
  startsWith(demographic_race_ethnicity_merged_data$site_name, "S") ~ "Siemens",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "P") ~ "Phillips Healthcare",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "G") ~ "GE Healthcare",
  TRUE ~ NA_character_)

#2.17 Create an age in years variable
demographic_race_ethnicity_merged_data$age_in_years <- floor((as.numeric(demographic_race_ethnicity_merged_data$interview_age)) / 12)

#2.18 Create a mean framewise displacement variable
#2.181 Subset the imaging data to only include subject IDs, assessment timepoint, and mean motion (fd) values
mean_fd_values <- comorbid_clinical_HC_group_connectivity_analysis_data %>% 
  dplyr::select(c(src_subject_id, eventname, rsfmri_c_ngd_meanmotion))

#2.182 Merge in the mean fd values to the demographic data
demographic_race_ethnicity_merged_data <- left_join(demographic_race_ethnicity_merged_data, mean_fd_values)

#2.2 Clean the CBCL data
#2.21 Remove the description row from the cbcl data
cbcl_data <- cbcl_data_raw[-1,]

#2.22 Remove all instances of NA, empty, or missing questions going into the calculation of scores of interest
cbcl_data <- cbcl_data %>% 
  filter(!is.na(cbcl_scr_dsm5_anxdisord_t) & cbcl_scr_dsm5_anxdisord_t != "" &
      !is.na(cbcl_scr_syn_anxdep_t) & cbcl_scr_syn_anxdep_t != "" &
      !is.na(cbcl_scr_syn_internal_t) & cbcl_scr_syn_internal_t != "" &
      !is.na(cbcl_scr_syn_external_t) & cbcl_scr_syn_external_t != "" &
      cbcl_scr_dsm5_anxdisord_nm == 0 &
      cbcl_scr_syn_anxdep_nm == 0 &
      cbcl_scr_syn_internal_nm == 0 &
      cbcl_scr_syn_external_nm == 0)

#2.23 Retain only baseline and 2 year follow up rows
cbcl_data_filtered <- cbcl_data %>% 
  filter(eventname == "baseline_year_1_arm_1" | eventname == "2_year_follow_up_y_arm_1")

#2.24 Retain only columns of interest
cbcl_data_filtered <- cbcl_data_filtered %>% 
  dplyr::select(c(src_subject_id, eventname, cbcl_scr_syn_internal_t, cbcl_scr_syn_external_t, cbcl_scr_syn_anxdep_t, cbcl_scr_dsm5_anxdisord_t))

#2.25 Merge the cbcl data of interest into the demographic data for use in constructing the summary stats tables
demographic_race_ethnicity_merged_data <- left_join(demographic_race_ethnicity_merged_data, cbcl_data_filtered)

#2.26 Ensure cbcl variables are numeric
demographic_race_ethnicity_merged_data$cbcl_scr_syn_internal_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_internal_t)
demographic_race_ethnicity_merged_data$cbcl_scr_syn_external_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_external_t)
demographic_race_ethnicity_merged_data$cbcl_scr_syn_anxdep_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_anxdep_t)
demographic_race_ethnicity_merged_data$cbcl_scr_dsm5_anxdisord_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_dsm5_anxdisord_t)


## Data Analysis ##

#1. Analyze connectivity differences in metrics of interest between the current Clinical Comorbid and HC groups 
#1.1 Establish the range of the dependent variables 
clinical_comorbid_HC_connectivity_analysis_dp_col_range <- 10:103

#1.2 Initialize an empty dataframe to store analysis values
clinical_comorbid_HC_connectivity_analysis_raw_results <- data.frame(
  column_name = character(), 
  IV = character(),
  estimate = numeric(), 
  std_error = numeric(), 
  t_value = numeric(), 
  f_value = numeric(), 
  df = numeric(), 
  residual_df = numeric(), 
  p_value = numeric(), 
  stringsAsFactors = FALSE)

#1.3 Run the linear model for the site-visit group 
for (column_number in clinical_comorbid_HC_connectivity_analysis_dp_col_range) {
  
  #1.31 Get the column name
  dependent_variable <- colnames(comorbid_clinical_HC_group_connectivity_analysis_data)[column_number]
  print(dependent_variable)
  
  #1.32 Fit the linear regression model
  clinical_comorbid_HC_connectivity_analysis_lm <- lmerTest::lmer(comorbid_clinical_HC_group_connectivity_analysis_data[, dependent_variable] ~ comorbidity_group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|rel_family_id), na.action = na.omit, data = comorbid_clinical_HC_group_connectivity_analysis_data)
  
  #1.33 Run the ANCOVA Model on the linear regression to extract omnibus group effect(s)
  clinical_comorbid_HC_connectivity_analysis_ANCOVA <- car::Anova(clinical_comorbid_HC_connectivity_analysis_lm, type = "II", test.statistic = "F")
  
  #1.34 Loop through fixed effects in the model to extract relevant results
  for (independent_variable in row.names(clinical_comorbid_HC_connectivity_analysis_ANCOVA)) {
    if (independent_variable == "rsfmri_c_ngd_meanmotion") {
      
      #1.341 Extract + store parameters of interest for continuous variables
      summary_lm <- summary(clinical_comorbid_HC_connectivity_analysis_lm)
      estimate <- summary_lm$coefficients[independent_variable, "Estimate"]
      std_error <- summary_lm$coefficients[independent_variable, "Std. Error"]
      t_value <- summary_lm$coefficients[independent_variable, "t value"]
      p_value <- summary_lm$coefficients[independent_variable, "Pr(>|t|)"]
      clinical_comorbid_HC_connectivity_analysis_raw_results <- rbind(clinical_comorbid_HC_connectivity_analysis_raw_results, data.frame(
        column_name = dependent_variable,
        IV = independent_variable,
        estimate = estimate,
        std_error = std_error,
        t_value = t_value,
        f_value = NA,
        df = NA,
        residual_df = NA,
        p_value = p_value))
    } else {
      
      #1.342 Extract + store parameters of interest for categorical variables
      f_value <- clinical_comorbid_HC_connectivity_analysis_ANCOVA[independent_variable, "F"]
      df <- clinical_comorbid_HC_connectivity_analysis_ANCOVA[independent_variable, "Df"]
      residual_df <- clinical_comorbid_HC_connectivity_analysis_ANCOVA[independent_variable, "Df.res"]
      p_value <- clinical_comorbid_HC_connectivity_analysis_ANCOVA[independent_variable, "Pr(>F)"] 
      clinical_comorbid_HC_connectivity_analysis_raw_results <- rbind(clinical_comorbid_HC_connectivity_analysis_raw_results, data.frame(
        column_name = dependent_variable,
        IV = independent_variable,
        estimate = NA,
        std_error = NA,
        t_value = NA,
        f_value = f_value,
        df = df,
        residual_df = residual_df,
        p_value = p_value
      ))
    }
  }
}


#2. Demographic Summary Statistics
#2.1 Generate summary stats for each distinct clinical group 
#2.11 Create a dataframe containing summary statistics of interest to the paper
grouped_clinical_comorbid_HC_sample_characteristics <- demographic_race_ethnicity_merged_data %>% 
  group_by(comorbidity_group) %>% 
  summarize(group_n = n(),
            mean_age = mean(age_in_years, na.rm = TRUE),
            sd_age = sd(age_in_years),
            n_female = sum(sex == "F", na.rm = TRUE),
            percent_female = ((n_female/group_n) * 100),
            mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
            sd_fd = sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
            n_white = sum(race_ethnicity == "White", na.rm = TRUE),
            percent_white = ((n_white/group_n) * 100),
            n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
            percent_black = ((n_black/group_n) * 100),
            n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
            percent_hispanic = ((n_hispanic/group_n) * 100),
            n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
            percent_asian = ((n_asian/group_n) * 100),
            n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
            percent_other = ((n_other/group_n) * 100),
            n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
            percent_ge_scanner = ((n_ge_scanner/group_n) * 100),
            n_phillips_scanner = sum(scanner_model == "Phillips Healthcare", na.rm = TRUE),
            percent_phillips_scanner = ((n_phillips_scanner/group_n) * 100),
            n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
            percent_siemens_scanner = ((n_siemens_scanner/group_n) * 100),
            mean_internalizing_cbcl = mean(cbcl_scr_syn_internal_t, na.rm = TRUE),
            sd_internalizing_cbcl = sd(cbcl_scr_syn_internal_t, na.rm = TRUE),
            mean_externalizing_cbcl = mean(cbcl_scr_syn_external_t, na.rm = TRUE),
            sd_externalizing_cbcl = sd(cbcl_scr_syn_external_t, na.rm = TRUE),
            mean_anxdep_cbcl = mean(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
            sd_anxdep_cbcl = sd(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
            mean_dsm5_anx_cbcl = mean(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE),
            sd_dsm5_anx_cbcl = sd(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE)
            )

#2.121 Pivot the demographic summary stats into a long version of the table
grouped_clinical_comorbid_HC_sample_characteristics_long <- grouped_clinical_comorbid_HC_sample_characteristics %>%
  pivot_longer(
    cols = -comorbidity_group,
    names_to = "variable",      
    values_to = "value")

#2.122 Pivot the long version of the table back to a wider format, such that each clinical group gets their own row
grouped_clinical_comorbid_HC_sample_characteristics_cleaned <- grouped_clinical_comorbid_HC_sample_characteristics_long %>%
  pivot_wider(
    names_from = comorbidity_group,
    values_from = value)


#2.2 Generate comparative summary stats for the same variables across the whole sample
#2.21 Create a dataframe containing summary statistics of interest to the paper
whole_sample_characteristics <- demographic_race_ethnicity_merged_data %>% 
  summarize(total_n = n(),
            mean_age = mean(age_in_years, na.rm = TRUE),
            sd_age = sd(age_in_years),
            n_female = sum(sex == "F", na.rm = TRUE),
            percent_female = ((n_female/total_n) * 100),
            mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
            sd_fd = sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
            n_white = sum(race_ethnicity == "White", na.rm = TRUE),
            percent_white = ((n_white/total_n) * 100),
            n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
            percent_black = ((n_black/total_n) * 100),
            n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
            percent_hispanic = ((n_hispanic/total_n) * 100),
            n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
            percent_asian = ((n_asian/total_n) * 100),
            n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
            percent_other = ((n_other/total_n) * 100),
            n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
            percent_ge_scanner = ((n_ge_scanner/total_n) * 100),
            n_phillips_scanner = sum(scanner_model == "Phillips Healthcare", na.rm = TRUE),
            percent_phillips_scanner = ((n_phillips_scanner/total_n) * 100),
            n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
            percent_siemens_scanner = ((n_siemens_scanner/total_n) * 100),
            mean_internalizing_cbcl = mean(cbcl_scr_syn_internal_t, na.rm = TRUE),
            sd_internalizing_cbcl = sd(cbcl_scr_syn_internal_t, na.rm = TRUE),
            mean_externalizing_cbcl = mean(cbcl_scr_syn_external_t, na.rm = TRUE),
            sd_externalizing_cbcl = sd(cbcl_scr_syn_external_t, na.rm = TRUE),
            mean_anxdep_cbcl = mean(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
            sd_anxdep_cbcl = sd(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
            mean_dsm5_anx_cbcl = mean(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE),
            sd_dsm5_anx_cbcl = sd(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE)
  )

#2.221 Pivot the whole sample summary stats into a long version of the table
whole_sample_characteristics_long <- whole_sample_characteristics %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "whole_sample")

#2.3 Merge the grouped and whole sample summary stat tables
clinical_comorbid_HC_sample_characteristics <- left_join(grouped_clinical_comorbid_HC_sample_characteristics_cleaned, whole_sample_characteristics_long)


## Clean Analysis Results ##

#1. Isolate the Clinical Comorbid vs HC group analysis results, FDR correct, determine/store significant results, and re-join with the rest of the results
#1.1 Subset the data based on the relevant DX vs CN IV variable strings
clinical_comorbid_HC_connectivity_analysis_p_adjust <- subset(clinical_comorbid_HC_connectivity_analysis_raw_results, grepl("group", clinical_comorbid_HC_connectivity_analysis_raw_results$IV))

#1.21 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
clinical_comorbid_HC_connectivity_analysis_p_adjust$p_adjusted <- p.adjust(clinical_comorbid_HC_connectivity_analysis_p_adjust$p_value, method = "fdr")

#1.22 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
clinical_comorbid_HC_connectivity_analysis_significant_results <- subset(clinical_comorbid_HC_connectivity_analysis_p_adjust, clinical_comorbid_HC_connectivity_analysis_p_adjust$p_adjusted <= 0.05)

#1.23 Specify the significant FC metrics
clinical_comorbid_HC_connectivity_analysis_significant_fc_metrics <- clinical_comorbid_HC_connectivity_analysis_significant_results$column_name

#1.31 Subset the original results dataframe to include only the rows where 'column_name' matches significant FC metrics
clinical_comorbid_HC_connectivity_analysis_significant_results_full_models <- clinical_comorbid_HC_connectivity_analysis_raw_results[clinical_comorbid_HC_connectivity_analysis_raw_results$column_name %in% clinical_comorbid_HC_connectivity_analysis_significant_fc_metrics, ]

#1.32 Merge the FDR corrected p values back into the full significant model results where applicable
clinical_comorbid_HC_connectivity_analysis_significant_results_full_models <- left_join(clinical_comorbid_HC_connectivity_analysis_significant_results_full_models, clinical_comorbid_HC_connectivity_analysis_significant_results)

#1.33 Merge the FDR corrected p values back into the full all model results where applicable
clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values <- left_join(clinical_comorbid_HC_connectivity_analysis_raw_results, clinical_comorbid_HC_connectivity_analysis_p_adjust)

#1.4 Pivot the full model significant results data to be wider (for copying into results tables)
#1.41 Convert necessary columns to numeric to avoid type issues
clinical_comorbid_HC_connectivity_analysis_significant_results_full_models <- clinical_comorbid_HC_connectivity_analysis_significant_results_full_models %>%
  mutate(across(c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted), as.numeric))

#1.42 Pivot the full model all results data to be wider (for copying into results tables)
clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <- clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = column_name, 
    names_from = IV, 
    values_from = c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted),
    names_glue = "{IV}_{.value}")

#1.5 Subset columns of interest and store them in the desired order
clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <-
  clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted %>%
  dplyr::select(
    c(column_name,
      comorbidity_group_f_value,
      comorbidity_group_df,
      comorbidity_group_residual_df,
      comorbidity_group_p_value,
      comorbidity_group_p_adjusted,
      rsfmri_c_ngd_meanmotion_estimate,
      rsfmri_c_ngd_meanmotion_t_value,
      rsfmri_c_ngd_meanmotion_std_error,
      rsfmri_c_ngd_meanmotion_p_value,
      sex_f_value,
      sex_df,
      sex_residual_df,
      sex_p_value,
      eventname_f_value,
      eventname_df,
      eventname_residual_df,
      eventname_p_value))


## Post-Hoc Contrasts ##

#1. Run pairwise EMM contrasts on each level of the comorbidity_group variable for nominally significant associations
#1.1 Create a list of relevant dependent variables
post_hoc_dependent_variables <- c("rsfmri_c_ngd_vta_ngd_vta", 
                                  "rsfmri_cor_ngd_cerc_scs_ptlh",
                                  "rsfmri_cor_ngd_df_scs_cderh",
                                  "rsfmri_cor_ngd_df_scs_aarh",
                                  "rsfmri_c_ngd_cgc_ngd_fo",
                                  "rsfmri_cor_ngd_cerc_scs_aglh",
                                  "rsfmri_cor_ngd_cerc_scs_hprh",
                                  "rsfmri_cor_ngd_vta_scs_aalh",
                                  "rsfmri_cor_ngd_cerc_scs_cdelh",
                                  "rsfmri_cor_ngd_sa_scs_ptlh",
                                  "rsfmri_c_ngd_dla_ngd_vta")

#1.2 Create a dataframe to store results
emmeans_results <- list()

#1.3 Iterate through each dependent variable
for (dv in post_hoc_dependent_variables) {
  
  #1.31 Fit the model for the current dependent variable
  post_hoc_model <- lmerTest::lmer(
    comorbid_clinical_HC_group_connectivity_analysis_data[, dv] ~ 
      comorbidity_group + rsfmri_c_ngd_meanmotion + sex + eventname + 
      (1 | site_name) + (1 | rel_family_id), 
    na.action = na.omit, 
    data = comorbid_clinical_HC_group_connectivity_analysis_data
  )
  
  #1.32 Compute estimated marginal means for the comorbidity group
  emm <- emmeans(post_hoc_model, ~ comorbidity_group)
  
  #1.33 Perform pairwise contrasts
  pairwise_contrasts <- contrast(emm, method = "pairwise", adjust = "fdr")
  
  #1.34 Store results
  emmeans_results[[dv]] <- list(
    "EMMs" = emm,
    "Pairwise_Contrasts" = pairwise_contrasts
  )
  
  #1.35 Print summary of results for the current DV
  print(paste("Results for:", dv))
  
  #1.36 Print the estimated marginal means
  print("Estimated Marginal Means:")
  print(summary(emm))
  
  #1.37 Print the pairwise contrasts
  print("Pairwise Contrasts:")
  print(summary(pairwise_contrasts))
  
}

#2. Create violin boxplots for each of the dv of interest (nominally significant associations)
#2.1 Combine all dependent variables into a single dataframe for plotting
comorbid_clinical_HC_group_connectivity_plot_data <- comorbid_clinical_HC_group_connectivity_analysis_data %>%
  pivot_longer(
    cols = all_of(post_hoc_dependent_variables),
    names_to = "Dependent_Variable",
    values_to = "Value") %>% 
  dplyr::select(c(src_subject_id, eventname, comorbidity_group, Dependent_Variable, Value))

#2.21 Custom titles for dependent variables
comorbid_clinical_HC_group_plot_titles <- c(
  "rsfmri_c_ngd_vta_ngd_vta" = "Within-VAN",
  "rsfmri_cor_ngd_cerc_scs_ptlh" = "CON - Left Putamen",
  "rsfmri_cor_ngd_df_scs_cderh" = "DMN - Right Caudate",
  "rsfmri_cor_ngd_df_scs_aarh" = "DMN - Right Accumbens",
  "rsfmri_c_ngd_cgc_ngd_fo" = "CON - FPN",
  "rsfmri_cor_ngd_cerc_scs_aglh" = "CON - Left Amygdala",
  "rsfmri_cor_ngd_cerc_scs_hprh" = "CON - Right Hippocampus",
  "rsfmri_cor_ngd_vta_scs_aalh" = "VAN - Left Accumbens",
  "rsfmri_cor_ngd_cerc_scs_cdelh" = "CON - Left Caudate",
  "rsfmri_cor_ngd_sa_scs_ptlh" = "SN - Left Putamen",
  "rsfmri_c_ngd_dla_ngd_vta" = "DAN - VAN")

#2.22 Add a column for readable variable names
comorbid_clinical_HC_group_connectivity_plot_data <- comorbid_clinical_HC_group_connectivity_plot_data %>%
  mutate(Dependent_Variable_Label = comorbid_clinical_HC_group_plot_titles[Dependent_Variable])

#2.22 Add a column for readable variable names
comorbid_clinical_HC_group_connectivity_plot_data <- comorbid_clinical_HC_group_connectivity_plot_data %>%
  mutate(Dependent_Variable_Label = comorbid_clinical_HC_group_plot_titles[Dependent_Variable])

#2.3 Create custom names for each of the clinical comorbidity groups for final plot
#2.31 Change the comorbidity group to character type to be able to alter values
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group <- as.character(comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group)

#2.321 GAD Only
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group[comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group == "GAD_Only"] <- "GAD Only"

#2.322 GAD + Any Comorbidity
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group[comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group == "GAD_with_comorbidities"] <- "GAD + Any Comorbidity"

#2.323 MDD no GAD
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group[comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group == "MDD_Only"] <- "MDD no GAD"

#2.324 Separation Anx no GAD
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group[comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group == "Separation_Anxiety_Only"] <- "Separation Anx no GAD"

#2.325 Social Anx no GAD
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group[comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group == "Social_Anxiety_Only"] <- "Social Anx no GAD"

#2.33 Change the comorbidity group variable back to factor type and set the HC group as the reference level
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group <- as.factor(comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group)
comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group <- relevel(comorbid_clinical_HC_group_connectivity_plot_data$comorbidity_group, ref = "HC")

#2.41 Extended significance annotations for multiple facets
significance_data <- data.frame(
  Dependent_Variable_Label = c(
    # Within-VAN
    "Within-VAN", "Within-VAN", "Within-VAN", "Within-VAN",
    # CON - Left Putamen
    "CON - Left Putamen", "CON - Left Putamen", "CON - Left Putamen",
    # DMN - Right Caudate
    "DMN - Right Caudate", "DMN - Right Caudate", "DMN - Right Caudate",
    # DMN - Right Accumbens
    "DMN - Right Accumbens", "DMN - Right Accumbens", "DMN - Right Accumbens",
    "DMN - Right Accumbens", "DMN - Right Accumbens", "DMN - Right Accumbens",
    "DMN - Right Accumbens", "DMN - Right Accumbens",
    # CON - FPN
    "CON - FPN", "CON - FPN", "CON - FPN"
  ),
  x_start = c(
    # Within-VAN
    "GAD Only", "GAD Only", "GAD + Any Comorbidity", "GAD + Any Comorbidity",
    # CON - Left Putamen
    "HC", "GAD Only", "MDD no GAD",
    # DMN - Right Caudate
    "HC", "GAD Only", "GAD + Any Comorbidity",
    # DMN - Right Accumbens
    "HC", "HC", "GAD Only", "GAD Only", "GAD + Any Comorbidity",
    "GAD + Any Comorbidity", "MDD no GAD", "Separation Anx no GAD",
    # CON - FPN
    "HC", "HC", "MDD no GAD"
  ),
  x_end = c(
    # Within-VAN
    "MDD no GAD", "Separation Anx no GAD", "MDD no GAD", "Separation Anx no GAD",
    # CON - Left Putamen
    "Social Anx no GAD", "MDD no GAD", "Social Anx no GAD",
    # DMN - Right Caudate
    "MDD no GAD", "MDD no GAD", "MDD no GAD",
    # DMN - Right Accumbens
    "MDD no GAD", "Separation Anx no GAD", "MDD no GAD", "Separation Anx no GAD",
    "MDD no GAD", "Separation Anx no GAD", "Social Anx no GAD", "Social Anx no GAD",
    # CON - FPN
    "MDD no GAD", "Social Anx no GAD", "Social Anx no GAD"
  ),
  y_position = c(
    # Within-VAN
    0.4, 0.45, 0.5, 0.55,
    # CON - Left Putamen
    0.55, 0.6, 0.65,
    # DMN - Right Caudate
    0.4, 0.45, 0.5,
    # DMN - Right Accumbens
    0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
    # CON - FPN
    0.4, 0.45, 0.5
  ),
  label = ""
)

#2.42 Calculate dynamic y_position for significance bars with offsets
max_values <- comorbid_clinical_HC_group_connectivity_plot_data %>% 
  group_by(Dependent_Variable_Label) %>% 
  summarise(max_value = max(Value, na.rm = TRUE), .groups = "drop")

#2.43 Adjust significance data for offsets in the DMN - Right Accumbens & CON - Left Putamen
significance_data <- significance_data %>% 
  left_join(max_values, by = "Dependent_Variable_Label") %>% 
  ungroup() %>% 
  group_by(Dependent_Variable_Label) %>% 
  mutate(
    # For specific problematic groups, dynamically offset significance bars
    y_position = pmin(
      max_value + 0.06 + (row_number() - 1) * 0.08, # Increment offsets for visual clarity
      0.8  # Constrain bars within visualization range
    )
  )

#2.51 Create the faceted violin plots
# Create the faceted violin plots without significance bars
facet_violin_plot <- ggplot(comorbid_clinical_HC_group_connectivity_plot_data,   
                            aes(x = comorbidity_group, y = Value, fill = comorbidity_group)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) +
  facet_wrap(~ Dependent_Variable_Label, scales = "free_y") +
  labs(x = "Comorbidity Group", 
       y = "Connectivity Value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(-0.6, 0.8))


#2.52 Dynamically add significance bars for each facet group
for (i in unique(significance_data$Dependent_Variable_Label)) {
  current_significance <- subset(significance_data, Dependent_Variable_Label == i)
  facet_violin_plot <- facet_violin_plot +
    geom_signif(
      data = current_significance,
      aes(xmin = x_start, xmax = x_end, annotations = label, y_position = y_position),
      manual = TRUE,
      inherit.aes = FALSE
    )
}


## Output ##

#1. Save the full model results as a csv file
write.csv(clinical_comorbid_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted, "./results/clinical_comorbid_HC_connectivity_analysis_results.csv", row.names = FALSE)

#2. Write the demographic summary stats as a csv file
write.csv(clinical_comorbid_HC_sample_characteristics, "./results/clinical_comorbid_HC_connectivity_grouped_demographic_stats.csv", row.names = FALSE)

#3. Save the nominally significant grouped connectivity plot as a PDF
