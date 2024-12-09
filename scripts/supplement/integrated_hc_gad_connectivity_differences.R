## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
library(tidyr)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(car)
library(ggseg)
library(ggseg3d)
library(ggsegGordon)
library(gridGraphics) 
library(patchwork)
options(digits = 8, scipen = 999) 

# Read in required data 
# Merged GAD + HC Sample
GAD_HC_sample <- read.csv("./data_processed/supplement/GAD_HC_resampled_merged_groups.csv")

# Family ID Data
family_ID_data <-  read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(src_subject_id, rel_family_id))

# Resting state fMRI Data
cleaned_qcd_rsfMRI_data <- read.csv("./data_processed/rsfMRI_data_qcd_unsubset.csv")

# Read in demographic data
demographic_data_raw <- read.csv("./data_raw/ABCD_parent_demographic_data.csv")

# Read in race & ethnicity data
ethnicity_data <- read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, race_ethnicity))

# Read in CBCL data
cbcl_data_raw <- read.csv("./data_raw/abcd_cbcls01.csv")

# Read in GAD reporter data
integrated_gad_reporter_data <- read.csv("./data_processed/supplement/integrated_gad_reporter_data.csv")


## Data Wrangling ##

#1.1 Merge the family ID data into the GAD + HC sample
#1.11 Join the family ID data with the GAD + HC sample
GAD_HC_sample <- left_join(GAD_HC_sample, family_ID_data)

#1.12 Retain columns of interest in the desired order
GAD_HC_sample_for_merging <- GAD_HC_sample %>% 
  dplyr::select(c(src_subject_id, rel_family_id, interview_age, sex, site_name, eventname, group, GAD_timepoint, analysis_group, subgroup))

#1.2 Join the rsfMRI data with the HC + GAD sample
GAD_HC_group_connectivity_analysis_data <- left_join(GAD_HC_sample_for_merging, cleaned_qcd_rsfMRI_data)

#1.3 Prep data for analyses
#1.31 Convert relevant columns to numeric type
GAD_HC_group_connectivity_analysis_data <- mutate_at(GAD_HC_group_connectivity_analysis_data, vars(11:105), as.numeric)

#1.321 Convert relevant columns to factor type
GAD_HC_group_connectivity_analysis_data <- mutate_at(GAD_HC_group_connectivity_analysis_data, vars(c(rel_family_id, sex, site_name, eventname, group, analysis_group)), as.factor)

#1.322 Set the reference level of the group variable for analysis purposes
GAD_HC_group_connectivity_analysis_data$group <- relevel(GAD_HC_group_connectivity_analysis_data$group, ref = "HC")

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
demographic_race_ethnicity_merged_data <- left_join(GAD_HC_group_connectivity_analysis_data, demographic_ethnicity_vars, by = c("src_subject_id" = "subjectkey"))

#2.16 Create a variable for scanner model 
#2.161 Make the site name variable character type to allow for conditional logic operations
demographic_race_ethnicity_merged_data$site_name <- as.character(demographic_race_ethnicity_merged_data$site_name)

#2.162 Use conditional logic to create the scanner model variable
demographic_race_ethnicity_merged_data$scanner_model <- dplyr::case_when(
  startsWith(demographic_race_ethnicity_merged_data$site_name, "S") ~ "Siemens",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "P") ~ "Phillips Healthcare",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "G") ~ "GE Healthcare",
  TRUE ~ NA_character_)

#2.17 Create an age in years variable
demographic_race_ethnicity_merged_data$age_in_years <- floor((as.numeric(demographic_race_ethnicity_merged_data$interview_age)) / 12)

#2.18 Retain only variables of interest to the summary statistics
demographic_race_ethnicity_merged_data <- demographic_race_ethnicity_merged_data %>% 
  dplyr::select(c(src_subject_id, age_in_years, sex, eventname, group, rsfmri_c_ngd_meanmotion, race_ethnicity, scanner_model))

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

#2.3 Merge in the GAD reporter data
#2.31 Left join the 
demographic_race_ethnicity_merged_data <- left_join(demographic_race_ethnicity_merged_data, integrated_gad_reporter_data)

#2.321 Change any instances of NA to "None" (for HC subjects) for demographic summary stats generation
demographic_race_ethnicity_merged_data$GAD_Source[is.na(demographic_race_ethnicity_merged_data$GAD_Source)] <- "None"

#2.322 Change any instances of "None" to "HC" for clean demographic summary stats tabling
demographic_race_ethnicity_merged_data$GAD_Source[demographic_race_ethnicity_merged_data$GAD_Source == "None"] <- "HC"


## Analysis ##

#1. Analyze connectivity differences in metrics of interest between the current GAD and HC groups 
#1.1 Establish the range of the dependent variables 
GAD_HC_connectivity_analysis_dp_col_range <- 12:105

#1.2 Initialize an empty dataframe to store analysis values
GAD_HC_connectivity_analysis_raw_results <- data.frame(
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
for (column_number in GAD_HC_connectivity_analysis_dp_col_range) {
  
  #1.31 Get the column name
  dependent_variable <- colnames(GAD_HC_group_connectivity_analysis_data)[column_number]
  print(dependent_variable)
  
  #1.32 Fit the linear regression model
  GAD_HC_connectivity_analysis_lm <- lmerTest::lmer(GAD_HC_group_connectivity_analysis_data[, dependent_variable] ~ group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|rel_family_id), na.action = na.omit, data = GAD_HC_group_connectivity_analysis_data)
  
  #1.33 Run the ANCOVA Model on the linear regression to extract omnibus group effect(s)
  GAD_HC_connectivity_analysis_ANCOVA <- car::Anova(GAD_HC_connectivity_analysis_lm, type = "II", test.statistic = "F")
  
  #1.34 Loop through fixed effects in the model to extract relevant results
  for (independent_variable in row.names(GAD_HC_connectivity_analysis_ANCOVA)) {
    if (independent_variable == "rsfmri_c_ngd_meanmotion") {
      
      #1.341 Extract + store parameters of interest for continuous variables
      summary_lm <- summary(GAD_HC_connectivity_analysis_lm)
      estimate <- summary_lm$coefficients[independent_variable, "Estimate"]
      std_error <- summary_lm$coefficients[independent_variable, "Std. Error"]
      t_value <- summary_lm$coefficients[independent_variable, "t value"]
      p_value <- summary_lm$coefficients[independent_variable, "Pr(>|t|)"]
      GAD_HC_connectivity_analysis_raw_results <- rbind(GAD_HC_connectivity_analysis_raw_results, data.frame(
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
      f_value <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "F"]
      df <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Df"]
      residual_df <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Df.res"]
      p_value <- GAD_HC_connectivity_analysis_ANCOVA[independent_variable, "Pr(>F)"] 
      GAD_HC_connectivity_analysis_raw_results <- rbind(GAD_HC_connectivity_analysis_raw_results, data.frame(
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


## Clean Analysis Results ##

#1. Isolate the GAD vs HC group analysis results, FDR correct, determine/store significant results, and re-join with the rest of the results
#1.1 Subset the data based on the relevant DX vs CN IV variable strings
GAD_HC_connectivity_analysis_p_adjust <- subset(GAD_HC_connectivity_analysis_raw_results, grepl("group", GAD_HC_connectivity_analysis_raw_results$IV))

#1.21 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
GAD_HC_connectivity_analysis_p_adjust$p_adjusted <- p.adjust(GAD_HC_connectivity_analysis_p_adjust$p_value, method = "fdr")

#1.22 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
GAD_HC_connectivity_analysis_significant_results <- subset(GAD_HC_connectivity_analysis_p_adjust, GAD_HC_connectivity_analysis_p_adjust$p_adjusted <= 0.05)

#1.23 Specify the significant FC metrics
GAD_HC_connectivity_analysis_significant_fc_metrics <- GAD_HC_connectivity_analysis_significant_results$column_name
  
#1.31 Subset the original results dataframe to include only the rows where 'column_name' matches significant FC metrics
GAD_HC_connectivity_analysis_significant_results_full_models <- GAD_HC_connectivity_analysis_raw_results[GAD_HC_connectivity_analysis_raw_results$column_name %in% GAD_HC_connectivity_analysis_significant_fc_metrics, ]

#1.32 Merge the FDR corrected p values back into the full significant model results where applicable
GAD_HC_connectivity_analysis_significant_results_full_models <- left_join(GAD_HC_connectivity_analysis_significant_results_full_models, GAD_HC_connectivity_analysis_significant_results)

#1.33 Merge the FDR corrected p values back into the full all model results where applicable
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values <- left_join(GAD_HC_connectivity_analysis_raw_results, GAD_HC_connectivity_analysis_p_adjust)

#1.4 Pivot the full model significant results data to be wider (for copying into results tables)
#1.41 Convert necessary columns to numeric to avoid type issues
GAD_HC_connectivity_analysis_significant_results_full_models <- GAD_HC_connectivity_analysis_significant_results_full_models %>%
  mutate(across(c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted), as.numeric))

#1.42 Pivot the full model all results data to be wider (for copying into results tables)
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <- GAD_HC_connectivity_analysis_results_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = column_name, 
    names_from = IV, 
    values_from = c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted),
    names_glue = "{IV}_{.value}")

#1.5 Subset columns of interest and store them in the desired order
GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted <-
  GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted %>%
  dplyr::select(
    c(column_name,
      group_f_value,
      group_df,
      group_residual_df,
      group_p_value,
      group_p_adjusted,
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


#2. Demographic Summary Statistics
#2.1 Generate summary stats for each distinct clinical group 
#2.11 Create a dataframe containing summary statistics of interest to the paper
grouped_integrated_GAD_HC_sample_characteristics <- demographic_race_ethnicity_merged_data %>% 
  group_by(GAD_Source) %>% 
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
grouped_integrated_GAD_HC_sample_characteristics_long <- grouped_integrated_GAD_HC_sample_characteristics %>%
  pivot_longer(
    cols = -GAD_Source,
    names_to = "variable",      
    values_to = "value")

#2.122 Pivot the long version of the table back to a wider format, such that each clinical group gets their own row
grouped_integrated_GAD_HC_sample_characteristics_cleaned <- grouped_integrated_GAD_HC_sample_characteristics_long %>%
  pivot_wider(
    names_from = GAD_Source,
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
integrated_GAD_HC_sample_characteristics <- left_join(grouped_integrated_GAD_HC_sample_characteristics_cleaned, whole_sample_characteristics_long)


## Post-Hoc Contrasts ##

#1. Create violin boxplots for each of the dv of interest (nominally significant associations)
#1.1 Create a list of relevant dependent variables
post_hoc_dependent_variables <- c("rsfmri_cor_ngd_sa_scs_ptlh",
                                  "rsfmri_c_ngd_cgc_ngd_vta",
                                  "rsfmri_cor_ngd_dsa_scs_agrh",
                                  "rsfmri_cor_ngd_cerc_scs_ptlh",
                                  "rsfmri_cor_ngd_cerc_scs_cdelh",
                                  "rsfmri_c_ngd_vta_ngd_vta")

#1.2 Combine all dependent variables into a single dataframe for plotting
GAD_HC_group_connectivity_plot_data <- GAD_HC_group_connectivity_analysis_data %>%
  pivot_longer(
    cols = all_of(post_hoc_dependent_variables),
    names_to = "Dependent_Variable",
    values_to = "Value") %>% 
  dplyr::select(c(src_subject_id, eventname, group, Dependent_Variable, Value))

#1.31 Custom titles for dependent variables
GAD_HC_group_connectivity_plot_titles <- c(
  "rsfmri_cor_ngd_sa_scs_ptlh" = "SN - Left Putamen",
  "rsfmri_c_ngd_cgc_ngd_vta" = "CON- VAN",
  "rsfmri_cor_ngd_dsa_scs_agrh" = "DAN - Right Amygdala",
  "rsfmri_cor_ngd_cerc_scs_ptlh" = "CON- Left Putamen",
  "rsfmri_cor_ngd_cerc_scs_cdelh" = "CON - Left Caudate",
  "rsfmri_c_ngd_vta_ngd_vta" = "Within-VAN")

#1.32 Add a column for readable variable names
GAD_HC_group_connectivity_plot_data <- GAD_HC_group_connectivity_plot_data %>%
  mutate(Dependent_Variable_Label = GAD_HC_group_connectivity_plot_titles[Dependent_Variable])

#1.4 Create the faceted violin plots
facet_violin_plot <- ggplot(GAD_HC_group_connectivity_plot_data,   
                            aes(x = group, y = Value, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, fill = "white") +
  facet_wrap(~ Dependent_Variable_Label, scales = "free_y") +
  labs(x = "Diagnostic Group", 
       y = "Connectivity Value") +
  theme_minimal() +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), # Remove the full border
    axis.line.x = element_line(color = "black", linewidth = 1), # Add x-axis line
    axis.line.y = element_line(color = "black", linewidth = 1), # Add y-axis line
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none")


## Output ##

#1. Save the full model results as a csv file
write.csv(GAD_HC_connectivity_analysis_results_merged_adjusted_p_values_pivoted, "./results/GAD_HC_connectivity_analysis_results.csv", row.names = FALSE)

#2. Save the demographic summary stats table as a csv file
write.csv(integrated_GAD_HC_sample_characteristics, "./results/integrated_GAD_HC_sample_characteristics.csv", row.names = FALSE)

#3. Save the nominally significant grouped connectivity plot as a PNG
ggsave("./results/integrated_GAD_HC_connectivity_plot.png", plot = facet_violin_plot, dpi = 720, width = 10, height = 6, bg = "white")
