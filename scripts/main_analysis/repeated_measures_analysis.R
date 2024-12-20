## Setup ##

# Load packages for loading, wrangling, mutating, and visualizing data + set environmental variables
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(stats)
library(lme4)
library(lmerTest)
library(polycor)
library(Cairo)
library(modelbased)
library(multcomp)
library(emmeans)
library(ggsignif)
library(easystats)
options(digits = 6, scipen = 999)

# Read in repeated measures analysis data
repeated_measures_grouped_imaging_data <- read.csv("./data_processed/main_analysis/repeated_measures_grouped_imaging_data.csv")

# Read in demographic data
demographic_data_raw <- read.csv("./data_raw/ABCD_parent_demographic_data.csv")

# Read in race & ethnicity data
ethnicity_data <- read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, race_ethnicity))

# Read in CBCL data
cbcl_data_raw <- read.csv("./data_raw/abcd_cbcls01.csv")


## Data Prep ##

#1.1 Ensure relevant numerical values are numeric
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  mutate_at(vars(9:14), as.numeric)

#1.2 Ensure relevant factor values are factor type
#1.21 Assessment Timepoint
repeated_measures_grouped_imaging_data$eventname <- as.factor(repeated_measures_grouped_imaging_data$eventname)
repeated_measures_grouped_imaging_data$eventname <- relevel(repeated_measures_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")

#1.22 Analysis (sub) group
repeated_measures_grouped_imaging_data$group <- as.factor(repeated_measures_grouped_imaging_data$group)
repeated_measures_grouped_imaging_data$group <- relevel(repeated_measures_grouped_imaging_data$group, ref = "Control")

#1.23 Biological sex
repeated_measures_grouped_imaging_data$sex <- as.factor(repeated_measures_grouped_imaging_data$sex)

#1.24 Scanner (site) name
repeated_measures_grouped_imaging_data$site_name <- as.factor(repeated_measures_grouped_imaging_data$site_name)

#1.25 Subject ID
repeated_measures_grouped_imaging_data$subjectkey <- as.factor(repeated_measures_grouped_imaging_data$subjectkey)

#1.26 Family ID
repeated_measures_grouped_imaging_data$family_id <- as.factor(repeated_measures_grouped_imaging_data$family_id)


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
demographic_race_ethnicity_merged_data <- left_join(repeated_measures_grouped_imaging_data, demographic_ethnicity_vars, by = "subjectkey")

#2.16 Create a variable for scanner model 
demographic_race_ethnicity_merged_data$scanner_model <- dplyr::case_when(
  startsWith(demographic_race_ethnicity_merged_data$site_name, "S") ~ "Siemens",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "P") ~ "Phillips Healthcare",
  startsWith(demographic_race_ethnicity_merged_data$site_name, "G") ~ "GE Healthcare",
  TRUE ~ NA_character_)

#2.17 Remove variables of non-interest to the demographic stats
demographic_race_ethnicity_merged_data <- demographic_race_ethnicity_merged_data %>% 
  dplyr::select(-c(rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))

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
demographic_race_ethnicity_merged_data <- left_join(demographic_race_ethnicity_merged_data, cbcl_data_filtered, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))

#2.26 Ensure cbcl variables are numeric
demographic_race_ethnicity_merged_data$cbcl_scr_syn_internal_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_internal_t)
demographic_race_ethnicity_merged_data$cbcl_scr_syn_external_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_external_t)
demographic_race_ethnicity_merged_data$cbcl_scr_syn_anxdep_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_syn_anxdep_t)
demographic_race_ethnicity_merged_data$cbcl_scr_dsm5_anxdisord_t <- as.numeric(demographic_race_ethnicity_merged_data$cbcl_scr_dsm5_anxdisord_t)


## Data Analysis ##

#1. Run repeated measures model to test if there were significant connectivity differences between groups across the baseline and followup scans
#1.1 Run the repeated measures multilevel model for each significant fMRI metric from the GAD vs Control analysis
#1.111 CON & LH Caudate Model
CON_LH_Caudate_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_cdelh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_grouped_imaging_data)

#1.112 Summary of CON & LH Caudate Model
summary(CON_LH_Caudate_repeated_measures_MEM)

#1.113 ANCOVA of CON & LH Caudate Model
CON_LH_Caudate_repeated_measures_ANCOVA <- car::Anova(CON_LH_Caudate_repeated_measures_MEM, type = "III", test.statistic = "F")


#1.121 CON & LH Amygdala Model
CON_LH_Amygdala_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_aglh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_grouped_imaging_data)

#1.122 Summary of CON & LH Amygdala Model
summary(CON_LH_Amygdala_repeated_measures_MEM)

#1.123 ANCOVA of CON & LH Amygdala Model
CON_LH_Amygdala_repeated_measures_ANCOVA <- car::Anova(CON_LH_Amygdala_repeated_measures_MEM, type = "III", test.statistic="F")


#1.131 Within-VAN Model
Within_VAN_repeated_measures_MEM <- lmerTest::lmer(rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_grouped_imaging_data)

#1.1321 Summary of Within-VAN Model
summary(Within_VAN_repeated_measures_MEM)

#1.1322 Generate the model dashboard to check for potential confounds
easystats::model_dashboard(Within_VAN_repeated_measures_MEM)

#1.133 ANCOVA of Within-VAN Model
Within_VAN_ANCOVA <- car::Anova(Within_VAN_repeated_measures_MEM, type = "III", test.statistic="F")

#1.134 Within-VAN Post-Hoc 
#1.1341 Estimated marginal means within groups over time points
Within_VAN_EMM <- emmeans(Within_VAN_repeated_measures_MEM, pairwise ~ eventname | group)

#1.1342 Generate the emmeans contrast for comparing slopes between groups over time (AKA relationship between estimated means between timepoints)
Within_VAN_PostHoc_Contrasts <- contrast(Within_VAN_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL, adjust = "none")

#1.1343 Generate standard error derived confidence intervals for the model results 
Within_VAN_PostHoc_Contrasts_CI <- confint(Within_VAN_PostHoc_Contrasts)

#1.135 Within-VAN Plot
#1.1531 Create the plot data
Within_VAN_emm_plot_data <- as.data.frame(Within_VAN_EMM$emmeans)

#1.1532 Alter the control group name to HC to match naming conventions throughout the paper
Within_VAN_emm_plot_data$group <- as.character(Within_VAN_emm_plot_data$group)
Within_VAN_emm_plot_data$group[Within_VAN_emm_plot_data$group == "Control"] <- "HC"
Within_VAN_emm_plot_data$group <- as.factor(Within_VAN_emm_plot_data$group)
Within_VAN_emm_plot_data$group <- relevel(Within_VAN_emm_plot_data$group, ref = "HC")

#1.1533 Calculate the number of subjects in each group
subject_counts <- repeated_measures_grouped_imaging_data %>%
  group_by(group) %>%
  summarise(n = n_distinct(subjectkey))
subject_counts$group <- as.character(subject_counts$group)
subject_counts$group[subject_counts$group == "Control"] <- "HC"
subject_counts$group <- as.factor(subject_counts$group)
subject_counts$group <- relevel(subject_counts$group, ref = "HC")

#1.1534 Create a vector of labels with group names and subject counts
group_labels <- paste0(subject_counts$group, " (n =", subject_counts$n, ")")

#1.1535 Plot the Within-VAN data
ggplot(Within_VAN_emm_plot_data, aes(x = eventname, y = emmean, group = group, color = group, fill = group)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, aes(fill = group)) +
  geom_ribbon(aes(ymin = emmean - SE, ymax = emmean + SE), alpha = 0.07, linetype = "blank") +
  labs(x = "Timepoint", y = "Within-VAN Connectivity") +
  scale_color_manual(name = "Diagnostic Group", values = c("HC" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
  scale_fill_manual(name = "Diagnostic Group", values = c("HC" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
  scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.19, 0.26), breaks = seq(0.18, 0.28, by = 0.02)) +
  theme_minimal() +
  annotate(geom = "text", family = "Arial") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))

#1.1536 Save the plot as a high-resolution image (PNG + PDF format)
ggsave("./results/Within_VAN_MLM_Plot.png", dpi = 720, width = 8, height = 6, bg = "white")
ggsave("./results/Within_VAN_MLM_Plot.pdf", dpi = 720, width = 8, height = 6, bg = "white", device = "pdf")


#1.141 DFN LH Putamen Model
DFN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_df_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_grouped_imaging_data)

#1.142 Summary of DFN LH Putamen Model
summary(DFN_LH_Putamen_repeated_measures_MEM)

#1.143 ANCOVA of DFN LH Putamen Model
DFN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(DFN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")


#1.151 SN LH Putamen Model
SN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_sa_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = repeated_measures_grouped_imaging_data)

#1.152 Summary of SN LH Putamen Model
summary(SN_LH_Putamen_repeated_measures_MEM)

#1.153 ANCOVA of SN LH Putamen Model
SN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(SN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")


#1.16 Create a dataframe containing the p-values from each omnibus ANCOVA (FC_Metric ~ Group*Timepoint) and FDR correct for reporting purposes
#1.161 Create a dataframe containing relevant column names and associated p-values
repeated_measures_MEM_p_values <- data.frame(
  DV = c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh"), 
  DV_Summary = c("CON - Left Caudate",
    "CON - Left Amygdala",
    "Within-VAN",
    "DMN - Left Putamen",
    "SA - Left Putamen"),
  IV = c("DX_Timepoint",
    "DX_Timepoint",
    "DX_Timepoint",
    "DX_Timepoint",
    "DX_Timepoint"),
  p_value = c(0.2752583, 
              0.44029831, 
              0.035288, 
              0.81878568, 
              0.34250017))

#1.162 FDR Correct p-values
repeated_measures_MEM_p_values$p_adjusted <- p.adjust(repeated_measures_MEM_p_values$p_value, method = "fdr")


#2. Demographic Summary Statistics
#2.1 Generate summary stats for each distinct clinical group 
#2.11 Create a dataframe containing summary statistics of interest to the paper and that are fixed between timepoints
grouped_repeated_measures_sample_characteristics <- demographic_race_ethnicity_merged_data %>%
  group_by(group) %>%
  summarize(
    group_n = n(),
    n_female = sum(sex == "F", na.rm = TRUE),
    percent_female = (n_female / group_n) * 100,
    n_white = sum(race_ethnicity == "White", na.rm = TRUE),
    percent_white = (n_white / group_n) * 100,
    n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
    percent_black = (n_black / group_n) * 100,
    n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
    percent_hispanic = (n_hispanic / group_n) * 100,
    n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
    percent_asian = (n_asian / group_n) * 100,
    n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
    percent_other = (n_other / group_n) * 100,
    n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
    percent_ge_scanner = (n_ge_scanner / group_n) * 100,
    n_phillips_scanner = sum(scanner_model == "Phillips Healthcare", na.rm = TRUE),
    percent_phillips_scanner = (n_phillips_scanner / group_n) * 100,
    n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
    percent_siemens_scanner = (n_siemens_scanner / group_n) * 100
  )

#2.12 Create a grouped by timepoint dataframe containing relevant sample characteristics that fluctuate between timepoints
grouped_by_timepoint_repeated_measures_sample_characteristics <- demographic_race_ethnicity_merged_data %>%
  group_by(group, eventname) %>%
  summarize(
    mean_age = mean(age_in_years, na.rm = TRUE),
    sd_age = sd(age_in_years, na.rm = TRUE),
    mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
    sd_fd = sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
    mean_internalizing_cbcl = mean(cbcl_scr_syn_internal_t, na.rm = TRUE),
    sd_internalizing_cbcl = sd(cbcl_scr_syn_internal_t, na.rm = TRUE),
    mean_externalizing_cbcl = mean(cbcl_scr_syn_external_t, na.rm = TRUE),
    sd_externalizing_cbcl = sd(cbcl_scr_syn_external_t, na.rm = TRUE),
    mean_anxdep_cbcl = mean(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
    sd_anxdep_cbcl = sd(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
    mean_dsm5_anx_cbcl = mean(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE),
    sd_dsm5_anx_cbcl = sd(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE)
  )

#2.13 Combine the grouped sample characteristics
grouped_repeated_measures_sample_characteristics <- grouped_by_timepoint_repeated_measures_sample_characteristics %>%
  left_join(grouped_repeated_measures_sample_characteristics, by = "group")

#2.14 Reshape the grouped sample charactersitics to wide format
grouped_repeated_measures_sample_characteristics_wide<- grouped_repeated_measures_sample_characteristics %>%
  pivot_longer(
    cols = -c(group, eventname),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = c(group, eventname),
    values_from = Value
  )

#2.15 Reorderthe grouped sample charactersitics columns for the desired paper layout
grouped_repeated_measures_sample_characteristics_cleaned <- grouped_repeated_measures_sample_characteristics_wide %>%
  select(Variable, everything()) 


#2.2 Generate comparative summary stats for the same variables across the whole sample
#2.21 Whole sample statistics 
whole_sample_repeated_measures_sample_characteristics <- demographic_race_ethnicity_merged_data %>%
  summarize(
    group_n = n(),
    n_female = sum(sex == "F", na.rm = TRUE),
    percent_female = (n_female / group_n) * 100,
    n_white = sum(race_ethnicity == "White", na.rm = TRUE),
    percent_white = (n_white / group_n) * 100,
    n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
    percent_black = (n_black / group_n) * 100,
    n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
    percent_hispanic = (n_hispanic / group_n) * 100,
    n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
    percent_asian = (n_asian / group_n) * 100,
    n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
    percent_other = (n_other / group_n) * 100,
    n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
    percent_ge_scanner = (n_ge_scanner / group_n) * 100,
    n_phillips_scanner = sum(scanner_model == "Phillips Healthcare", na.rm = TRUE),
    percent_phillips_scanner = (n_phillips_scanner / group_n) * 100,
    n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
    percent_siemens_scanner = (n_siemens_scanner / group_n) * 100
  )

#2.22 Create a whole sample by timepoint dataframe containing relevant sample characteristics that fluctuate between timepoints
whole_sample_by_timepoint_repeated_measures_sample_characteristics <- demographic_race_ethnicity_merged_data %>%
  group_by(eventname) %>%
  summarize(
    mean_age = mean(age_in_years, na.rm = TRUE),
    sd_age = sd(age_in_years, na.rm = TRUE),
    mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
    sd_fd = sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
    mean_internalizing_cbcl = mean(cbcl_scr_syn_internal_t, na.rm = TRUE),
    sd_internalizing_cbcl = sd(cbcl_scr_syn_internal_t, na.rm = TRUE),
    mean_externalizing_cbcl = mean(cbcl_scr_syn_external_t, na.rm = TRUE),
    sd_externalizing_cbcl = sd(cbcl_scr_syn_external_t, na.rm = TRUE),
    mean_anxdep_cbcl = mean(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
    sd_anxdep_cbcl = sd(cbcl_scr_syn_anxdep_t, na.rm = TRUE),
    mean_dsm5_anx_cbcl = mean(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE),
    sd_dsm5_anx_cbcl = sd(cbcl_scr_dsm5_anxdisord_t, na.rm = TRUE)
  )

#2.23 Merge the whole sample and whole sample by timepoint summary stats
#2.231 Create a dummy coded variable for merging purposes
#2.2311 Whole sample dataframe 
whole_sample_repeated_measures_sample_characteristics$for_merging <- "merge"

#2.2312 Whole sample by timepoint dataframe
whole_sample_by_timepoint_repeated_measures_sample_characteristics$for_merging <- rep("merge")

#2.232 Left join the whole sample and whole sample by timepoint summary stats dataframes
whole_sample_repeated_measures_sample_characteristics <- full_join(whole_sample_repeated_measures_sample_characteristics, whole_sample_by_timepoint_repeated_measures_sample_characteristics) %>% 
  dplyr::select(-for_merging)

#2.233 Add a group column denoting the whole sample for merging purposes
whole_sample_repeated_measures_sample_characteristics$group <- rep("Whole_Sample")

#2.24 Reshape the whole sample charactersitics to wide format
whole_sample_repeated_measures_sample_characteristics_wide <- whole_sample_repeated_measures_sample_characteristics %>%
  pivot_longer(
    cols = -c(group, eventname),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = c(group, eventname),
    values_from = Value
  )

#2.15 Reorder the grouped sample charactersitics columns for the desired paper layout
whole_sample_repeated_measures_sample_characteristics_cleaned <- whole_sample_repeated_measures_sample_characteristics_wide %>%
  select(Variable, everything()) 


#2.3 Merge the grouped and whole sample summary stat tables
repeated_measures_sample_characteristics <- full_join(grouped_repeated_measures_sample_characteristics_cleaned, whole_sample_repeated_measures_sample_characteristics_cleaned)


## Output ## 

#1. Write the demographic summary stats for the repeated measures analysis as a csv file
write.csv(repeated_measures_sample_characteristics, "./results/repeated_measures_sample_demographic_stats.csv", row.names = FALSE)
