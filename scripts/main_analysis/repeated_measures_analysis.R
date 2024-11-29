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


## Data Prep ##

#1. Ensure relevant numerical values are numeric
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  mutate_at(vars(9:14), as.numeric)

#2. Ensure relevant factor values are factor type
#2.1 Assessment Timepoint
repeated_measures_grouped_imaging_data$eventname <- as.factor(repeated_measures_grouped_imaging_data$eventname)
repeated_measures_grouped_imaging_data$eventname <- relevel(repeated_measures_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")

#2.2 Analysis (sub) group
repeated_measures_grouped_imaging_data$group <- as.factor(repeated_measures_grouped_imaging_data$group)
repeated_measures_grouped_imaging_data$group <- relevel(repeated_measures_grouped_imaging_data$group, ref = "Control")

#2.3 Biological sex
repeated_measures_grouped_imaging_data$sex <- as.factor(repeated_measures_grouped_imaging_data$sex)

#2.4 Scanner (site) name
repeated_measures_grouped_imaging_data$site_name <- as.factor(repeated_measures_grouped_imaging_data$site_name)

#2.5 Subject ID
repeated_measures_grouped_imaging_data$subjectkey <- as.factor(repeated_measures_grouped_imaging_data$subjectkey)

#2.6 Family ID
repeated_measures_grouped_imaging_data$family_id <- as.factor(repeated_measures_grouped_imaging_data$family_id)


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
  labs(x = "Timepoint", y = "Within-VAN FC") +
  scale_color_manual(name = "Diagnostic Group", values = c("HC" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
  scale_fill_manual(name = "Diagnostic Group", values = c("HC" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
  scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.19, 0.26), breaks = seq(0.18, 0.28, by = 0.02)) +
  theme_minimal() +
  annotate(family = ) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        plot.margin = unit(c(1, 2, 1, 1), "lines"))

#1.1536 Save the plot as a high-resolution image (PNG format)
ggsave(".results/Within_VAN_MLM_Plot.png", dpi = 720, width = 8, height = 6, bg = "white")


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