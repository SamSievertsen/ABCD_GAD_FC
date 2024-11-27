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



## Data Analysis ##

#1. Run repeated measures model to test if there were significant connectivity differences between groups across the baseline and followup scans
#1.1 Run the repeated measures multilevel model for each significant fMRI metric from the GAD vs Control analysis
#1.111 CON & LH Caudate Model
CON_LH_Caudate_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_cdelh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(CON_LH_Caudate_repeated_measures_MEM)
CON_LH_Caudate_repeated_measures_ANCOVA <- car::Anova(CON_LH_Caudate_repeated_measures_MEM, type = "III", test.statistic="F")

#1.121 CON & LH Amygdala Model
CON_LH_Amygdala_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_aglh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(CON_LH_Amygdala_repeated_measures_MEM)
CON_LH_Amygdala_repeated_measures_ANCOVA <- car::Anova(CON_LH_Amygdala_repeated_measures_MEM, type = "III", test.statistic="F")

#1.131 Within-VAN Model
Within_VAN_repeated_measures_MEM <- lmerTest::lmer(rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(Within_VAN_repeated_measures_MEM)
Within_VAN_ANCOVA <- car::Anova(Within_VAN_repeated_measures_MEM, type = "III", test.statistic="F")
Within_VAN_ANCOVA

#1.132 Within-VAN Post-Hoc 
#1.1321 Estimated marginal means within groups over time points
Within_VAN_EMM <- emmeans(Within_VAN_repeated_measures_MEM, pairwise ~ eventname | group)

#1.1522 Generate the emmeans contrast for comparing slopes between groups over time (AKA relationship between estimated means between timepoints)
Within_VAN_PostHoc_Contrasts <- contrast(Within_VAN_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL, adjust = "none")

#1.153 Within-VAN Plot
#1.1531 Create the plot data
Within_VAN_emm_plot_data <- as.data.frame(Within_VAN_EMM$emmeans)

#1.1532 Alter the control group name to HC to match naming conventions throughout the paper
Within_VAN_emm_plot_data$group <- as.character(Within_VAN_emm_plot_data$group)
Within_VAN_emm_plot_data$group[Within_VAN_emm_plot_data$group == "Control"] <- "HC"
Within_VAN_emm_plot_data$group <- as.factor(Within_VAN_emm_plot_data$group)
Within_VAN_emm_plot_data$group <- relevel(Within_VAN_emm_plot_data$group, ref = "HC")
#1.1533 Calculate the number of subjects in each group
subject_counts <- analysis_five_grouped_imaging_data %>%
  group_by(group) %>%
  summarise(n = n_distinct(subjectkey))
subject_counts$group <- as.character(subject_counts$group)
subject_counts$group[subject_counts$group == "Control"] <- "HC"
subject_counts$group <- as.factor(subject_counts$group)
subject_counts$group <- relevel(subject_counts$group, ref = "HC")
#1.1534 Create a vector of labels with group names and subject counts
group_labels <- paste0(subject_counts$group, " (n =", subject_counts$n, ")")
#1.1535 Plot the data
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
# Save the plot as a high-resolution image (PNG format)
ggsave("Within_VAN_MLM_Plot.png", dpi = 720, width = 8, height = 6, bg = "white")

#1.141 DFN LH Putamen Model
DFN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_df_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(DFN_LH_Putamen_repeated_measures_MEM)
DFN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(DFN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")

#1.15 SN LH Putamen Model
SN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_sa_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(SN_LH_Putamen_repeated_measures_MEM)
SN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(SN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")

#1.18 Create a dataframe containing the p-values from each omnibus ANCOVA (FC_Metric ~ Group*Timepoint) and FDR correct for reporting purposes
#1.181 Create a dataframe containing relevant column names and associated p-values
ABCD_GAD_Repeated_Measures_MEM_P_Values <- data.frame(
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
#1.182 FDR Correct p-values
ABCD_GAD_Repeated_Measures_MEM_P_Values$p_adjusted <- p.adjust(ABCD_GAD_Repeated_Measures_MEM_P_Values$p_value, method = "fdr")



#1.25 Within-VAN Model
# Within_VAN_repeated_measures_PostHoc <- glht(Within_VAN_repeated_measures_MEM, linfct = mcp(group = "Tukey"), test = adjusted("fdr"))
# summary(Within_VAN_repeated_measures_PostHoc)
#1.2511 Generate the emmeans for obtaining marginal mean estimate differences within groups over time points
Within_VAN_EMM <- emmeans(Within_VAN_repeated_measures_MEM, pairwise ~ eventname | group, adjust = "fdr")
print(Within_VAN_EMM)
#1.2512 Generate standard error derived confidence intervals for the model results 
Within_VAN_EMM_CI <- confint(Within_VAN_EMM, adjust = "sidak")
Within_VAN_EMM_CI
#1.252 Generate the emmeans for obtaining marginal mean estimate differences between groups over time points for potential comparison
emmeans::emmeans(Within_VAN_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "none")
#1.2531 Generate the emmeans contrast for comparing slopes between groups over time (AKA relationship between estimated means between timepoints) as well as the FDR corrected p-values
Within_VAN_PostHoc_Contrasts <- contrast(Within_VAN_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL)
Within_VAN_PostHoc_Contrasts_FDR_Corrected <- contrast(Within_VAN_EMM[[1]], interaction = c("poly", "pairwise"), by = NULL, adjust = "fdr")
print(Within_VAN_PostHoc_Contrasts)
print(Within_VAN_PostHoc_Contrasts_FDR_Corrected)
#1.2533 Generate standard error derived confidence intervals for the model results 
Within_VAN_PostHoc_Contrasts_CI <- confint(Within_VAN_PostHoc_Contrasts)
Within_VAN_PostHoc_Contrasts_CI
#1.254 Generate the model dashboard to check for potential confounds
easystats::model_dashboard(Within_VAN_repeated_measures_MEM)

