#Setup: Load packages for loading, wrangling, mutating, and visualizing data
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
library(sampling)
library(rsample)
library(splitstackshape)
library(polycor)
library(Cairo)
library(circlize)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(chorddiag)
library(modelbased)
library(rstatix)
library(multcomp)
library(emmeans)
library(boot)
library(ggsignif)
library(easystats)
library(extrafont)
options(digits = 4)

# Set-up: Merge in resting state fMRI data
ABCD_rsfMRI_Data <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/ABCD_rsfMRI_Data.csv")
# Set-up: Remove the first row of the df to account for the descriptions of the variables in the pertinent diagnoses and raw data
ABCD_rsfMRI_Data_no_desc <- ABCD_rsfMRI_Data[-1,]
# Set-up: Create the site_name variable containing the scanner and site info
ABCD_rsfMRI_Data_no_desc$site_name <- substr(ABCD_rsfMRI_Data_no_desc$rsfmri_c_ngd_visitid, 1, 4)
# Set-up: Read in the Family ID data 
abcd_family_id_data <- read.delim("C:/Users/Sam Sievertsen/Desktop/ABCD/ABCC_Package_1203705_Tabulated-Data/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, rel_family_id))
# Set-up: Change the name of the family ID column
abcd_family_id_data <- abcd_family_id_data %>% rename(family_id = rel_family_id)
# Set-up: Merge the family ID data with the resting state data
ABCD_rsfMRI_Data_no_desc <- left_join(ABCD_rsfMRI_Data_no_desc, abcd_family_id_data)
# Set-up: Read in analysis four (analysis 2 group connectivity difference) subjectkey's and merge with resting state data
analysis_five_GAD_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_GAD_subjectkeys_groups.csv")
analysis_five_CN_sub_groups <- read.csv("C:/Users/Sam Sievertsen/Desktop/SamResearch/ABCD_Data/analysis_four_CN_subjectkeys_groups.csv")
analysis_five_sample_groups <- full_join(analysis_five_CN_sub_groups, analysis_five_GAD_sub_groups)
# Set-up: Join the grouped subjectkeys with the rsfMRI imaging data. Keep only the vars of interest 
analysis_five_grouped_imaging_data <- merge(analysis_five_sample_groups, ABCD_rsfMRI_Data_no_desc)
analysis_five_grouped_imaging_data$age_in_years <- floor(as.numeric(as.character(analysis_five_grouped_imaging_data$interview_age)) / 12)
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data %>% 
  dplyr::select(c(subjectkey, eventname, group, interview_age, age_in_years, sex, site_name, family_id, rsfmri_c_ngd_meanmotion, rsfmri_cor_ngd_cerc_scs_cdelh, rsfmri_cor_ngd_cerc_scs_aglh, rsfmri_c_ngd_vta_ngd_vta, rsfmri_cor_ngd_df_scs_ptlh, rsfmri_cor_ngd_sa_scs_ptlh))
# Set-up: Change the values of the "group" variable to represent their experimental identity
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "control"] <- "Control"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "baseline_GAD"] <- "GAD Remitter"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "followup_GAD"] <- "GAD Converter"
analysis_five_grouped_imaging_data$group[analysis_five_grouped_imaging_data$group == "GAD_Both"] <- "Continuous GAD"
# Set-up: make numerical values numeric
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data %>%
  mutate_at(vars(9:14), as.numeric)
# Set-up: make factor values factor type
analysis_five_grouped_imaging_data$eventname <- as.factor(analysis_five_grouped_imaging_data$eventname)
analysis_five_grouped_imaging_data$group <- as.factor(analysis_five_grouped_imaging_data$group)
analysis_five_grouped_imaging_data$sex <- as.factor(analysis_five_grouped_imaging_data$sex)
analysis_five_grouped_imaging_data$site_name <- as.factor(analysis_five_grouped_imaging_data$site_name)
analysis_five_grouped_imaging_data$subjectkey <- as.factor(analysis_five_grouped_imaging_data$subjectkey)
analysis_five_grouped_imaging_data$family_id <- as.factor(analysis_five_grouped_imaging_data$family_id)
# Set-up: Remove any rows with NA or empty values and retain subjects who have both baseline and followup scans only (to perform subtraction)
analysis_five_grouped_imaging_data <- na.omit(analysis_five_grouped_imaging_data)
analysis_five_grouped_imaging_data <- analysis_five_grouped_imaging_data[rowSums(analysis_five_grouped_imaging_data == "") == 0, ]
analysis_five_grouped_imaging_data <-  analysis_five_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(any(eventname == "baseline_year_1_arm_1") &
    any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()
# Set-up: Set the reference level of the diagnostic group to "controls" and the eventname to "baseline" 
analysis_five_grouped_imaging_data$group <- relevel(analysis_five_grouped_imaging_data$group, ref = "Control")
analysis_five_grouped_imaging_data$eventname <- relevel(analysis_five_grouped_imaging_data$eventname, ref = "baseline_year_1_arm_1")

#1. Run repeated measures model to test if there were significant connectivity differences between groups across the baseline and followup scans
#1.1 Run the repeated measures multilevel model for each significant fMRI metric from the GAD vs Control analysis
#1.111 CON & LH Caudate Model
CON_LH_Caudate_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_cdelh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(CON_LH_Caudate_repeated_measures_MEM)
CON_LH_Caudate_repeated_measures_ANCOVA <- car::Anova(CON_LH_Caudate_repeated_measures_MEM, type = "III", test.statistic="F")
#1.112 CON & LH Caudate Plot
# CON_LH_Caudate_effects_model <- effects::effect("group:eventname", CON_LH_Caudate_repeated_measures_MEM)
# plot(CON_LH_Caudate_effects_model)
# CON_LH_Caudate_effects_data <- as.data.frame(CON_LH_Caudate_effects_model)
# ggplot(CON_LH_Caudate_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "CON - LH Caudate FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

# #1.121 CON & LH Putamen Model
# CON_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey), data = analysis_five_grouped_imaging_data)
# summary(CON_LH_Putamen_repeated_measures_MEM)
# anova(CON_LH_Putamen_repeated_measures_MEM, type = "III")
# #1.122 CON & LH Putamen Plot
# CON_LH_Putamen_effects_model <- effects::effect("group:eventname", CON_LH_Putamen_repeated_measures_MEM)
# plot(CON_LH_Putamen_effects_model)
# CON_LH_Putamen_effects_data <- as.data.frame(CON_LH_Putamen_effects_model)
# ggplot(CON_LH_Putamen_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "CON - LH Putamen FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

#1.131 CON & LH Amygdala Model
CON_LH_Amygdala_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_aglh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(CON_LH_Amygdala_repeated_measures_MEM)
CON_LH_Amygdala_repeated_measures_ANCOVA <- car::Anova(CON_LH_Amygdala_repeated_measures_MEM, type = "III", test.statistic="F")
# #1.132 CON & LH Amygdala Plot
# CON_LH_Amygdala_effects_model <- effects::effect("group:eventname", CON_LH_Amygdala_repeated_measures_MEM)
# plot(CON_LH_Amygdala_effects_model)
# CON_LH_Amygdala_effects_data <- as.data.frame(CON_LH_Amygdala_effects_model)
# ggplot(CON_LH_Amygdala_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "CON - LH Amygdala FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

# #1.141 CON & RH Hippocampus Model
# CON_RH_Hippocampus_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_cerc_scs_hprh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey), data = analysis_five_grouped_imaging_data)
# summary(CON_RH_Hippocampus_repeated_measures_MEM)
# anova(CON_RH_Hippocampus_repeated_measures_MEM, type = "III")
# #1.142 CON & RH Hippocampus Plot
# CON_RH_Hippocampus_effects_model <- effects::effect("group:eventname", CON_RH_Hippocampus_repeated_measures_MEM)
# plot(CON_RH_Hippocampus_effects_model)
# CON_RH_Hippocampus_effects_data <- as.data.frame(CON_RH_Hippocampus_effects_model)
# ggplot(CON_RH_Hippocampus_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "CON - RH Hippocampus FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

#1.151 Within-VAN Model
Within_VAN_repeated_measures_MEM <- lmerTest::lmer(rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(Within_VAN_repeated_measures_MEM)
Within_VAN_ANCOVA <- car::Anova(Within_VAN_repeated_measures_MEM, type = "III", test.statistic="F")
Within_VAN_ANCOVA
#1.152 Within-VAN Post-Hoc 
#1.1521 Estimated marginal means within groups over time points
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
# ggplot(Within_VAN_emm_plot_data, aes(x = eventname, y = emmean, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = emmean - SE, ymax = emmean + SE), alpha = 0.07, linetype = "blank") +
#   labs(x = "Timepoint", y = "Within-VAN FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = group_labels) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0.19, 0.26), breaks = seq(0.18, 0.28, by = 0.02)) +  # Adjust y-axis
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())
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

#1.161 DFN LH Putamen Model
DFN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_df_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(DFN_LH_Putamen_repeated_measures_MEM)
DFN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(DFN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")
# #1.162 DFN LH Putamen Plot
# DFN_LH_Putamen_effects_model <- effects::effect("group:eventname", DFN_LH_Putamen_repeated_measures_MEM)
# plot(DFN_LH_Putamen_effects_model)
# DFN_LH_Putamen_effects_data <- as.data.frame(DFN_LH_Putamen_effects_model)
# ggplot(DFN_LH_Putamen_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "DFN - LH Putamen FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

#1.17 SN LH Putamen Model
SN_LH_Putamen_repeated_measures_MEM <- lmerTest::lmer(rsfmri_cor_ngd_sa_scs_ptlh ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id) + (1|subjectkey), data = analysis_five_grouped_imaging_data)
summary(SN_LH_Putamen_repeated_measures_MEM)
SN_LH_Putamen_repeated_measures_ANCOVA <- car::Anova(SN_LH_Putamen_repeated_measures_MEM, type = "III", test.statistic="F")
# #1.172 SN LH Putamen Plot
# SN_LH_Putamen_effects_model <- effects::effect("group:eventname", SN_LH_Putamen_repeated_measures_MEM)
# plot(SN_LH_Putamen_effects_model)
# SN_LH_Putamen_effects_data <- as.data.frame(SN_LH_Putamen_effects_model)
# ggplot(SN_LH_Putamen_effects_data, aes(x = eventname, y = fit, group = group, color = group, fill = group)) +
#   geom_line(alpha = 0.5) +
#   geom_point(size = 2, aes(fill = group)) +
#   geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.05, linetype = "blank") + 
#   labs(x = "Timepoint", y = "SN - LH Putamen FC") +
#   scale_color_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9")) +
#   scale_fill_manual(name = "Diagnostic Group", values = c("Control" = "#218380", "Continuous GAD" = "#D81159", "GAD Converter" = "#FDB833", "GAD Remitter" = "#5390D9"), labels = c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter")) +
#   scale_x_discrete(labels = c("baseline_year_1_arm_1" = "Baseline", "2_year_follow_up_y_arm_1" = "2 Year Followup"), expand = c(0.01, 0.01)) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.ticks = element_blank())

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

#1.2 Post Hoc Group Comparisons for Repeated Measures MEMs
#1.21 CON & LH Caudate Model
emmeans::lsmeans(CON_LH_Caudate_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
emmeans::lsmeans(CON_LH_Caudate_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")), adjust = "FDR")
easystats::model_dashboard(CON_LH_Caudate_repeated_measures_MEM)
# #1.22 CON & LH Putamen Model
# emmeans::lsmeans(CON_LH_Putamen_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
# emmeans::lsmeans(CON_LH_Putamen_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")))
# easystats::model_dashboard(CON_LH_Putamen_repeated_measures_MEM)
#1.23 CON & LH Amygdala Model
emmeans::lsmeans(CON_LH_Amygdala_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
emmeans::lsmeans(CON_LH_Amygdala_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")), adjust = "FDR")
easystats::model_dashboard(CON_LH_Amygdala_repeated_measures_MEM)
# #1.24 CON & RH Hippocampus Model
# emmeans::lsmeans(CON_RH_Hippocampus_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
# emmeans::lsmeans(CON_RH_Hippocampus_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")), adjust = "FDR")
# easystats::model_dashboard(CON_RH_Hippocampus_repeated_measures_MEM)

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

#1.26 DFN LH Putamen Model
emmeans::lsmeans(DFN_LH_Putamen_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
emmeans::lsmeans(DFN_LH_Putamen_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")), adjust = "FDR")
easystats::model_dashboard(DFN_LH_Putamen_repeated_measures_MEM)
#1.27 SN LH Putamen Model
emmeans::lsmeans(SN_LH_Putamen_repeated_measures_MEM, pairwise ~ group | eventname, adjust = "FDR")
emmeans::lsmeans(SN_LH_Putamen_repeated_measures_MEM, pairwise ~ group & eventname, pairs = list(c("control", "baseline_GAD", "followup_GAD", "GAD_Both")), adjust = "FDR")
easystats::model_dashboard(SN_LH_Putamen_repeated_measures_MEM)

#1.3 Bootstrapping confidence intervals for within-group point estimate difference results
# Define a function to compute estimated marginal means and post-hoc contrasts using bootstrapping
suppressMessages({bootstrap_analysis <- function(model_formula, data, n_bootstraps = 10000) {
  # Function to fit the model, compute EMMs, and post-hoc contrasts
  bootstrap_func <- function(data, indices) {
    # Initialize a flag to check if at least one subject has been sampled for each group
    all_groups_sampled <- FALSE
    # Repeat sampling until at least one subject is sampled for each group
    while (!all_groups_sampled) {
      # Sample subjects with replacement
      sampled_subjects <- sample(unique(data$subjectkey), size = 2074, replace = TRUE)
      # Create an empty dataframe to store the sampled data
      sampled_data <- data.frame(matrix(ncol = ncol(data), nrow = 0))
      # For each sampled subject, select one row for each time point
      for (subj in sampled_subjects) {
        # Select one row for baseline
        baseline_row <- data[data$subjectkey == subj & data$eventname == "baseline_year_1_arm_1", ]
        # Select one row for year 1 follow-up
        year2_followup_row <- data[data$subjectkey == subj & data$eventname == "2_year_follow_up_y_arm_1", ]
        # Combine the selected rows
        combined_row <- rbind(baseline_row, year2_followup_row)
        # Add the combined row to the sampled data
        sampled_data <- rbind(sampled_data, combined_row)
      }
      # Check if at least one subject has been sampled for each group
      all_groups_sampled <- all(c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter") %in% sampled_data$group)
      # If at least one subject is sampled for each group, proceed with model fitting
      if (all_groups_sampled) {
        # Print the number of unique values in the "group" variable
        cat("Number of unique values in 'group' variable:", length(unique(sampled_data$group)), "\n")
        group_levels <- levels(sampled_data$group)
        cat("Unique values in 'group' variable:", group_levels, "\n")
        # Print the dimensions of sampled_data
        print(dim(sampled_data))
        print(str(sampled_data))
        # Fit the model
        model <- lmerTest::lmer(model_formula, data = sampled_data)
        # Compute estimated marginal means
        emm <- emmeans(model, pairwise ~ eventname | group, adjust = "FDR")
        emm_df <- data.frame(emm$contrasts)
        emm_vector <- c("Control" = NA, "Continuous GAD" = NA, "GAD Converter" = NA, "GAD Remitter" = NA)
      } else {
        cat("Resampling required: Not all groups have at least one sampled subject.\n")
      }
    }
    # Loop through the groups in the specified order
    for (group_name in names(emm_vector)) {
      # Check if there is an estimate for the current group
      if (group_name %in% emm_df$group) {
        # Get the estimate for the current group
        estimate <- emm_df$estimate[emm_df$group == group_name]
        # Store the estimate in the appropriate position in the emm_vector
        emm_vector[group_name] <- estimate
      } else {
        # Print a message if the group is not found in the df
        message(paste("No estimate found for", group_name))
      }
    }
    emm_vector <- emm_df$estimate
    return(emm_vector)
  }
  # Perform bootstrapping
  boot_results <- boot(data, bootstrap_func, R = n_bootstraps)
  return(boot_results)
}
})
# Models and their formulas
models <- list(Within_VAN_repeated_measures_MEM = "rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey)")
# Perform bootstrap analysis for each model
suppressMessages({results <- list()
for (model_name in names(models)) {
  model_formula <- models[[model_name]]
  boot_results <- bootstrap_analysis(model_formula, analysis_five_grouped_imaging_data)
  results[[model_name]] <- boot_results
}})
results[[model_name]] <- boot_results
boot_ci <- lapply(1:4, function(i) {
  boot.ci(boot_results, index = i)
})

#1.4 Bootstrapping confidence intervals for between-group slope estimate difference results
# Define a function to compute estimated marginal means and post-hoc contrasts using bootstrapping
suppressMessages({bootstrap_analysis <- function(model_formula, data, n_bootstraps = 10000) {
  # Function to fit the model, compute EMMs, and post-hoc contrasts
  bootstrap_func <- function(data, indices) {
    # Initialize a flag to check if at least one subject has been sampled for each group
    all_groups_sampled <- FALSE
    # Repeat sampling until at least one subject is sampled for each group
    while (!all_groups_sampled) {
      # Sample subjects with replacement
      sampled_subjects <- sample(unique(data$subjectkey), size = 2074, replace = TRUE)
      # Create an empty dataframe to store the sampled data
      sampled_data <- data.frame(matrix(ncol = ncol(data), nrow = 0))
      # For each sampled subject, select one row for each time point
      for (subj in sampled_subjects) {
        # Select one row for baseline
        baseline_row <- data[data$subjectkey == subj & data$eventname == "baseline_year_1_arm_1", ]
        # Select one row for year 1 follow-up
        year2_followup_row <- data[data$subjectkey == subj & data$eventname == "2_year_follow_up_y_arm_1", ]
        # Combine the selected rows
        combined_row <- rbind(baseline_row, year2_followup_row)
        # Add the combined row to the sampled data
        sampled_data <- rbind(sampled_data, combined_row)
      }
      # Check if at least one subject has been sampled for each group
      all_groups_sampled <- all(c("Control", "Continuous GAD", "GAD Converter", "GAD Remitter") %in% sampled_data$group)
      # If at least one subject is sampled for each group, proceed with model fitting
      if (all_groups_sampled) {
        # Print the number of unique values in the "group" variable
        cat("Number of unique values in 'group' variable:", length(unique(sampled_data$group)), "\n")
        group_levels <- levels(sampled_data$group)
        cat("Unique values in 'group' variable:", group_levels, "\n")
        # Print the dimensions of sampled_data
        print(dim(sampled_data))
        print(str(sampled_data))
        # Fit the model
        model <- lmerTest::lmer(model_formula, data = sampled_data)
        # Compute estimated marginal means
        emm <- emmeans(model, pairwise ~ eventname | group, adjust = "FDR")
        # Compute post-hoc contrasts
        contrasts <- contrast(emm[[1]], interaction = c("poly", "pairwise"), by = NULL)
        # Store the contrast slope_contrast_results in a df
        contrasts_df <- data.frame(contrasts)
        # Create a named vector to store the pairwise estimates derived from the contrast
        contrasts_vector <- c("Control - Continuous GAD" = NA, "Control - GAD Converter" = NA, "Control - GAD Remitter" = NA, "Continuous GAD - GAD Converter" = NA, "Continuous GAD - GAD Remitter" = NA, "GAD Converter - GAD Remitter" = NA)
      } else {
        cat("Resampling required: Not all groups have at least one sampled subject.\n")
      }
    }
    # Loop through the groups in the specified order
    for (group_name in names(contrasts_vector)) {
      # Check if there is an estimate for the current group
      if (group_name %in% contrasts_df$group_pairwise) {
        # Get the estimate for the current group
        estimate <- contrasts_df$estimate[contrasts_df$group_pairwise == group_name]
        # Store the estimate in the appropriate position in the contrasts_vector
        contrasts_vector[group_name] <- estimate
      } else {
        # Print a message if the group is not found in the df
        message(paste("No estimate found for", group_name))
      }
    }
    contrasts_vector <- contrasts_df$estimate
    return(contrasts_vector)
  }
  # Perform bootstrapping
  boot_slope_contrast_results <- boot(data, bootstrap_func, R = n_bootstraps)
  return(boot_slope_contrast_results)
}
})
# Models and their formulas
models <- list(Within_VAN_repeated_measures_MEM = "rsfmri_c_ngd_vta_ngd_vta ~ group*eventname + rsfmri_c_ngd_meanmotion + sex + site_name + (1|subjectkey)")
# Perform bootstrap analysis for each model
suppressMessages({slope_contrast_results <- list()
for (model_name in names(models)) {
  model_formula <- models[[model_name]]
  boot_slope_contrast_results <- bootstrap_analysis(model_formula, analysis_five_grouped_imaging_data)
  slope_contrast_results[[model_name]] <- boot_slope_contrast_results
}})
# Store all bootstrapped model slope_contrast_results in the boot list object
slope_contrast_results[[model_name]] <- boot_slope_contrast_results
# Bootstrap confidence intervals based on the distribution of estimates
slope_contrast_boot_ci <- lapply(1:6, function(i) {
  boot.ci(boot_slope_contrast_results, index = i)
})


