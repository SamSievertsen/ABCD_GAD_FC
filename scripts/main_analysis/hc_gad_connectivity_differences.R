## Setup ## 

# Load packages for loading, wrangling, mutating, and visualizing data
library(dplyr)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(polycor)
library(Cairo)
library(ggseg)
library(ggseg3d)
library(ggsegGordon)
library(gridGraphics) 
library(patchwork)

# Read in analysis data
site_visit_analysis_data <- read.csv("./data_processed/main_analysis/site_visit_analysis_data.csv")

# Read in demographic data
demographic_data_raw <- read.csv("./data_raw/ABCD_parent_demographic_data.csv")

# Read in race & ethnicity data
ethnicity_data_raw <- read.delim("./data_raw/acspsw03.txt") %>% 
  filter(eventname == "baseline_year_1_arm_1") %>% 
  dplyr::select(c(subjectkey, race_ethnicity))


## Data Wrangling ## 

#1. Transform all variables to the appropriate data type for modeling
#1.1 All numeric type variables
site_visit_analysis_data <- mutate_at(site_visit_analysis_data, vars(10:104), as.numeric)

#1.2 Factor type variables
#1.21 Biological sex
site_visit_analysis_data$sex <- as.factor(site_visit_analysis_data$sex)

#1.22 Timepoint
site_visit_analysis_data$eventname <- as.factor(site_visit_analysis_data$eventname)

#1.23 Scanner (site)
site_visit_analysis_data$site_name <- as.factor(site_visit_analysis_data$site_name)

#1.24 Family ID
site_visit_analysis_data$rel_family_id <- as.factor(site_visit_analysis_data$rel_family_id)

#1.25 Analysis group (HC vs GAD)
site_visit_analysis_data$analysis_group <- as.factor(site_visit_analysis_data$analysis_group)
site_visit_analysis_data$analysis_group <- relevel(site_visit_analysis_data$analysis_group, ref = "control")


## Analysis ##

#1. Run the linear model for the site-visit group using the relevant dataset. 

#1.1 Establish the range of the dependent variables 
site_visit_dp_col_range <- 11:104

#1.2 Create an empty dataframe to store analysis values
site_visit_results_df <- data.frame(
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
for (col_num in site_visit_dp_col_range) {
  
  #1.31 Get the column name
  analysis_col_name <- colnames(site_visit_analysis_data)[col_num]
  print(analysis_col_name)
  
  #1.32 Run the linear regression model
  site_visit_lm_analysis <- lmerTest::lmer(site_visit_analysis_data[, analysis_col_name] ~ analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|rel_family_id), na.action = na.omit, data = site_visit_analysis_data)
  
  #1.33 Run the ANCOVA Model
  site_visit_ANCOVA_analysis <- car::Anova(site_visit_lm_analysis, type = "II", test.statistic = "F")
  
  #1.34 Loop through fixed effects in the model
  for (iv_name in row.names(site_visit_ANCOVA_analysis)) {
    if (iv_name == "rsfmri_c_ngd_meanmotion") {
      
      #1.341 For continuous variables (rsfmri_c_ngd_meanmotion)
      summary_lm <- summary(site_visit_lm_analysis)
      estimate <- summary_lm$coefficients[iv_name, "Estimate"]
      std_error <- summary_lm$coefficients[iv_name, "Std. Error"]
      t_value <- summary_lm$coefficients[iv_name, "t value"]
      p_value <- summary_lm$coefficients[iv_name, "Pr(>|t|)"]
      site_visit_results_df <- rbind(site_visit_results_df, data.frame(
        column_name = analysis_col_name,
        IV = iv_name,
        estimate = estimate,
        std_error = std_error,
        t_value = t_value,
        f_value = NA,
        df = NA,
        residual_df = NA,
        p_value = p_value))
      
    } else {
      
      #1.342 For categorical variables (analysis_group, sex, eventname)
      f_value <- site_visit_ANCOVA_analysis[iv_name, "F"]
      df <- site_visit_ANCOVA_analysis[iv_name, "Df"]
      residual_df <- site_visit_ANCOVA_analysis[iv_name, "Df.res"]
      p_value <- site_visit_ANCOVA_analysis[iv_name, "Pr(>F)"] 
      site_visit_results_df <- rbind(site_visit_results_df, data.frame(
        column_name = analysis_col_name,
        IV = iv_name,
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

#1.41 Subset the data based on the relevant DX vs CN IV variable strings
site_visit_p_adjust_subset <- subset(site_visit_results_df, grepl("group", site_visit_results_df$IV))

#1.42 Create a new column conducting an FDR (p adjustment) on the derived p values, which controls the expected proportion of false positives among all the rejected hypotheses
site_visit_p_adjust_subset$p_adjusted <- p.adjust(site_visit_p_adjust_subset$p_value, method = "fdr")

#1.43 Create a new dataframe where only significant results from the analysis variable (DX vs CN) are stored
site_visit_significant_results <- subset(site_visit_p_adjust_subset, site_visit_p_adjust_subset$p_adjusted <= 0.05)

#1.44 Write the new results to an xlsx file
write_xlsx(site_visit_significant_results, "./results/site_visit_significant_results.xlsx")

#1.45 For creating the finalized analysis table(s), create a subset of the analysis results containing only data from FC metrics that displaying significant group differences
#1.451 Specify the significant FC metrics
significant_fc_metrics <-
  c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh")

#1.452 Subset the dataframe to include only the rows where 'column_name' matches significant FC metrics
site_visit_significant_results_full_models <- site_visit_results_df[site_visit_results_df$column_name %in% significant_fc_metrics, ]

#1.4531 Merge the FDR corrected P values back into the full significant model results where applicable
site_visit_significant_results_full_models <- left_join(site_visit_significant_results_full_models, site_visit_significant_results)

#1.4532 Merge the FDR corrected P values back into the full all model results where applicable
site_visit_results_df_merged_adjusted_p_values <- left_join(site_visit_results_df, site_visit_p_adjust_subset)

#1.4541 Pivot the full model significant results data to be wider (for copying into results tables)
# Convert necessary columns to numeric to avoid type issues
site_visit_significant_results_full_models <- site_visit_significant_results_full_models %>%
  mutate(across(c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted), as.numeric))

#1.4542 Pivot the full model all results data to be wider (for copying into results tables)
site_visit_results_df_merged_adjusted_p_values_pivoted <- site_visit_results_df_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = column_name, 
    names_from = IV, 
    values_from = c(estimate, std_error, t_value, f_value, df, residual_df, p_value, p_adjusted),
    names_glue = "{IV}_{.value}")

#1.455 Subset columns of interest and store them in the desired order
site_visit_results_df_merged_adjusted_p_values_pivoted <-
  site_visit_results_df_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(column_name,
      analysis_group_f_value,
      analysis_group_df,
      analysis_group_residual_df,
      analysis_group_p_value,
      analysis_group_p_adjusted,
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
      eventname_p_value
    )
  )

#1.455 Write the full results from both the whole analysis and significant models to an XLSX file
write_xlsx(site_visit_results_df_merged_adjusted_p_values_pivoted, "./results/site_visit_all_full_model_results.xlsx")


#2. Summary stat compilation + visualization of significant analysis results
#2.11 Store the variable names for the significant variables
rsfmri_cor_ngd_cerc_scs_cdelh	= "cingulo-opercular-network_left-caudate"
rsfmri_cor_ngd_cerc_scs_aglh = "cingulo-opercular-network_left-amygdala"
rsfmri_c_ngd_vta_ngd_vta = "within_ventral-attention-network"
rsfmri_cor_ngd_df_scs_ptlh = "default-mode-network_left-putamen"
rsfmri_cor_ngd_sa_scs_ptlh = "salience-network_left-putamen"

#2.12 Pull Descriptive stats to determine nature/direction of relationships between groups and rsfMRI connectivity stats of interest
significant_analysis_results_table <-  site_visit_analysis_data %>%
  group_by(analysis_group) %>%
  summarise(
    "mean_cingulo-opercular-network_left-caudate" = mean(as.numeric(rsfmri_cor_ngd_cerc_scs_cdelh)),
    "mean_cingulo-opercular-network_left-amygdala" = mean(as.numeric(rsfmri_cor_ngd_cerc_scs_aglh)),
    "mean_within_ventral-attention-network" = mean(as.numeric(rsfmri_c_ngd_vta_ngd_vta)),
    "mean_default-mode-network_left-putamen" = mean(as.numeric(rsfmri_cor_ngd_df_scs_ptlh)),
    "mean_salience-network_left-putamen" = mean(as.numeric(rsfmri_cor_ngd_sa_scs_ptlh)),
    "sd_cingulo-opercular-network_left-caudate" = sd(as.numeric(rsfmri_cor_ngd_cerc_scs_cdelh)),
    "sd_cingulo-opercular-network_left-amygdala" = sd(as.numeric(rsfmri_cor_ngd_cerc_scs_aglh)),
    "sd_within_ventral-attention-network" = sd(as.numeric(rsfmri_c_ngd_vta_ngd_vta)),
    "sd_default-mode-network_left-putamen" = sd(as.numeric(rsfmri_cor_ngd_df_scs_ptlh)),
    "sd_salience-network_left-putamen" = sd(as.numeric(rsfmri_cor_ngd_sa_scs_ptlh)))

#2.13 Write the significant result means to a csv file
write.csv(significant_analysis_results_table, "./results/significant_site_visit_analysis_results_table.csv", row.names = FALSE)

#2.2 Create the group difference + anatomical plots for significant variables of interest
#2.20 Create a version of the data specifically for plotting
site_visit_analysis_data_for_plotting <- site_visit_analysis_data

#2.201 Change the analysis group variable to character type for plotting purposes
site_visit_analysis_data_for_plotting$analysis_group <- as.character(site_visit_analysis_data_for_plotting$analysis_group)

#2.202 Change the control group to the HC group for plotting purposes
site_visit_analysis_data_for_plotting$analysis_group[site_visit_analysis_data_for_plotting$analysis_group == "control"] <- "HC"

#2.203 Change the analysis group back to factor type so as to set the HC group as the reference level
site_visit_analysis_data_for_plotting$analysis_group <- as.factor(site_visit_analysis_data_for_plotting$analysis_group)
site_visit_analysis_data_for_plotting$analysis_group <- relevel(site_visit_analysis_data_for_plotting$analysis_group, ref = "HC")


#2.21 CO Network - Left Amygdala Plot
gordon_atlas <- data.frame(brain_regions(gordon))
print(gordon_atlas)

#2.211 Co Network Regions
CO_network_data <- data.frame(
  region = brain_regions(gordon)[25:64],
  reg_col = brain_regions(gordon)[25:64])

#2.2121 CO Network Plot  
cingulo_opercular_network_plot <-
  ggseg(CO_network_data,
        atlas = gordon,
        colour = "black",
        size = 0.1,
        show.legend = FALSE,
        position = "stacked",
        mapping = aes(fill = reg_col)) +
  scale_fill_brain2(gordon$palette[CO_network_data$region]) +
  labs(x = "", y = "") +
  theme(text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank())

#2.2122 Save the CO Network Plot
ggsave("cingulo_opercular_network_plot.png", plot = cingulo_opercular_network_plot, bg = "white", dpi = 300)

#2.213 Amygdala Regions  
aseg_atlas <- data.frame(brain_regions(aseg))
print(aseg_atlas) 
Amygdala_data <- data.frame(
  region = brain_regions(aseg)[3],
  reg_col = brain_regions(aseg)[3],
  hemi = "left")

#2.214 Amygdala Plot
Amygdala_plot <- ggseg(Amygdala_data, atlas = aseg, 
                       colour = "black",
                       size = 0.1,
                       show.legend = FALSE,
                       view = "coronal",
                       mapping = aes(fill = reg_col)) +
  scale_fill_brain2(aseg$palette[Amygdala_data$region]) + 
  labs(x = "Left Amygdala") +
  theme(text = element_text(face = "bold"),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

#2.215 Save the Amygdala Plot
ggsave("Amygdala_plot.png", plot = Amygdala_plot, bg = "white", dpi = 300)

#2.216 Amygdala + CO Network violin plot
Amygdala_CO_Network_violin_plot <- ggplot(data = site_visit_analysis_data_for_plotting, aes(x = analysis_group, y = rsfmri_cor_ngd_cerc_scs_aglh)) +
  geom_violin(aes(fill = analysis_group)) +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("") +
  xlab("") +
  ggtitle("CON - Left Amygdala") +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none")

#2.217 Save the Amygdala + CO Network violin plot
ggsave("Amygdala_CO_Network_violin_plot.png", plot = Amygdala_CO_Network_violin_plot, bg = "white", height = 3, width = 4.50)

#2.22 CO Network - Left Caudate
#2.221 Caudate Regions  
print(aseg_atlas) 
Caudate_data <- data.frame(
  region = brain_regions(aseg)[5],
  reg_col = brain_regions(aseg)[5],
  hemi = "left")

#2.222 Caudate Plot
Caudate_plot <- ggseg(Caudate_data, atlas = aseg, 
                      colour = "black",
                      size = 0.1,
                      show.legend = FALSE,
                      view = "coronal",
                      mapping = aes(fill = reg_col)) +
  scale_fill_brain2(aseg$palette[Caudate_data$region]) + 
  theme(panel.background = element_rect(fill = "transparent")) +  
  labs(x = "Left Caudate") +
  theme(text = element_text(face = "bold"),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_blank())

#2.223 Save the caudate plot
ggsave("Caudate_plot.png", plot = Caudate_plot, bg = "white", dpi = 300)

#2.224 LH Caudate + CO Network violin plot
Caudate_CO_Network_violin_plot <- ggplot(data = site_visit_analysis_data_for_plotting, aes(x = analysis_group, y = rsfmri_cor_ngd_cerc_scs_cdelh)) +
  geom_violin(aes(fill = analysis_group)) +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none")

#2.225 Save the LH Caudate + CO Network violin plot
ggsave("Caudate_CO_Network_violin_plot.png", plot = Caudate_CO_Network_violin_plot, bg = "white", height = 3, width = 4.5, dpi = 300)

#2.23 DMN Network - Left Putamen
print(gordon_atlas)

#2.231 DMN Network Regions
DFM_network_data <- data.frame(
  region = brain_regions(gordon)[65:105],
  reg_col = brain_regions(gordon)[65:105])

#2.232 DMN Network Plot  
default_mode_network_plot <- ggseg(
  DFM_network_data,
  atlas = gordon,
  colour = "black",
  size = 0.1,
  show.legend = FALSE,
  position = "stacked",
  mapping = aes(fill = reg_col)) +
  scale_fill_brain2(gordon$palette[DFM_network_data$region]) +
  labs(x = "", y = "") +
  theme(
    text = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

#2.233 Save the DMN plot
ggsave("default_mode_network_plot.png", plot = default_mode_network_plot, bg = "white")

#2.234 Putamen Regions  
Putamen_data <- data.frame(
  region = brain_regions(aseg)[16],
  reg_col = brain_regions(aseg)[16],
  hemi = "left")

#2.2351 Putamen Plot
Putamen_plot <- ggseg(Putamen_data, atlas = aseg, 
                      colour = "black",
                      size = 0.1,
                      show.legend = FALSE,
                      view = "coronal",
                      mapping = aes(fill = reg_col)) +
  scale_fill_brain2(aseg$palette[Putamen_data$region]) + 
  theme(panel.background = element_rect(fill = "transparent")) +  
  labs(x = "", y = "") +
  theme(text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_blank())

#2.2352 Save the Putamen Plot
ggsave("Putamen_plot.png", plot = Putamen_plot, bg = "white", dpi = 300)

#2.236 LH Putamen + DMN violin plot
Putamen_DM_Network_violin_plot <- ggplot(data = site_visit_analysis_data_for_plotting, aes(x = analysis_group, y = rsfmri_cor_ngd_df_scs_ptlh)) +
  geom_violin(aes(fill = analysis_group)) +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none")

#2.237 Save the LH Putamen + DMN violin plot
ggsave("Putamen_DM_Network_violin_plot.png", plot = Putamen_DM_Network_violin_plot, bg = "white", width = 4.5, height = 3, dpi = 300)

#2.24 SA Network - Left Putamen
#2.241 SA Network Regions
SA_network_data <- data.frame(
  region = brain_regions(gordon)[218:221],
  reg_col = brain_regions(gordon)[218:221])

#2.242 SA Network Plot  
salience_network_plot <- ggseg(SA_network_data, atlas = gordon, 
                               colour = "black",
                               size = 0.1,
                               show.legend = FALSE,
                               position = "stacked",
                               mapping = aes(fill = reg_col)) +
  scale_fill_brain2(gordon$palette[SA_network_data$region]) + 
  labs(x = "", y = "") +
  theme(element_blank(),
        axis.line = element_blank(), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_blank())

#2.243 Save the SA Network plot
ggsave("salience_network_plot.png", plot = salience_network_plot, bg = "white")

#2.244 LH Putamen + SA Network violin plot
Putamen_SA_Network_violin_plot <- ggplot(data = site_visit_analysis_data_for_plotting, aes(x = analysis_group, y = rsfmri_cor_ngd_sa_scs_ptlh)) +
  geom_violin(aes(fill = analysis_group)) +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("") +
  xlab("") +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none")

#2.245 Save the LH Putamen + SA Network violin plot
ggsave("Putamen_SA_Network_violin_plot.png", plot = Putamen_SA_Network_violin_plot, bg = "white", width = 4.5, height = 3, dpi = 300)

#2.25 Within VAN
#2.251 VAN Regions
VAN_network_data <- data.frame(
  region = brain_regions(gordon)[268:290],
  reg_col = brain_regions(gordon)[268:290])

#2.252 VAN Plot  
ventral_attention_network_plot <- ggseg(VAN_network_data, atlas = gordon, 
                                        colour = "black",
                                        size = 0.1,
                                        show.legend = FALSE,
                                        position = "stacked",
                                        mapping = aes(fill = reg_col)) +
  scale_fill_brain2(gordon$palette[VAN_network_data$region]) + 
  labs(x = "", y = "") +
  theme(text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_blank())

#2.253 Save the VAN plot
ggsave("ventral_attention_network_plot.png", plot = ventral_attention_network_plot, bg = "white")

#2.254 Within VAN violin plot
Within_VAN_Network_violin_plot <- ggplot(data = site_visit_analysis_data_for_plotting, aes(x = analysis_group, y = rsfmri_c_ngd_vta_ngd_vta)) +
  geom_violin(aes(fill = analysis_group)) +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("Functional Connectivity") +
  xlab("") +
  scale_fill_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) + 
  theme_bw() +
  theme(
    axis.line = element_line(color = "black", linewidth = 1, linetype = "solid"),
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.position = "none")

#2.255 Save the Within VAN violin plot
ggsave("Within_VAN_Network_violin_plot.png", plot = Within_VAN_Network_violin_plot, bg = "white", width = 4.5, height = 3, dpi = 300)


#3. Compile demographic information for case and control groups
#3.1 Create a new data frame containing relevant summary stat variables
#3.11 Remove the first row the df to account for the descriptions of the variables in the data
demographic_Data_no_desc <- demographic_data_raw[-1,]

#3.121 Create a modified "race" and "ethnicity" column based on the raw data
demographic_Data_no_desc$race <-
  ifelse(
    rowSums(demographic_Data_no_desc[, grepl("demo_race_a_p___", names(demographic_Data_no_desc))] == 1) > 1,
    "multi_racial",
    ifelse(
      demographic_Data_no_desc$demo_race_a_p___10 == 1,
      "White",
      ifelse(
        demographic_Data_no_desc$demo_race_a_p___11 == 1,
        "Black_African_American",
        ifelse(
          demographic_Data_no_desc$demo_race_a_p___12 == 1,
          "American_Indian_Native_American",
          ifelse(
            demographic_Data_no_desc$demo_race_a_p___13 == 1,
            "Alaska_Native",
            ifelse(
              demographic_Data_no_desc$demo_race_a_p___14 == 1,
              "Native_Hawaiian",
              ifelse(
                demographic_Data_no_desc$demo_race_a_p___15 == 1,
                "Guamanian",
                ifelse(
                  demographic_Data_no_desc$demo_race_a_p___16 == 1,
                  "Samoan",
                  ifelse(
                    demographic_Data_no_desc$demo_race_a_p___17 == 1,
                    "Other_Pacific_Islander",
                    ifelse(
                      demographic_Data_no_desc$demo_race_a_p___18 == 1,
                      "Asian_Indian",
                      ifelse(
                        demographic_Data_no_desc$demo_race_a_p___19 == 1,
                        "Chinese",
                        ifelse(
                          demographic_Data_no_desc$demo_race_a_p___20 == 1,
                          "Filipino",
                          ifelse(
                            demographic_Data_no_desc$demo_race_a_p___21 == 1,
                            "Japanese",
                            ifelse(
                              demographic_Data_no_desc$demo_race_a_p___22 == 1,
                              "Korean",
                              ifelse(
                                demographic_Data_no_desc$demo_race_a_p___23 == 1,
                                "Vietnamese",
                                ifelse(
                                  demographic_Data_no_desc$demo_race_a_p___24 == 1,
                                  "Other_Asian",
                                  ifelse(
                                    demographic_Data_no_desc$demo_race_a_p___25 == 1,
                                    "Other_Race",
                                    ifelse(
                                      demographic_Data_no_desc$demo_race_a_p___77 == 1,
                                      "Refuse_To_Answer",
                                      ifelse(demographic_Data_no_desc$demo_race_a_p___99 == 1, "Unsure", NA)
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

#3.122 Change all coded race/ethnicity values to their proper names
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == 1] <- "White"
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == 2] <- "Black"
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == 3] <- "Hispanic"
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == 4] <- "Asian"
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == 5] <- "Other"
ethnicity_data_raw$race_ethnicity[ethnicity_data_raw$race_ethnicity == ""] <- "Unsure"

#3.13 Create a new dataframe containing only relevant demographic data 
demographic_vars <- demographic_Data_no_desc %>% dplyr::select(subjectkey, race)
demographic_vars <- left_join(demographic_vars, ethnicity_data_raw)

#3.14 Merge the demographic race variable with the imaging and clinical data
demo_df_raw <- merge(demographic_vars, site_visit_analysis_data, by = "subjectkey")

#3.211 Count the N unique subjects' race
Analysis_two_race_N <- demo_df_raw %>%
  group_by(race) %>%
  summarize(n = n(), 
            percent = (n/nrow(demo_df_raw))*100)

#3.212 Create ethnicity N's
#3.2121 Whole sample race/ethnicity N + percent
Analysis_two_ethnicity_N <- demo_df_raw %>%
  group_by(race_ethnicity) %>%
  summarize(n = n(), 
            percent = (n/nrow(demo_df_raw))*100)

#3.2121 Grouped race/ethnicity N + percent
Analysis_two_dx_grouped_ethnicity_N <- demo_df_raw %>%
  group_by(analysis_group, race_ethnicity) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(analysis_group) %>%
  mutate(total_in_group = sum(n),
         percent = (n / total_in_group) * 100)

#3.22 Sex
Analysis_two_sex_N <- demo_df_raw %>% 
  group_by(sex) %>% 
  summarize(n = n())

#3.23 Site
Analysis_two_site_N <- demo_df_raw %>% 
  group_by(site_name) %>% 
  summarize(n = n())

#3.24 Analysis
Analysis_two_column_N <- demo_df_raw %>% 
  group_by(analysis_group, eventname) %>% 
  summarize(n = n())

#3.3 Generate Summary Stats for the age in years of subjects at each imaging time point
#3.31 Age N
Analysis_two_age_N <- demo_df_raw %>%
  group_by(age_in_years) %>% 
  summarize(n = n())

#3.32 Age SD
Analysis_two_age_N$sd <- sd(demo_df_raw$age_in_years)

#3.33 Age Mean
Analysis_two_age_N$mean <- mean(demo_df_raw$age_in_years)