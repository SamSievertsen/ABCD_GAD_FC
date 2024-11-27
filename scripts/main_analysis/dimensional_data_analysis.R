## Setup ##

# Load packages for loading, wrangling, mutating, and visualizing data and set environmental variables
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(stringr)
library(data.table)
library(tidyr)
library(purrr)
library(lme4)
library(lmerTest)
library(polycor)
library(Cairo)
library(viridis)
library(patchwork)
library(modelbased)
library(cowplot)
suppressWarnings({library(gridExtra)})
suppressWarnings({library(sjPlot)})
suppressWarnings({library(sjmisc)})
suppressWarnings({library(sjlabelled)})
options(scipen = 999, digits = 6)

# Read in analysis data
dimensional_analysis_data <- read.csv("./data_processed/main_analysis/dimensional_analysis_data.csv")


## Data Wrangling ##

#1. Change data values to their required type for analysis
#1.1 Change rsfMRI and CBCL/BPM vaues to numeric 
dimensional_analysis_data <- dimensional_analysis_data %>%
  mutate_at(vars(10:21), as.numeric)

#1.2 Factor type values

#1.21 Subject ID
dimensional_analysis_data$subjectkey <- as.factor(dimensional_analysis_data$subjectkey)

#1.22 Family ID
dimensional_analysis_data$family_id <- as.factor(dimensional_analysis_data$family_id)

#1.23 Sex
dimensional_analysis_data$sex <- as.factor(dimensional_analysis_data$sex)

#1.24 Assessment Timepoint
dimensional_analysis_data$eventname <- as.factor(dimensional_analysis_data$eventname)

#1.25 Scanner site (site name)
dimensional_analysis_data$site_name <- as.factor(dimensional_analysis_data$site_name)

#1.26 Analysis group
dimensional_analysis_data$analysis_group <- as.factor(dimensional_analysis_data$analysis_group)


#2. Create additional dataframes for GAD specific and HC specific analyses
#2.1 GAD group data
#2.11 Subset the dimensional analysis data to just include the GAD group
gad_group_dimensional_analysis_data <- subset(dimensional_analysis_data, analysis_group == "GAD")

#2.12 Ensure the relevant variables are the correct data type
#2.121 Assessment timepoint
gad_group_dimensional_analysis_data$eventname <- as.factor(gad_group_dimensional_analysis_data$eventname)

#2.122 Sex
gad_group_dimensional_analysis_data$sex <- as.factor(gad_group_dimensional_analysis_data$sex)


#2.2 HC group data
#2.21 Subset the dimensional analysis data to just include the HC group
hc_group_dimensional_analysis_data <- subset(dimensional_analysis_data, analysis_group == "control")

#2.12 Ensure the relevant variables are the correct data type
#2.121 Assessment timepoint
hc_group_dimensional_analysis_data$eventname <- as.factor(hc_group_dimensional_analysis_data$eventname)

#2.122 Sex
hc_group_dimensional_analysis_data$sex <- as.factor(hc_group_dimensional_analysis_data$sex)


#3. Create a copy of the imaging CBCL data for plotting 
dimensional_analysis_data_for_plotting <- dimensional_analysis_data

#3.11 Alter the analysis group variable to match paper formatting
#3.12 Change the variable to character type
dimensional_analysis_data_for_plotting$analysis_group <- as.character(dimensional_analysis_data_for_plotting$analysis_group)

#3.13 Change the "control" values to "HC" values for aesthetic purposes
dimensional_analysis_data_for_plotting$analysis_group[dimensional_analysis_data_for_plotting$analysis_group == "control"] <- "HC"

#3.14 Change the analysis group variable back to factor type
dimensional_analysis_data_for_plotting$analysis_group <- as.factor(dimensional_analysis_data_for_plotting$analysis_group)

#3.15 Set the renamed HC group values as the reference level
dimensional_analysis_data_for_plotting$analysis_group <- relevel(dimensional_analysis_data_for_plotting$analysis_group, ref = "HC")


## CBCL Data Analysis ## 

#1. Test whether symptom level between groups is associated with connectivity 
#1.11 Create a columns range for the CBCL columns
cbcl_col_range <- 16:19
connectivity_col_range <- 11:15

#1.12 Create an empty dataframe to store analysis values
cbcl_analysis_symptom_group_int_results_df <- data.frame(
  cbcl_column_name = character(),
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

#1.13 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  
  #1.131 Get the CBCL metric column name
  cbcl_col_name <- colnames(dimensional_analysis_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #1.132 Get the connectivity metric column name
    connectivity_col_name <- colnames(dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    
    #1.133 Run the linear regression model
    cbcl_analysis_symptom_group_int_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, "*analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = dimensional_analysis_data)
    
    #1.134 Run the ANCOVA Model
    cbcl_analysis_symptom_group_int_ANCOVA_analysis <- car::Anova(cbcl_analysis_symptom_group_int_mlm_analysis, type = "III", test.statistic = "F")
    
    #1.135 Loop through fixed effects in the model
    for (iv_name in row.names(cbcl_analysis_symptom_group_int_ANCOVA_analysis)) {
      summary_lm <- summary(cbcl_analysis_symptom_group_int_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #1.1351 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          cbcl_analysis_symptom_group_int_results_df <- rbind(
            cbcl_analysis_symptom_group_int_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
              column_name = connectivity_col_name,
              IV = iv_name,
              estimate = estimate,
              std_error = std_error,
              t_value = t_value,
              f_value = NA,
              df = NA,
              residual_df = NA,
              p_value = p_value))
        }
      } else {
        
        #1.1352 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- cbcl_analysis_symptom_group_int_ANCOVA_analysis[iv_name, "F"]
          df <- cbcl_analysis_symptom_group_int_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- cbcl_analysis_symptom_group_int_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- cbcl_analysis_symptom_group_int_ANCOVA_analysis[iv_name, "Pr(>F)"]
          cbcl_analysis_symptom_group_int_results_df <- rbind(
            cbcl_analysis_symptom_group_int_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
              column_name = connectivity_col_name,
              IV = iv_name,
              estimate = NA,
              std_error = NA,
              t_value = NA,
              f_value = f_value,
              df = df,
              residual_df = residual_df,
              p_value = p_value))
        }
      }
    }
  }
}

#1.21 Subset the data based on the relevant IV variable strings for FDR correction of p values
#1.211 Subset the data for interaction terms containing ":analysis_group" in the IV variable.
cbcl_analysis_symptom_group_int_p_adjust_subset <- subset(cbcl_analysis_symptom_group_int_results_df, grepl(":analysis_group", IV))

#1.212 Subset the data for main effects containing "cbcl_" in the IV variable.
cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset <- subset(cbcl_analysis_symptom_group_int_results_df, grepl("cbcl_", IV))

#1.213 Further filter the CBCL subset to include only rows where the "f_value" is NA.
cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset <- cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset[is.na(cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset$f_value), ]

#1.221 Create a new column conducting an FDR (p adjustment) on the derived p values
#1.2211 Interaction term p values
cbcl_analysis_symptom_group_int_p_adjust_subset$p_adjusted <- p.adjust(cbcl_analysis_symptom_group_int_p_adjust_subset$p_value, method = "fdr")

#1.2212 Main effect p values
cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset$p_adjusted <- p.adjust(cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset$p_value, method = "fdr")

#1.222 Merge the FDR corrected p-values together
cbcl_analysis_symptom_group_int_p_adjust_subset <- full_join(cbcl_analysis_symptom_group_int_p_adjust_subset, cbcl_analysis_symptom_group_int_cbcl_p_adjust_subset)

#1.23 Create a new dataframe wherein only significant results are stored
cbcl_analysis_symptom_group_int_significant_results <- subset(cbcl_analysis_symptom_group_int_p_adjust_subset, p_adjusted <= 0.05)

#1.24 Join the FDR corrected p-values with the rest of the model results 
cbcl_analysis_symptom_group_int_merged_adjusted_p_values <- left_join(cbcl_analysis_symptom_group_int_results_df, cbcl_analysis_symptom_group_int_p_adjust_subset)

#1.25 Remove all numbers from the strings in the column IV (for later merging purposes)
cbcl_analysis_symptom_group_int_merged_adjusted_p_values$IV <- gsub("\\d+", "", cbcl_analysis_symptom_group_int_merged_adjusted_p_values$IV)

#1.261 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
cbcl_analysis_symptom_group_int_merged_adjusted_p_values <- cbcl_analysis_symptom_group_int_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "analysis_group" ~ "Diagnosis Group",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint",
      grepl("^cbcl_.*:analysis_group", IV) ~ "CBCL Score*Diagnosis Group Interaction",
      TRUE ~ IV))

#1.262 Remove rows with Intercept values
cbcl_analysis_symptom_group_int_merged_adjusted_p_values <- cbcl_analysis_symptom_group_int_merged_adjusted_p_values[cbcl_analysis_symptom_group_int_merged_adjusted_p_values$IV != "(Intercept)", ]

#1.27 Pivot the full model results to be in wide format
cbcl_analysis_symptom_group_int_merged_adjusted_p_values_pivoted <- cbcl_analysis_symptom_group_int_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#1.28 Reorder columns to match the order in the paper
cbcl_analysis_symptom_group_int_merged_adjusted_p_values_pivoted <-
  cbcl_analysis_symptom_group_int_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(
      cbcl_column_name,
      column_name,
      `CBCL Score*Diagnosis Group Interaction_f_value`,
      `CBCL Score*Diagnosis Group Interaction_df`,
      `CBCL Score*Diagnosis Group Interaction_residual_df`,
      `CBCL Score*Diagnosis Group Interaction_p_value`,
      `CBCL Score*Diagnosis Group Interaction_p_adjusted`,
      `CBCL Score_estimate`,
      `CBCL Score_t_value`,
      `CBCL Score_std_error`,
      `CBCL Score_p_value`,
      `CBCL Score_p_adjusted`,
      `Diagnosis Group_f_value`,
      `Diagnosis Group_df`,
      `Diagnosis Group_residual_df`,
      `Diagnosis Group_p_value`,
      `FD Motion_estimate`,
      `FD Motion_t_value`,
      `FD Motion_std_error`,
      `FD Motion_p_value`,
      `Sex_f_value`,
      `Sex_df`,
      `Sex_residual_df`,
      `Sex_p_value`,
      `Timepoint_f_value`,
      `Timepoint_df`,
      `Timepoint_residual_df`,
      `Timepoint_p_value`))

#1.29 Write the full model results as a csv file
write.csv(cbcl_analysis_symptom_group_int_merged_adjusted_p_values_pivoted, "./results/cbcl_merged_intx_analysis_results.csv", row.names = FALSE)


#1.3 Plot the FC ~ CBCL*DX Group Data
#1.31 Define the Merged group CBCL scores and connectivity metrics
#1.1311 Define a vector containing the CBCL score variable names for the merged group analysis
Merged_cbcl_scores <- c("cbcl_scr_syn_anxdep_t", "cbcl_scr_syn_external_t", "cbcl_scr_syn_internal_t", "cbcl_scr_dsm5_anxdisord_t")

#1.1312 Define a vector containing descriptive labels corresponding to the CBCL score variables
Merged_cbcl_labels <- c("Anxious Depressed", "Externalizing", "Internalizing", "DSM-5 Anxiety")

#1.1313 Define a vector containing the connectivity metric variable names for the merged group analysis
Merged_connectivity_metrics <- c("rsfmri_cor_ngd_cerc_scs_cdelh", "rsfmri_cor_ngd_cerc_scs_aglh", "rsfmri_c_ngd_vta_ngd_vta", "rsfmri_cor_ngd_df_scs_ptlh", "rsfmri_cor_ngd_sa_scs_ptlh")

#1.1314 Define a vector containing descriptive labels corresponding to the connectivity metrics
Merged_connectivity_labels <- c("CON - Left Caudate", "CON - Left Amygdala", "Within-VAN", "DMN - Left Putamen", "SN - Left Putamen")

#1.32 Create an empty list to store the plots
Merged_plots_list <- list()

#1.33 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(Merged_cbcl_scores)) {
  for (j in seq_along(Merged_connectivity_metrics)) {
    cbcl_score <- Merged_cbcl_scores[i]
    cbcl_label <- Merged_cbcl_labels[i]
    connectivity_metric <- Merged_connectivity_metrics[j]
    connectivity_label <- Merged_connectivity_labels[j]

    #1.332 Extract the beta estimate from the dataframe
    beta_value <- cbcl_analysis_symptom_group_int_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric, cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)
    
    #1.333 Create each merged scatterplot
    Merged_scatterplot <- ggplot(dimensional_analysis_data_for_plotting, 
                                 aes_string(x = cbcl_score, y = connectivity_metric, color = "analysis_group")) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      scale_color_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8, margin = margin(0.01, unit = "cm"), angle = 0),
            axis.text.y = element_text(size = 8, margin = margin(0.01, unit = "cm"), angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))

    #1.334 Annotate each plot with beta value
    Merged_scatterplot <- Merged_scatterplot +
      annotate("text",
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,     # Adjust horizontal alignment
               vjust = 1.1)   # Adjust vertical alignment
    
    #1.335 Store each plot in the list of plots to facet
    Merged_plots_list[[length(Merged_plots_list) + 1]] <- Merged_scatterplot
  }
}

#1.34 Combine all plots into a facet plot
Merged_facet_plot <- do.call(grid.arrange, c(Merged_plots_list, nrow = 5))

#1.35 Create a single plot with a legend from one of the plots
Merged_facet_legend_plot <- ggplot(dimensional_analysis_data_for_plotting, 
                                   aes_string(x = Merged_cbcl_scores[1], 
                                              y = Merged_connectivity_metrics[1], 
                                              color = "analysis_group")) +
  geom_point() +
  scale_color_manual(values = c("HC" = "#218380", "GAD" = "#D81159")) +
  theme(legend.position = "bottom") +
  labs(color = "Diagnostic Group")

#1.36 Extract the legend from the single plot
Merged_facet_legend <- cowplot::get_plot_component(Merged_facet_legend_plot, 'guide-box-bottom', return_all = TRUE)

#1.371 Combine the facet plot and the extracted legend
Merged_facet_plot_with_legend <- cowplot::plot_grid(Merged_facet_plot, Merged_facet_legend, ncol = 1, rel_heights = c(1, 0.1))

#1.372 Print the final plot
print(Merged_facet_plot_with_legend)


#1.38 Save the final facet plot including the legend
ggsave("./results/cbcl_merged_intx_analysis_facet_plot_with_legend.pdf", Merged_facet_plot_with_legend, bg = "white", width = 8.5, height = 11, units = "in", device = "pdf")


#2. Determine if there are any within GAD group differences in connectivity based on CBCL symptom dimension expression 
#2.1 Create a columns range for the CBCL and connectivity columns
cbcl_col_range <- 16:19
connectivity_col_range <- 11:15

#2.2 Create an empty dataframe to store analysis values
cbcl_analysis_GAD_results_df <- data.frame(cbcl_column_name = character(),
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

#2.23 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  
  #2.231 Get the CBCL metric column name
  cbcl_col_name <- colnames(gad_group_dimensional_analysis_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #2.232 Get the connectivity metric column name
    connectivity_col_name <- colnames(gad_group_dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    
    #2.233 Run the linear regression model
    cbcl_analysis_GAD_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, " + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = gad_group_dimensional_analysis_data)
    
    #2.234 Run the ANCOVA Model
    cbcl_analysis_GAD_ANCOVA_analysis <- car::Anova(cbcl_analysis_GAD_mlm_analysis, type = "II", test.statistic = "F")
    
    #2.235 Loop through fixed effects in the model
    for (iv_name in row.names(cbcl_analysis_GAD_ANCOVA_analysis)) {
      summary_lm <- summary(cbcl_analysis_GAD_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #2.2351 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          cbcl_analysis_GAD_results_df <- rbind(
            cbcl_analysis_GAD_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        
        #2.2352 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- cbcl_analysis_GAD_ANCOVA_analysis[iv_name, "F"]
          df <- cbcl_analysis_GAD_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- cbcl_analysis_GAD_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- cbcl_analysis_GAD_ANCOVA_analysis[iv_name, "Pr(>F)"]
          cbcl_analysis_GAD_results_df <- rbind(
            cbcl_analysis_GAD_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}

#2.31 Subset the data based on the relevant IV variable strings for FDR correction of p values
cbcl_analysis_GAD_cbcl_p_adjust_subset <- subset(cbcl_analysis_GAD_results_df, grepl("cbcl_", IV))

#2.32 Create a new column conducting an FDR (p adjustment) on the derived p values
cbcl_analysis_GAD_cbcl_p_adjust_subset$p_adjusted <- p.adjust(cbcl_analysis_GAD_cbcl_p_adjust_subset$p_value, method = "fdr")

#2.33 Create a new dataframe where only significant results are stored
cbcl_analysis_GAD_significant_results <- subset(cbcl_analysis_GAD_cbcl_p_adjust_subset, p_adjusted <= 0.05)

#2.34 Join the FDR corrected p-values with the rest of the model results 
cbcl_analysis_GAD_merged_adjusted_p_values <- left_join(cbcl_analysis_GAD_results_df, cbcl_analysis_GAD_cbcl_p_adjust_subset)

#2.35 Remove all numbers from the strings in the column IV (for later merging purposes)
cbcl_analysis_GAD_merged_adjusted_p_values$IV <- gsub("\\d+", "", cbcl_analysis_GAD_merged_adjusted_p_values$IV)

#2.361 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
cbcl_analysis_GAD_merged_adjusted_p_values <- cbcl_analysis_GAD_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))

#2.362 Remove rows with Intercept values
cbcl_analysis_GAD_merged_adjusted_p_values <- cbcl_analysis_GAD_merged_adjusted_p_values[cbcl_analysis_GAD_merged_adjusted_p_values$IV != "(Intercept)", ]

#2.37 Pivot the full model results to be wider
cbcl_analysis_GAD_merged_adjusted_p_values_pivoted <- cbcl_analysis_GAD_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#2.38 Reorder columns to match the order in the paper results
cbcl_analysis_GAD_merged_adjusted_p_values_pivoted <- cbcl_analysis_GAD_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `CBCL Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))

#2.4 Write the full model results as a csv file
write.csv(cbcl_analysis_GAD_merged_adjusted_p_values_pivoted, "./results/cbcl_gad_group_analysis_results.csv", row.names = FALSE)

#2.5 Plot the GAD group dimensional (CBCL) analysis results
#2.511 Convert the beta values of interest to numeric for in plot formatting
cbcl_analysis_GAD_merged_adjusted_p_values_pivoted$`CBCL Score_estimate` <- as.numeric(cbcl_analysis_GAD_merged_adjusted_p_values_pivoted$`CBCL Score_estimate`)

#2.512 Define the GAD group CBCL scores and connectivity metrics
GAD_cbcl_scores <-
  c("cbcl_scr_syn_anxdep_t",
    "cbcl_scr_syn_external_t",
    "cbcl_scr_syn_internal_t", 
    "cbcl_scr_dsm5_anxdisord_t")
GAD_cbcl_labels <-
  c("Anxious Depressed",
    "Externalizing",
    "Internalizing",
    "DSM-5 Anxiety")
GAD_connectivity_metrics <-
  c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh")
GAD_connectivity_labels <-
  c("CON - Left Caudate",
    "CON - Left Amygdala",
    "Within-VAN",
    "DMN - Left Putamen",
    "SN - Left Putamen")

#2.52 Create an empty list to store the plots
GAD_plots_list <- list()

#2.53 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(GAD_cbcl_scores)) {
  for (j in seq_along(GAD_connectivity_metrics)) {
    cbcl_score <- GAD_cbcl_scores[i]
    cbcl_label <- GAD_cbcl_labels[i]
    connectivity_metric <- GAD_connectivity_metrics[j]
    connectivity_label <- GAD_connectivity_labels[j]
    
    #2.531 Extract the beta estimate from the dataframe
    beta_value <-
      cbcl_analysis_GAD_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric & 
             cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)

    #2.532 Create GAD group scatterplot
    GAD_scatterplot <-
      ggplot(gad_group_dimensional_analysis_data, aes_string(x = cbcl_score, y = connectivity_metric)) +
      geom_point(color = "#D81159") +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black",
                                     size = 1,
                                     linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8,
                                       margin = margin(0.01, unit = "cm"),
                                       angle = 0),
            axis.text.y = element_text(
              size = 8,
              margin = margin(0.01, unit = "cm"),
              angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))
    
    #2.533 Annotate with beta value
    GAD_scatterplot <- GAD_scatterplot +
      annotate("text",
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,   # Adjust horizontal alignment
               vjust = 1.1)
    
    #2.534 Store the plot in the list
    GAD_plots_list[[length(GAD_plots_list) + 1]] <-
      GAD_scatterplot
  }
}

#2.54 Combine all plots into a facet plot
GAD_facet_plot <-
  do.call(grid.arrange, c(GAD_plots_list, nrow = 5))

#2.55 Save the GAD group facet plot
ggsave("./results/cbcl_gad_group_analysis_facet_plot.pdf", GAD_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300, device = "pdf")


#3. Determine if there are any within Control group differences in connectivity based on CBCL symptom dimension expression 
#3.1 Create an empty dataframe to store analysis values
cbcl_analysis_control_results_df <-
  data.frame(cbcl_column_name = character(),
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

#3.2 Run the lm for each combination of CBCL metric and connectivity metric
for (cbcl_col_num in cbcl_col_range) {
  
  #3.21 Get the CBCL metric column name
  cbcl_col_name <- colnames(hc_group_dimensional_analysis_data)[cbcl_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #3.22 Get the connectivity metric column name
    connectivity_col_name <- colnames(hc_group_dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", cbcl_col_name))
    
    #3.23 Run the linear regression model
    cbcl_analysis_control_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", cbcl_col_name, " + rsfmri_c_ngd_meanmotion + sex + eventname + (1|site_name) + (1|family_id)"), na.action = na.omit, data = hc_group_dimensional_analysis_data)
    
    #3.24 Run the ANCOVA Model
    cbcl_analysis_control_ANCOVA_analysis <- car::Anova(cbcl_analysis_control_mlm_analysis, type = "II", test.statistic = "F")
    
    #3.25 Loop through fixed effects in the model
    for (iv_name in row.names(cbcl_analysis_control_ANCOVA_analysis)) {
      summary_lm <- summary(cbcl_analysis_control_mlm_analysis)
      if (iv_name == cbcl_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #3.251 For continuous variables (cbcl_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          cbcl_analysis_control_results_df <- rbind(
            cbcl_analysis_control_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        
        #3.252 For categorical variables (analysis_group, analysis_group*cbcl_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- cbcl_analysis_control_ANCOVA_analysis[iv_name, "F"]
          df <- cbcl_analysis_control_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- cbcl_analysis_control_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- cbcl_analysis_control_ANCOVA_analysis[iv_name, "Pr(>F)"]
          cbcl_analysis_control_results_df <- rbind(
            cbcl_analysis_control_results_df,
            data.frame(cbcl_column_name = cbcl_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}

#3.31 Subset the data based on the relevant IV variable strings for FDR correction of p values
cbcl_analysis_control_cbcl_p_adjust_subset <- subset(cbcl_analysis_control_results_df, grepl("cbcl_", IV))

#3.32 Create a new column conducting an FDR (p adjustment) on the derived p values
cbcl_analysis_control_cbcl_p_adjust_subset$p_adjusted <- p.adjust(cbcl_analysis_control_cbcl_p_adjust_subset$p_value, method = "fdr")

#3.33 Create a new dataframe wherein only significant results are stored
cbcl_analysis_control_significant_results <- subset(cbcl_analysis_control_cbcl_p_adjust_subset, p_adjusted <= 0.05)

#3.34 Join the FDR corrected p-values with the rest of the model results 
cbcl_analysis_control_merged_adjusted_p_values <- left_join(cbcl_analysis_control_results_df, cbcl_analysis_control_cbcl_p_adjust_subset)

#3.35 Remove all numbers from the strings in the column IV (for later merging purposes)
cbcl_analysis_control_merged_adjusted_p_values$IV <- gsub("\\d+", "", cbcl_analysis_control_merged_adjusted_p_values$IV)

#3.361 Rename IV values to match paper's format (e.g., CBCL Score, Diagnosis Group, FD Motion, Sex, Timepoint)
cbcl_analysis_control_merged_adjusted_p_values <- cbcl_analysis_control_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^cbcl_", IV) & !grepl("analysis_group", IV) ~ "CBCL Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))

#3.362 Remove rows with Intercept values
cbcl_analysis_control_merged_adjusted_p_values <- cbcl_analysis_control_merged_adjusted_p_values[cbcl_analysis_control_merged_adjusted_p_values$IV != "(Intercept)", ]

#3.37 Pivot the full model results to be wider
cbcl_analysis_control_merged_adjusted_p_values_pivoted <- cbcl_analysis_control_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(cbcl_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#3.38 Reorder columns to match the order in your paper
cbcl_analysis_control_merged_adjusted_p_values_pivoted <- cbcl_analysis_control_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(cbcl_column_name, column_name, `CBCL Score_estimate`, `CBCL Score_t_value`, `CBCL Score_std_error`, `CBCL Score_p_value`, `CBCL Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`, `Timepoint_f_value`, `Timepoint_df`, `Timepoint_residual_df`, `Timepoint_p_value`))

#3.354 Write the full model results as a csv file
write.csv(cbcl_analysis_control_merged_adjusted_p_values_pivoted, "./results/cbcl_hc_group_analysis_results.csv", row.names = FALSE)


#3.4 Plot the results of the hc group cbcl analysis
#3.41 Convert the beta values of interest to numeric for in plot formatting
cbcl_analysis_control_merged_adjusted_p_values_pivoted$`CBCL Score_estimate` <- as.numeric(cbcl_analysis_control_merged_adjusted_p_values_pivoted$`CBCL Score_estimate`)

#3.51 Define the control group CBCL scores and connectivity metrics
control_cbcl_scores <-
  c("cbcl_scr_syn_anxdep_t",
    "cbcl_scr_syn_external_t",
    "cbcl_scr_syn_internal_t",
    "cbcl_scr_dsm5_anxdisord_t")
control_cbcl_labels <-
  c("Anxious Depressed",
    "Externalizing",
    "Internalizing",
    "DSM-5 Anxiety")
control_connectivity_metrics <-
  c("rsfmri_cor_ngd_cerc_scs_cdelh",
    "rsfmri_cor_ngd_cerc_scs_aglh",
    "rsfmri_c_ngd_vta_ngd_vta",
    "rsfmri_cor_ngd_df_scs_ptlh",
    "rsfmri_cor_ngd_sa_scs_ptlh")
control_connectivity_labels <-
  c("CON - Left Caudate",
    "CON - Left Amygdala",
    "Within-VAN",
    "DMN - Left Putamen",
    "SN - Left Putamen")

#3.52 Create an empty list to store the plots
control_plots_list <- list()

#3.53 Loop through each combination of CBCL score and connectivity metric
for (i in seq_along(control_cbcl_scores)) {
  for (j in seq_along(control_connectivity_metrics)) {
    cbcl_score <- control_cbcl_scores[i]
    cbcl_label <- control_cbcl_labels[i]
    connectivity_metric <- control_connectivity_metrics[j]
    connectivity_label <- control_connectivity_labels[j]
    
    #3.531 Extract the beta estimate from the dataframe
    beta_value <-
      cbcl_analysis_control_merged_adjusted_p_values_pivoted %>%
      filter(column_name == connectivity_metric,
             cbcl_column_name == cbcl_score) %>%
      pull(`CBCL Score_estimate`)

    #3.532 Create control_scatterplot
    control_scatterplot <-
      ggplot(hc_group_dimensional_analysis_data, aes_string(x = cbcl_score, y = connectivity_metric)) +
      geom_point(color = "#218380") +
      geom_smooth(method = "lm", se = TRUE) +
      coord_cartesian(ylim = c(-0.6, 0.6)) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black",
                                     size = 1,
                                     linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x = element_text(size = 8,
                                       margin = margin(0.01, unit = "cm"),
                                       angle = 0),
            axis.text.y = element_text(
              size = 8,
              margin = margin(0.01, unit = "cm"),
              angle = 0),
            axis.title.x = element_text(size = 8, margin = margin(t = 8)),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none") +
      labs(x = paste("CBCL", cbcl_label, "T Score"),
           y = paste(connectivity_label))
    
    #3.533 Annotate with beta value
    control_scatterplot <- control_scatterplot +
      annotate("text",
               x = Inf,       # Top right corner for x-axis
               y = Inf,       # Top right corner for y-axis
               label = paste("b =", round(as.numeric(beta_value), 4)),
               size = 3,      # Adjust size as needed
               hjust = 1,   # Adjust horizontal alignment
               vjust = 1.1)

    #3.534 Store the plot in the list
    control_plots_list[[length(control_plots_list) + 1]] <-
      control_scatterplot
  }
}

#3.54 Combine all plots into a facet plot
control_facet_plot <- do.call(grid.arrange, c(control_plots_list, nrow = 5))

#3.55 Save the plot
ggsave("./results/cbcl_hc_group_analysis_facet_plot.pdf", control_facet_plot, bg = "white", width = 8.5, height = 11, units = "in", dpi = 300, device = "pdf")


## BPM Data Analysis ##

#1. Test whether symptom is associated with connectivity in the whole group
#1.1 Create a column range for the bpm and connectivity columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14

#1.2 Create an empty dataframe to store analysis values
symptom_group_int_bpm_con_results_df <- data.frame(IV = character(), column_name = character(), bpm_col_name = character(), f_value = numeric(), df = numeric(), residual_df = numeric(), p_value = numeric())

#1.3 Run the lm for each combination of bpm metric and connectivity metric
symptom_group_int_bpm_con_results_df <- data.frame(
  bpm_column_name = character(),
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

#1.3 Run the lm for each combination of bpm metric and connectivity metric
for (bpm_col_num in bpm_col_range) {
  
  #1.31 Get the bpm metric column name
  bpm_col_name <- colnames(dimensional_analysis_data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #1.32 Get the connectivity metric column name
    connectivity_col_name <- colnames(dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    
    #1.33 Run the linear regression model
    symptom_group_int_bpm_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, "*analysis_group + rsfmri_c_ngd_meanmotion + sex + (1|site_name) + (1|family_id)"), na.action = na.omit, data = dimensional_analysis_data)
    
    #1.34 Run the ANCOVA Model
    symptom_group_int_bpm_con_ANCOVA_analysis <- car::Anova(symptom_group_int_bpm_con_mlm_analysis, type = "III", test.statistic = "F")
    
    #1.35 Loop through fixed effects in the model
    for (iv_name in row.names(symptom_group_int_bpm_con_ANCOVA_analysis)) {
      summary_lm <- summary(symptom_group_int_bpm_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #1.351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          symptom_group_int_bpm_con_results_df <- rbind(
            symptom_group_int_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        
        #1.352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- symptom_group_int_bpm_con_ANCOVA_analysis[iv_name, "F"]
          df <- symptom_group_int_bpm_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- symptom_group_int_bpm_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- symptom_group_int_bpm_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          symptom_group_int_bpm_con_results_df <- rbind(
            symptom_group_int_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}

#1.4 Subset the data based on the relevant IV variable strings for FDR correction of p values
symptom_group_int_bpm_con_p_adjust_subset <- subset(symptom_group_int_bpm_con_results_df, grepl(":analysis_group", IV))
symptom_group_int_bpm_con_bpm_p_adjust_subset <- subset(symptom_group_int_bpm_con_results_df, grepl("bpm_", IV))
symptom_group_int_bpm_con_bpm_p_adjust_subset <- symptom_group_int_bpm_con_bpm_p_adjust_subset[is.na(symptom_group_int_bpm_con_bpm_p_adjust_subset$f_value), ]

#1.51 Create a new column conducting an FDR (p adjustment) on the derived p values
symptom_group_int_bpm_con_p_adjust_subset$p_adjusted <- p.adjust(symptom_group_int_bpm_con_p_adjust_subset$p_value, method = "fdr")
symptom_group_int_bpm_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(symptom_group_int_bpm_con_bpm_p_adjust_subset$p_value, method = "fdr")

#1.52 Merge the FDR corrected p-values together
symptom_group_int_bpm_con_p_adjust_subset <- full_join(symptom_group_int_bpm_con_p_adjust_subset, symptom_group_int_bpm_con_bpm_p_adjust_subset)

#1.6 Create a new dataframe where only significant results are stored
symptom_group_int_bpm_con_significant_results <- subset(symptom_group_int_bpm_con_p_adjust_subset, p_adjusted <= 0.05)

#1.7 Join the FDR corrected p-values with the rest of the model results 
symptom_group_int_bpm_con_merged_adjusted_p_values <- left_join(symptom_group_int_bpm_con_results_df, symptom_group_int_bpm_con_p_adjust_subset)

#1.81 Remove all numbers from the strings in the column IV (for later merging purposes)
symptom_group_int_bpm_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", symptom_group_int_bpm_con_merged_adjusted_p_values$IV)

#1.82 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
symptom_group_int_bpm_con_merged_adjusted_p_values <- symptom_group_int_bpm_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "analysis_group" ~ "Diagnosis Group",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint",
      grepl("^bpm_.*:analysis_group", IV) ~ "bpm Score*Diagnosis Group Interaction",
      TRUE ~ IV))

#1.83 Remove rows with Intercept values
symptom_group_int_bpm_con_merged_adjusted_p_values <- symptom_group_int_bpm_con_merged_adjusted_p_values[symptom_group_int_bpm_con_merged_adjusted_p_values$IV != "(Intercept)", ]

#1.91 Pivot the full model results to be wide format
symptom_group_int_bpm_con_merged_adjusted_p_values_pivoted <- symptom_group_int_bpm_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#1.92 Reorder columns to match the order in the paper
symptom_group_int_bpm_con_merged_adjusted_p_values_pivoted <- symptom_group_int_bpm_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score*Diagnosis Group Interaction_f_value`, `bpm Score*Diagnosis Group Interaction_df`, `bpm Score*Diagnosis Group Interaction_residual_df`, `bpm Score*Diagnosis Group Interaction_p_value`, `bpm Score*Diagnosis Group Interaction_p_adjusted`, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `Diagnosis Group_f_value`, `Diagnosis Group_df`, `Diagnosis Group_residual_df`, `Diagnosis Group_p_value`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))

#1.93 Write the full model results as a csv file
write.csv(symptom_group_int_bpm_con_merged_adjusted_p_values_pivoted, "./results/bpm_merged_group_results.csv", row.names = FALSE)


#2. Test for differences in predicting connectivity based on bpm symptoms in the control group
#2.1 Create a columns range for the bpm columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14

#2.2 Create an empty dataframe to store analysis values
control_group_bpm_con_results_df <- data.frame(bpm_column_name = character(),
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

#2.3 Run the lm for each combination of bpm metric and connectivity metric
for (bpm_col_num in bpm_col_range) {
  
  #2.31 Get the bpm metric column name
  bpm_col_name <- colnames(hc_group_dimensional_analysis_data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #2.32 Get the connectivity metric column name
    connectivity_col_name <- colnames(hc_group_dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    
    #2.33 Run the linear regression model
    control_group_bpm_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, " + rsfmri_c_ngd_meanmotion + sex + (1|site_name)"), na.action = na.omit, data = hc_group_dimensional_analysis_data)
    
    #2.34 Run the ANCOVA Model
    control_group_bpm_con_ANCOVA_analysis <- car::Anova(control_group_bpm_con_mlm_analysis, type = "II", test.statistic = "F")
    
    #2.35 Loop through fixed effects in the model
    for (iv_name in row.names(control_group_bpm_con_ANCOVA_analysis)) {
      summary_lm <- summary(control_group_bpm_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #2.351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          control_group_bpm_con_results_df <- rbind(
            control_group_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        
        #2.352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- control_group_bpm_con_ANCOVA_analysis[iv_name, "F"]
          df <- control_group_bpm_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- control_group_bpm_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- control_group_bpm_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          control_group_bpm_con_results_df <- rbind(
            control_group_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}

#2.4 Subset the data based on the relevant IV variable strings for FDR correction of p values
control_group_bpm_con_bpm_p_adjust_subset <- subset(control_group_bpm_con_results_df, grepl("bpm_", IV))

#2.5 Create a new column conducting an FDR (p adjustment) on the derived p values
control_group_bpm_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(control_group_bpm_con_bpm_p_adjust_subset$p_value, method = "fdr")

#2.6 Create a new dataframe where only significant results are stored
control_group_bpm_con_significant_results <- subset(control_group_bpm_con_bpm_p_adjust_subset, p_adjusted <= 0.05)

#2.7 Join the FDR corrected p-values with the rest of the model results 
control_group_bpm_con_merged_adjusted_p_values <- left_join(control_group_bpm_con_results_df, control_group_bpm_con_bpm_p_adjust_subset)

#2.8 Remove all numbers from the strings in the column IV (for later merging purposes)
control_group_bpm_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", control_group_bpm_con_merged_adjusted_p_values$IV)

#2.91 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
control_group_bpm_con_merged_adjusted_p_values <- control_group_bpm_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))

#2.92 Remove rows with Intercept values
control_group_bpm_con_merged_adjusted_p_values <- control_group_bpm_con_merged_adjusted_p_values[control_group_bpm_con_merged_adjusted_p_values$IV != "(Intercept)", ]

#2.101 Pivot the full model results to be wider
control_group_bpm_con_merged_adjusted_p_values_pivoted <- control_group_bpm_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#2.102 Reorder columns to match the order in your paper
control_group_bpm_con_merged_adjusted_p_values_pivoted <- control_group_bpm_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))

#2.11 Write the full model results as a csv file
write.csv(control_group_bpm_con_merged_adjusted_p_values_pivoted, "./results/bpm_control_group_results.csv", row.names = FALSE)


#3. Test for differences in predicting connectivity based on bpm symptoms in the GAD group
#3.1 Create a columns range for the bpm columns
bpm_col_range <- 20:21
connectivity_col_range <- 10:14

#3.2 Create an empty dataframe to store analysis values
GAD_group_bpm_con_results_df <- data.frame(bpm_column_name = character(),
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

#3.3 Run the lm for each combination of bpm metric and connectivity metric
for (bpm_col_num in bpm_col_range) {
  
  #3.31 Get the bpm metric column name
  bpm_col_name <- colnames(gad_group_dimensional_analysis_data)[bpm_col_num]
  for (connectivity_col_num in connectivity_col_range) {
    
    #3.32 Get the connectivity metric column name
    connectivity_col_name <- colnames(gad_group_dimensional_analysis_data)[connectivity_col_num]
    print(paste("Modeling", connectivity_col_name, "~", bpm_col_name))
    
    #3.33 Run the linear regression model
    GAD_group_bpm_con_mlm_analysis <- lmerTest::lmer(formula = paste0(connectivity_col_name, " ~ ", bpm_col_name, " + rsfmri_c_ngd_meanmotion + sex + (1|site_name)"), na.action = na.omit, data = gad_group_dimensional_analysis_data)
    
    #3.34 Run the ANCOVA Model
    GAD_group_bpm_con_ANCOVA_analysis <- car::Anova(GAD_group_bpm_con_mlm_analysis, type = "II", test.statistic = "F")
    
    #3.35 Loop through fixed effects in the model
    for (iv_name in row.names(GAD_group_bpm_con_ANCOVA_analysis)) {
      summary_lm <- summary(GAD_group_bpm_con_mlm_analysis)
      if (iv_name == bpm_col_name || iv_name == "rsfmri_c_ngd_meanmotion") {
        
        #3.351 For continuous variables (bpm_col_name, rsfmri_c_ngd_meanmotion)
        matched_row <- grep(paste0("^", iv_name, "$"), rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          estimate <- summary_lm$coefficients[matched_row, "Estimate"]
          std_error <- summary_lm$coefficients[matched_row, "Std. Error"]
          t_value <- summary_lm$coefficients[matched_row, "t value"]
          p_value <- summary_lm$coefficients[matched_row, "Pr(>|t|)"]
          GAD_group_bpm_con_results_df <- rbind(
            GAD_group_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = estimate,
                       std_error = std_error,
                       t_value = t_value,
                       f_value = NA,
                       df = NA,
                       residual_df = NA,
                       p_value = p_value))
        }
      } else {
        
        #3.352 For categorical variables (analysis_group, analysis_group*bpm_col_name, sex, eventname)
        matched_row <- grep(iv_name, rownames(summary_lm$coefficients))
        print(paste("Matching", iv_name, "with row index", matched_row))
        if (length(matched_row) > 0) {
          f_value <- GAD_group_bpm_con_ANCOVA_analysis[iv_name, "F"]
          df <- GAD_group_bpm_con_ANCOVA_analysis[iv_name, "Df"]
          residual_df <- GAD_group_bpm_con_ANCOVA_analysis[iv_name, "Df.res"]
          p_value <- GAD_group_bpm_con_ANCOVA_analysis[iv_name, "Pr(>F)"]
          GAD_group_bpm_con_results_df <- rbind(
            GAD_group_bpm_con_results_df,
            data.frame(bpm_column_name = bpm_col_name, 
                       column_name = connectivity_col_name,
                       IV = iv_name,
                       estimate = NA,
                       std_error = NA,
                       t_value = NA,
                       f_value = f_value,
                       df = df,
                       residual_df = residual_df,
                       p_value = p_value))
        }
      }
    }
  }
}

#3.4 Subset the data based on the relevant IV variable strings for FDR correction of p values
GAD_group_bpm_con_bpm_p_adjust_subset <- subset(GAD_group_bpm_con_results_df, grepl("bpm_", IV))

#3.5 Create a new column conducting an FDR (p adjustment) on the derived p values
GAD_group_bpm_con_bpm_p_adjust_subset$p_adjusted <- p.adjust(GAD_group_bpm_con_bpm_p_adjust_subset$p_value, method = "fdr")

#3.6 Create a new dataframe where only significant results are stored
GAD_group_bpm_con_significant_results <- subset(GAD_group_bpm_con_bpm_p_adjust_subset, p_adjusted <= 0.05)

#3.7 Join the FDR corrected p-values with the rest of the model results 
GAD_group_bpm_con_merged_adjusted_p_values <- left_join(GAD_group_bpm_con_results_df, GAD_group_bpm_con_bpm_p_adjust_subset)

#3.8 Remove all numbers from the strings in the column IV (for later merging purposes)
GAD_group_bpm_con_merged_adjusted_p_values$IV <- gsub("\\d+", "", GAD_group_bpm_con_merged_adjusted_p_values$IV)

#3.91 Rename IV values to match paper's format (e.g., bpm Score, Diagnosis Group, FD Motion, Sex, Timepoint)
GAD_group_bpm_con_merged_adjusted_p_values <- GAD_group_bpm_con_merged_adjusted_p_values %>%
  mutate(
    IV = case_when(
      grepl("^bpm_", IV) & !grepl("analysis_group", IV) ~ "bpm Score",
      IV == "rsfmri_c_ngd_meanmotion" ~ "FD Motion",
      IV == "sex" ~ "Sex",
      IV == "eventname" ~ "Timepoint"))

#3.92 Remove rows with Intercept values
GAD_group_bpm_con_merged_adjusted_p_values <- GAD_group_bpm_con_merged_adjusted_p_values[GAD_group_bpm_con_merged_adjusted_p_values$IV != "(Intercept)", ]

#3.101 Pivot the full model results to be wider
GAD_group_bpm_con_merged_adjusted_p_values_pivoted <- GAD_group_bpm_con_merged_adjusted_p_values %>%
  pivot_wider(
    id_cols = c(bpm_column_name, column_name),
    names_from = IV,
    values_from = c(f_value, df, residual_df, p_value, p_adjusted, estimate, t_value, std_error),
    names_glue = "{IV}_{.value}")

#3.102 Reorder columns to match the order in your paper
GAD_group_bpm_con_merged_adjusted_p_values_pivoted <- GAD_group_bpm_con_merged_adjusted_p_values_pivoted %>%
  dplyr::select(c(bpm_column_name, column_name, `bpm Score_estimate`, `bpm Score_t_value`, `bpm Score_std_error`, `bpm Score_p_value`, `bpm Score_p_adjusted`, `FD Motion_estimate`, `FD Motion_t_value`, `FD Motion_std_error`, `FD Motion_p_value`, `Sex_f_value`, `Sex_df`, `Sex_residual_df`, `Sex_p_value`))

#3.11 Write the full model results as a csv file
write.csv(GAD_group_bpm_con_merged_adjusted_p_values_pivoted, "./results/bpm_GAD_group_results.csv", row.names = FALSE)