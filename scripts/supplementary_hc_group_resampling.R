## Set Up ##

# Load in necessary packages and configure environmental variables
library(dplyr)
options(digits = 8, scipen = 999) 

# Read in required data 
# Comorbidity + non-resampled HC group data
comorbidity_hc_group_not_resampled <- read.csv("./data_processed/supplemental_comorbidity_sample_hc_not_resampled.csv")


## Data Wrangling + Prep for Resampling ##

#1. Create a broader variable differentiating the HC and clinical (comorbidity) groups and clean the dataframe for resampling
#1.1 Create the broad clinical group variable to be used in resampling the HC group
broad_cleaned_comorbidity_hc_group_not_resampled <- comorbidity_hc_group_not_resampled %>% 
  mutate(broad_clinical_group = if_else(comorbidity_group == "HC", "HC", "Clinical_Group"))

#1.2 Retain only columns of interest to the resampling
broad_cleaned_comorbidity_hc_group_not_resampled <- broad_cleaned_comorbidity_hc_group_not_resampled %>% 
  dplyr::select(c(src_subject_id, eventname, site_name, comorbidity_group, broad_clinical_group))


#2. Create a dataframe specific to the HC and Clinical (comorbidity) groups
#2.1 Clinical (comorbidity) group
comorbidity_group <- broad_cleaned_comorbidity_hc_group_not_resampled %>% 
  filter(broad_clinical_group == "Clinical_Group")

#2.2 HC group
HC_group_not_resampled <- broad_cleaned_comorbidity_hc_group_not_resampled %>% 
  filter(broad_clinical_group == "HC")

#3. Determine if there are any subjects who have data at multiple sites, and if so, randomly sample the data point to keep for each subject
#3.1 Identify subjects who have different site_name values across different eventname's
multi_site_subjects <- HC_group_not_resampled %>%
  group_by(src_subject_id) %>%
  filter(n_distinct(site_name) > 1) %>%
  pull(src_subject_id) %>%
  unique()

#3.12 Randomly select one data point for each subject with assessment data at multiple sites, since all multi site subjects have both baseline and 2 year follow up data
#3.121 Initialize an empty dataframe 
multi_site_selected_data <- data.frame()

#3.122 Loop through the multi-site subjects and store their sampled data point
for (subject in multi_site_subjects) {
  
  #3.1221 Subset the data for the current subject
  subject_data <- subset(HC_group_not_resampled, src_subject_id == subject)
  
  #3.1222 Randomly select either the baseline or followup data point with a 50% probability of each
  if (runif(1) > 0.5) {
    multi_site_selected_data <-
      rbind(multi_site_selected_data, subject_data[subject_data$eventname == "baseline_year_1_arm_1",])
  } else {
    multi_site_selected_data <-
      rbind(multi_site_selected_data, subject_data[subject_data$eventname == "2_year_follow_up_y_arm_1",])
  }
}

#3.13 Replace the existing multi site subject data with the resampled data 
#3.131 Remove the current multi site subject data from the HC group 
HC_group_cleaned_not_resampled <- HC_group_not_resampled %>%
  filter(!(src_subject_id %in% multi_site_subjects))

#3.132 Add back in the selected multi site subject data
HC_group_cleaned_not_resampled <- full_join(HC_group_cleaned_not_resampled, multi_site_selected_data)

#4. Generate summary statistics outlining the distribution of comorbidity + HC subjects at each combination of site and eventname
#4.1 Create a dataframe containing the N and % of subjects at each combination of site and eventname, with the % of subjects at each timepoint within site rounded to 2 decimal places for purposes of maximizing the sample during matching
comorbidity_site_visit_distribution <- comorbidity_group %>% 
  group_by(site_name, eventname) %>%
  summarize(comorbidity_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(comorbidity_site_n = sum(comorbidity_site_visit_n),
         comorbidity_site_pct = round((comorbidity_site_visit_n / comorbidity_site_n), digits = 2))

#4.2 Likewise for the HC distribution, create a dataframe containing the N and % of subjects at each combination of site and eventname, with the % of subjects at each timepoint within site rounded to 2 decimal places for purposes of maximizing the sample during matching
HC_site_visit_distribution <- HC_group_cleaned_not_resampled %>% 
  group_by(site_name, eventname) %>%
  summarize(control_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(control_site_n = sum(control_site_visit_n),
         control_site_pct = round((control_site_visit_n / control_site_n), digits = 2))

#4.3 Merge the HC + comorbidity distributions at each site-visit combination in the comorbidity data
#4.31 Combine the comorbidity + HC distribution data
merged_HC_comorbidity_site_visit_distribution <- merge(comorbidity_site_visit_distribution, HC_site_visit_distribution, by = c("eventname", "site_name"), all.x =  TRUE)

#4.32 Add a new variable to the merged distribution stats that contains the maximum possible number of controls for each combination of site + eventname based on the comorbidity distribution and N controls
merged_HC_comorbidity_site_visit_distribution <- merged_HC_comorbidity_site_visit_distribution %>% 
  group_by(eventname, site_name) %>%
  mutate(max_controls = floor(comorbidity_site_pct * control_site_n))         


## HC Group Resampling ## 

#1. Create a function to recursively sample the healthy control data at every combination of site and time point
resample_healthy_controls <- function(site_id, HC_group_data, merged_HC_comorbidity_distribution, comorbidity_group_data) {
  
  #1.11 Resample the healthy control data for each specified site
  healthy_control_data <- HC_group_data %>%
    filter(site_name == site_id) %>%
    dplyr::select(src_subject_id, eventname, site_name)
  
  #1.12 Merge the healthy control data with the pertinent site x visit distribution
  merged_HC_comorbidity_site_visit_distribution <-
    merge(
      merged_HC_comorbidity_distribution,
      healthy_control_data, 
      by = c("eventname", "site_name"), all.y = TRUE)
  
  #1.13 Create the baseline and follow-up only subject dataframes
  #1.131 Baseline subjects only (no follow up data)
  healthy_control_baseline_only <- merged_HC_comorbidity_site_visit_distribution %>%
    group_by(src_subject_id) %>%
    filter(all(eventname == "baseline_year_1_arm_1") & n() == 1)
  
  #1.132 Followup subjects only (no baseline data)
  healthy_control_followup_only <- merged_HC_comorbidity_site_visit_distribution %>%
    group_by(src_subject_id) %>%
    filter(all(eventname == "2_year_follow_up_y_arm_1") & n() == 1)
  
  #1.133 Save the baseline + followup only subject counts
  baseline_only_n <- nrow(healthy_control_baseline_only)
  followup_only_n <- nrow(healthy_control_followup_only)
  
  #1.134 Save the subject IDs for baseline and follow-up only subjects
  baseline_control_subs_only <- healthy_control_baseline_only$src_subject_id
  followup_control_subs_only <- healthy_control_followup_only$src_subject_id
  
  #1.14 Store unique control IDs at the current site
  num_unique_control_subs <- length(unique(merged_HC_comorbidity_site_visit_distribution$src_subject_id))
  
  #1.21 Subset comorbidity subjects for this site
  comorbidity_subs <- comorbidity_group_data %>%
    filter(site_name == site_id)
  
  #1.22 Store the total number of comorbidity subjects at this site
  total_comorbidity_sub_n <- length(unique(comorbidity_subs$src_subject_id))
  
  #1.31 Further subset the comorbidity data into baseline and follow-up subjects
  baseline_comorbidity_data <- comorbidity_subs %>% filter(eventname == "baseline_year_1_arm_1")
  followup_comorbidity_data <- comorbidity_subs %>% filter(eventname == "2_year_follow_up_y_arm_1")
  
  #1.32 Calculate the comorbidity percentage distribution for baseline and follow up
  #1.321 Baseline comorbidity percentage distribution
  baseline_comorbidity_site_visit_percent <- merged_HC_comorbidity_site_visit_distribution %>% 
    filter(site_name == site_id & eventname == "baseline_year_1_arm_1") %>%
    pull(comorbidity_site_pct) %>% 
    unique()
  
  #1.322 Follow up comorbidity percentage distribution
  followup_comorbidity_site_visit_percent <- merged_HC_comorbidity_site_visit_distribution %>%
    filter(site_name == site_id & eventname == "2_year_follow_up_y_arm_1") %>%
    pull(comorbidity_site_pct) %>% 
    unique()
  
  #1.4 Handle the edge case present for some sites with only one time point of data collected
  if (is.null(followup_comorbidity_site_visit_percent) || 
      length(followup_comorbidity_site_visit_percent) == 0 || 
      all(is.na(followup_comorbidity_site_visit_percent))) {
    followup_comorbidity_site_visit_percent <- 0
  }
  
  #1.51 Create a list of IDs for + count of subjects that have data at both time points
  both_timepoint_sub_list <- healthy_control_data$src_subject_id[duplicated(healthy_control_data$src_subject_id)]
  both_timepoint_sub_list_n <- length(unique(both_timepoint_sub_list))
  
  #1.52 Create a dataframe with subjects who have data at both time points
  both_timepoint_subs <- healthy_control_data[healthy_control_data$src_subject_id %in% both_timepoint_sub_list, ]
  
  #1.61 Calculate the max N baseline and follow-up controls
  max_baseline_controls <- baseline_only_n + both_timepoint_sub_list_n
  max_followup_controls <- followup_only_n + both_timepoint_sub_list_n
  
  #1.62 Calculate the desired N of controls at both baseline and followup based on the comorbidity distribution
  baseline_controls_desired <- floor(num_unique_control_subs * baseline_comorbidity_site_visit_percent)
  followup_controls_desired <- floor(num_unique_control_subs * followup_comorbidity_site_visit_percent)
  
  #1.63 Check whether there are currently sufficient subjects to meet the max control N at both timepoints 
  sufficient_baseline_control_n <- max_baseline_controls > baseline_controls_desired
  sufficient_followup_control_n <- max_followup_controls > followup_controls_desired
  
  #1.64 Create an adjusted total control N at each timepoint if there are insufficient controls to sample the max possible number at each timepoint
  #1.641 In the event there are insufficient controls available at baseline 
  adjusted_total_control_n_if_insufficient_baseline <-
    floor(max_baseline_controls / baseline_comorbidity_site_visit_percent)
  
  #1.642 In the event there are insufficient controls available at followup 
  adjusted_total_control_n_if_insufficient_followup <-
    ifelse(
      followup_comorbidity_site_visit_percent > 0,
      floor(max_followup_controls / followup_comorbidity_site_visit_percent), 0)
  
  #1.65 Initialize a control sample counts dataframe
  updated_control_counts <- data.frame(num_baseline_controls = 0, num_followup_controls = 0)
  
  #1.66 Adjust the number of controls to sample based on relevant limitations
  #1.661 If not enough baseline controls, set to max and adjust follow-up accordingly
  if (!sufficient_baseline_control_n & sufficient_followup_control_n) {
    updated_control_counts$num_baseline_controls <-
      max_baseline_controls
    updated_control_counts$num_followup_controls <-
      floor(adjusted_total_control_n_if_insufficient_baseline * followup_comorbidity_site_visit_percent)
    
    #1.662 If not enough follow-up controls, set to max and adjust baseline accordingly
  } else if (sufficient_baseline_control_n & !sufficient_followup_control_n) {
    updated_control_counts$num_baseline_controls <-
      floor(adjusted_total_control_n_if_insufficient_followup * baseline_comorbidity_site_visit_percent)
    updated_control_counts$num_followup_controls <-
      max_followup_controls
    
    #1.663 If both baseline and follow-up control N's are sufficient, use the calculated desired counts  
  } else {
    updated_control_counts$num_baseline_controls <-
      baseline_controls_desired
    updated_control_counts$num_followup_controls <-
      followup_controls_desired
  }
  
  #1.71 Sample baseline controls
  baseline_control_sample <-
    healthy_control_data[healthy_control_data$src_subject_id %in% baseline_control_subs_only,]
  
  #1.711 Get the number of baseline controls to sample and the current related count
  num_baseline_controls_to_sample <-
    updated_control_counts$num_baseline_controls
  current_baseline_control_n <-
    length(unique(baseline_control_sample$src_subject_id))
  
  #1.712 Check if the current baseline controls match the desired count
  if (current_baseline_control_n == num_baseline_controls_to_sample) {
    message("Number of baseline controls matched expected number")
    
    #1.713 If more controls are present than needed, randomly sample the required number
  } else {
    if (current_baseline_control_n > num_baseline_controls_to_sample) {
      baseline_control_sample <-
        baseline_control_sample[sample(nrow(baseline_control_sample),
                                       num_baseline_controls_to_sample),]
      #1.714 If fewer controls are present than needed, add from the both timepoint list of subjects
    } else {
      baseline_control_sample_to_add <-
        healthy_control_data[healthy_control_data$src_subject_id %in% both_timepoint_sub_list &
                               healthy_control_data$eventname == "baseline_year_1_arm_1",]
      baseline_control_sample_to_add <-
        baseline_control_sample_to_add[sample(
          nrow(baseline_control_sample_to_add),
          num_baseline_controls_to_sample - current_baseline_control_n
        ),]
      baseline_control_sample <-
        rbind(baseline_control_sample, baseline_control_sample_to_add)
    }
  }
  
  #1.72 Sample followup controls
  followup_control_sample <-
    healthy_control_data[healthy_control_data$src_subject_id %in% followup_control_subs_only,]
  
  #1.721 Get the number of followup controls to sample and the current related count
  num_followup_controls_to_sample <-
    updated_control_counts$num_followup_controls
  current_followup_control_n <-
    length(unique(followup_control_sample$src_subject_id))
  subs_already_sampled <- baseline_control_sample$src_subject_id
  
  #1.722 Check if the current baseline controls match the desired count
  if (current_followup_control_n == num_followup_controls_to_sample) {
    message("Number of followup controls matched expected number")
    
    #1.723 If more controls are present than needed, randomly sample the required number
  } else {
    if (current_followup_control_n > num_followup_controls_to_sample) {
      followup_control_sample <-
        followup_control_sample[sample(nrow(followup_control_sample),
                                       num_followup_controls_to_sample),]
      
      #1.724 If fewer controls are present than needed, add from the both timepoint list of subjects
    } else {
      followup_control_sample_to_add <-
        healthy_control_data[healthy_control_data$src_subject_id %in% both_timepoint_sub_list &
                               !healthy_control_data$src_subject_id %in% subs_already_sampled &
                               healthy_control_data$eventname == "2_year_follow_up_y_arm_1",]
      followup_control_sample_to_add <-
        followup_control_sample_to_add[sample(
          nrow(followup_control_sample_to_add),
          num_followup_controls_to_sample - current_followup_control_n
        ),]
      followup_control_sample <-
        rbind(followup_control_sample, followup_control_sample_to_add)
    }
  }
  
  #1.73 Merge the baseline and follow-up resampled controls
  control_sample <- merge(baseline_control_sample, followup_control_sample, all = TRUE)
  
  #1.81 Verify unique subjects
  unique_control_sample <- length(unique(control_sample$src_subject_id)) == nrow(control_sample)
  
  #1.82 If the sampled subjects are unique, the function will exit accordingly 
  if (unique_control_sample) {
    print("Only unique controls sampled")
  } else {
    print("Duplicated controls present - troubleshoot sampling process")
  }
  
  #1.9 Return the sampled subjects as a dataframe to be stored in the list of subjects sampled at each site
  return(control_sample)
}

#2. Employ the resampling function to resample the healthy control group
#2.11 Get unique site names from HC_group_cleaned_not_resampled
unique_sites <- unique(HC_group_cleaned_not_resampled$site_name)

#2.12 Initialize an empty list to store samples for each site
healthy_control_samples_list <- list()

#2.2 Loop over each unique site and apply the resample function
for (site in unique_sites) {
  healthy_control_sample <- resample_healthy_controls(site, HC_group_cleaned_not_resampled, merged_HC_comorbidity_site_visit_distribution, comorbidity_group)
  
  #2.21 Store the result in the list with the site name as the key
  healthy_control_samples_list[[site]] <- healthy_control_sample
}

#2.3 Combine all samples into a single dataframe
combined_healthy_control_samples <- bind_rows(healthy_control_samples_list, .id = "site_name")

#2.41 Add back in the clinical + analysis group columns for the HC group
healthy_control_sample_resampled <- combined_healthy_control_samples %>% 
  mutate(GAD_Source = rep("None"),
         Social_Anxiety_Disorder = 0,
         Social_Anxiety_Disorder_Source = rep("None"),
         Separation_Anxiety_Disorder = 0,
         Separation_Anxiety_Disorder_Source = rep("None"),
         MDD = 0,
         MDD_Source = rep("None"),
         comorbidity_group = rep("HC"))

#2.42 Add back in the demographic information 
healthy_control_sample_resampled <- left_join(healthy_control_sample_resampled, merged_comorbidity_HC_sample_data_cleaned)

#3. Merge the comorbidity group with the newly sampled HC group
comorbidity_HC_sample <- full_join(comorbidity_group, healthy_control_sample_resampled)


## Verify Accuracy + Reliability of Resampling Procedure ##

#1. Determine whether the expected number/proportion of subjects were successfully sampled within the HC group at each combination of site + timepoint 
#1.1 Create a dataframe containing the distribution of sampled controls at each site + timepoint combination (both N and % distribution)
resampled_HC_site_visit_distribution <- healthy_control_sample_resampled %>% 
  group_by(site_name, eventname) %>%
  summarize(control_site_visit_n = n()) %>%
  group_by(site_name) %>%
  mutate(control_site_n = sum(control_site_visit_n),
         control_site_pct = round((control_site_visit_n / control_site_n), digits = 2))

#1.21 Create a dataframe comparing the distribution of comorbidity subjects and sampled controls at each site-event combination (both N and % distribution)
resampled_HC_site_visit_distribution_comparison <- merge(comorbidity_site_visit_distribution, resampled_HC_site_visit_distribution, by = c("site_name", "eventname"), all = TRUE)

#1.22 Add a column that tests whether the comorbidity distribution matches the sampled control subject distribution (allowing for a 0.01 difference due to rounding error by r's "round" function)
resampled_HC_site_visit_distribution_comparison$distribution_converged <- 
  ifelse(
    abs(round(resampled_HC_site_visit_distribution_comparison$control_site_pct, 2) - 
          round(resampled_HC_site_visit_distribution_comparison$comorbidity_site_pct, 2)) < 0.011,
    "Success", "Debug")

#1.3 Check whether the sampled control distributions at each site-event combination and whether that matches the distribution of the comorbidity subjects, and print a message detailing whether those were successful or not
print(
  if_else(
    length(unique(healthy_control_sample_resampled$src_subject_id)) == length(healthy_control_sample_resampled$src_subject_id) &&
      all(resampled_HC_site_visit_distribution_comparison$converged == "Success") ,
    "Control sample distribution is viable & contains only unique subjects",
    "Error: Controls not sampled correctly"))

#2. More scrupulously assess + cumulatively report that the expected number/proportion of subjects were successfully sampled within the HC group information for each site + timepoint
#2.11 Check whether the sampled control distributions at each site-event combination and whether that matches the distribution of the comorbidity subjects
resampled_HC_site_visit_distribution_comparison$final_cn_sample_unique <-
  if_else(
    length(unique(healthy_control_sample_resampled$src_subject_id)) == length(healthy_control_sample_resampled$src_subject_id) &&
      all(resampled_HC_site_visit_distribution_comparison$converged == "Success"), 
    TRUE, FALSE)

#2.12 Store a TRUE/FALSE value that signifies whether the sampled control distributions were successful or not
final_HC_sample_unique <-
  if_else(
    length(unique(healthy_control_sample_resampled$src_subject_id)) == length(healthy_control_sample_resampled$src_subject_id) &&
      all(resampled_HC_site_visit_distribution_comparison$converged == "Success") ,
    TRUE, FALSE)

#2.2 If the controls were sampled successfully, print a message detailing this; if they were not sampled correctly, print the subjects that were duplicated 
if_else(final_HC_sample_unique == TRUE,
        print("NA; Controls Sampled Successfully", print(healthy_control_sample_resampled[duplicated(healthy_control_sample_resampled$src_subject_id), "subjectkey"])))


## Output ## 

#1. Write the site + visit distribution comparison counts as a csv file
write.csv(resampled_HC_site_visit_distribution_comparison, "./data_processed/resampled_HC_site_visit_distribution_comparison.csv", row.names = FALSE)

#2. Write the final control and merged HC + comorbidity sample dataframes as csv files
#2.1 Write the final HC sample as a csv file
write.csv(healthy_control_sample_resampled, "./data_processed/healthy_control_group_resampled.csv", row.names = FALSE)

#2.2 Write the final merged HC + comorbidity sample as a csv file
write.csv(comorbidity_HC_sample, "./data_processed/comorbidity_HC_resampled_merged_groups.csv", row.names = FALSE)