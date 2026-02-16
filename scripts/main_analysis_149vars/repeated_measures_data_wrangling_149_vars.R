## Setup ##

library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tidyr)
options(digits = 6, scipen = 999)

## Paths ##

#1.1 Imaging data with all timepoints + all 149 vars
qcd_rsfmri_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_149vars.csv"

#1.2 Subject sample used for the project
sample_path <- "./data_processed/main_analysis/site_visit_analysis_data.csv"

#1.3 Significant DV list from the 149-var HC vs GAD step
sig_results_path <- "./results/results_149_vars/site_visit_significant_results_149vars.xlsx"

#1.4 Reproducible source of longitudinal trajectory group labels; this is the prior repeated measures dataset that includes: "Control", "GAD Converter", "GAD Remitter", and "Continuous GAD"
rm_group_source_path <- "./data_processed/main_analysis/repeated_measures_grouped_imaging_data.csv"

#1.5 Output
out_path <- "./data_processed/main_analysis/repeated_measures_grouped_imaging_data_149vars.csv"


## Data Wrangling ##

#2. Read in imaging data
qcd_rsfmri_data <- read.csv(qcd_rsfmri_path) %>%
  mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)),
    site_name = trimws(as.character(site_name)),
    sex = trimws(as.character(sex)))

#3. Read in the project sample, keep unique subject list + family/site/sex anchors
merged_sample <- read.csv(sample_path) %>%
  mutate(
    subjectkey = trimws(as.character(subjectkey)),
    site_name = trimws(as.character(site_name)),
    sex = trimws(as.character(sex)),
    rel_family_id = trimws(as.character(rel_family_id))) %>%
  distinct(subjectkey, .keep_all = TRUE) %>%
  rename(family_id = rel_family_id)

#3.1 Create an anchor table for optional validation checks (i.e., sex/site consistency)
merged_sample_anchors <- merged_sample %>%
  dplyr::select(
    subjectkey,
    family_id,
    sex_anchor = sex,
    site_anchor = site_name)

#4. Read in the significant results + DV list
sig_df <- readxl::read_xlsx(sig_results_path)

#4.1 Pull the column_name vector for the significant DVs from the corrected, expanded analysis
sig_dvs <- sig_df$column_name %>% as.character() %>% unique()

#4.2 Safety check for any missing DVs of interest
missing_dvs <- setdiff(sig_dvs, names(qcd_rsfmri_data))
if (length(missing_dvs) > 0) {
  stop(
    "These significant DV columns are missing from qcd_rsfmri_data: ",
    paste(missing_dvs, collapse = ", ")
  )
}

#5. Restrict imaging data to baseline + 2y follow-up only
qcd_rsfmri_data <- qcd_rsfmri_data %>%
  filter(eventname %in% c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"))

#6. Join sample subject list onto imaging data to recover BOTH timepoints per subject
#6.1 Retain only subjects in the matched sample, then attach family_id and anchors
repeated_measures_grouped_imaging_data <- qcd_rsfmri_data %>%
  semi_join(merged_sample_anchors, by = "subjectkey") %>%
  left_join(merged_sample_anchors, by = "subjectkey")

#6.2 Check whether imaging sex/site match sample anchors
#6.2.1 Create a df outlining sex mismatches w previous sample
sex_mismatch_n <- repeated_measures_grouped_imaging_data %>%
  filter(!is.na(sex), !is.na(sex_anchor), sex != sex_anchor) %>%
  nrow()

#6.2.2 Create a df outlining site mismatches w previous sample
site_mismatch_n <- repeated_measures_grouped_imaging_data %>%
  filter(!is.na(site_name), !is.na(site_anchor), site_name != site_anchor) %>%
  nrow()

#6.2.3 Print any sex mismatches
if (sex_mismatch_n > 0) {
  warning("Sex anchor mismatch detected in ", sex_mismatch_n, " rows after merge; investigate if unexpected.")
}

#6.2.4 Print any site mismatches
if (site_mismatch_n > 0) {
  warning("Site anchor mismatch detected in ", site_mismatch_n, " rows after merge; investigate if unexpected.")
}

#6.3 Drop anchor helper columns
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  dplyr::select(-sex_anchor, -site_anchor)

#6.4 Keep only subjects with both baseline and follow-up before joining trajectory labels
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(
    any(eventname == "baseline_year_1_arm_1") &
      any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()

#6.5 Apply rsfMRI QC inclusion if available (highly recommended)
if ("imgincl_rsfmri_include" %in% names(repeated_measures_grouped_imaging_data)) {
  repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
    mutate(imgincl_rsfmri_include = suppressWarnings(as.numeric(imgincl_rsfmri_include))) %>%
    filter(imgincl_rsfmri_include == 1)
}

#7. Create age in years
repeated_measures_grouped_imaging_data$age_in_years <- floor(
  as.numeric(as.character(repeated_measures_grouped_imaging_data$interview_age)) / 12)

#8. Attach longitudinal trajectory group labels from the reproducible repeated-measures dataset
#8.1 Read the prior repeated measures dataset and build a subject-level group lookup
rm_group_source <- read.csv(rm_group_source_path) %>%
  mutate(
    subjectkey = trimws(as.character(subjectkey)),
    group = trimws(as.character(group)))

#8.2 Create a subject-level mapping (subjectkey -> group) and enforce uniqueness
#8.2.1 Create the dx group lookup table
rm_group_lookup <- rm_group_source %>%
  dplyr::select(subjectkey, group) %>%
  distinct()

#8.2.2 Identify any duplicated subjects
dup_group_subjects <- rm_group_lookup %>%
  count(subjectkey) %>%
  filter(n > 1)

#8.2.3 Tag and message any duplicated subjects
if (nrow(dup_group_subjects) > 0) {
  stop(
    "Some subjectkeys map to multiple trajectory groups in rm_group_source_path; cannot safely join. ",
    "Example subjectkeys: ",
    paste(head(dup_group_subjects$subjectkey, 10), collapse = ", ")
  )
}

#8.3 Join the trajectory group labels onto the 149-var repeated measures dataset
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  inner_join(rm_group_lookup, by = "subjectkey")

#8.4 Ensure no missing group labels after join (should be zero after inner_join)
#8.4.1 Create a vector outlining the missingness of the group var in the merged data
missing_group_n <- repeated_measures_grouped_imaging_data %>%
  filter(is.na(group) | group == "") %>%
  nrow()

#8.4.2 If anyone is missing the group variable, list here
if (missing_group_n > 0) {
  stop("Unexpected: missing group labels remain after inner_join().")
}

#9. Keep only columns needed for repeated measures modeling
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  dplyr::select(
    subjectkey,
    eventname,
    group,
    interview_age,
    age_in_years,
    sex,
    site_name,
    family_id,
    rsfmri_c_ngd_meanmotion,
    all_of(sig_dvs)
  )

#10. Filter to subjects with both baseline and follow-up
repeated_measures_grouped_imaging_data <- repeated_measures_grouped_imaging_data %>%
  group_by(subjectkey) %>%
  filter(
    any(eventname == "baseline_year_1_arm_1") &
      any(eventname == "2_year_follow_up_y_arm_1")) %>%
  ungroup()

## Output ##

write.csv(repeated_measures_grouped_imaging_data, out_path, row.names = FALSE)
