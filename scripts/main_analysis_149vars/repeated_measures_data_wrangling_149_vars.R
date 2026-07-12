## Setup ##

# Load packages for reading, wrangling, and writing data
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tidyr)

# Set numeric display options
options(digits = 6, scipen = 999)


## Paths ##

#1.1 Corrected and QC-filtered imaging data with all 149 connectivity metrics
qcd_rsfmri_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_149vars.csv"

#1.2 Unfiltered raw rsfMRI data
raw_rsfmri_path <- "./data_raw/ABCD_rsfMRI_Data.csv"

#1.3 Official mapping of previously used to corrected rsfMRI column names
column_mapping_path <- "./data_raw/rsfmri_gpnet_aseg_correction.csv"

#1.4 Significant DV list from the corrected 149-variable HC versus GAD analysis
sig_results_path <- "./results/results_149_vars/site_visit_significant_results_149vars.xlsx"

#1.5 Published repeated measures dataset containing the longitudinal diagnostic trajectory groups
rm_group_source_path <- "./data_processed/main_analysis/repeated_measures_grouped_imaging_data.csv"

#1.6 Primary QC-complete repeated measures output
primary_out_path <- "./data_processed/main_analysis/repeated_measures_grouped_imaging_data_149vars.csv"

#1.7 published-cohort sensitivity output using corrected connectivity values
sensitivity_out_path <-
  paste0("./data_processed/main_analysis/",
    "repeated_measures_grouped_imaging_data_149vars_",
    "published_cohort_sensitivity.csv")

#1.8 Audit and sample-flow outputs
qc_row_audit_out_path <- "./results/results_149_vars/repeated_measures_qc_row_audit_149vars.csv"

qc_subject_audit_out_path <- "./results/results_149_vars/repeated_measures_qc_subject_audit_149vars.csv"

sample_flow_out_path <- "./results/results_149_vars/repeated_measures_sample_flow_149vars.csv"

dv_manifest_out_path <- "./results/results_149_vars/repeated_measures_dv_manifest_149vars.csv"


## Constants ##

#2.1 Define the assessment waves used in the repeated measures analyses
analysis_events <- c(
  "baseline_year_1_arm_1",
  "2_year_follow_up_y_arm_1")

#2.2 Define the expected primary sample counts
expected_primary_group_counts <- tibble::tibble(
  group = c(
    "Continuous GAD",
    "Control",
    "GAD Converter",
    "GAD Remitter"),
  expected_n = c(
    12L,
    1606L,
    46L,
    46L))

#2.3 Define the expected published-cohort sensitivity sample counts
expected_sensitivity_group_counts <- tibble::tibble(
  group = c(
    "Continuous GAD",
    "Control",
    "GAD Converter",
    "GAD Remitter"),
  expected_n = c(
    12L,
    1955L,
    55L,
    52L))


## Validation Function ##

#3. Validate the structure, completeness, and group counts of a repeated
# measures dataset
validate_repeated_measures_dataset <- function(
    data,
    dataset_name,
    expected_n_subjects,
    expected_n_rows,
    expected_group_counts,
    required_vars) {
  
  #3.1 Confirm that subject-event keys are unique
  duplicate_keys <- data %>%
    dplyr::count(subjectkey, eventname) %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_keys) > 0) {
    stop(
      dataset_name,
      " contains duplicated subjectkey-eventname rows.")
  }
  
  #3.2 Confirm that every participant has exactly the two expected waves
  subject_wave_check <- data %>%
    dplyr::group_by(subjectkey) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_events = dplyr::n_distinct(eventname),
      has_baseline =
        any(eventname == "baseline_year_1_arm_1"),
      has_followup =
        any(eventname == "2_year_follow_up_y_arm_1"),
      .groups = "drop")
  
  invalid_subject_waves <- subject_wave_check %>%
    dplyr::filter(
      n_rows != 2 |
        n_events != 2 |
        !has_baseline |
        !has_followup)
  
  if (nrow(invalid_subject_waves) > 0) {
    stop(
      dataset_name,
      " contains participants without exactly one baseline and one ",
      "two-year follow-up row.")
  }
  
  #3.3 Confirm the expected total numbers of participants and rows
  observed_n_subjects <- dplyr::n_distinct(data$subjectkey)
  observed_n_rows <- nrow(data)
  
  if (observed_n_subjects != expected_n_subjects) {
    stop(
      dataset_name,
      " contains ",
      observed_n_subjects,
      " participants; expected ",
      expected_n_subjects,
      ".")
  }
  
  if (observed_n_rows != expected_n_rows) {
    stop(
      dataset_name,
      " contains ",
      observed_n_rows,
      " rows; expected ",
      expected_n_rows,
      ".")
  }
  
  #3.4 Confirm the expected diagnostic trajectory group counts
  observed_group_counts <- data %>%
    dplyr::distinct(subjectkey, group) %>%
    dplyr::count(group, name = "observed_n")
  
  group_count_check <- expected_group_counts %>%
    dplyr::full_join(
      observed_group_counts,
      by = "group") %>%
    dplyr::mutate(
      expected_n = tidyr::replace_na(expected_n, 0L),
      observed_n = tidyr::replace_na(observed_n, 0L),
      matches_expected = expected_n == observed_n)
  
  if (any(!group_count_check$matches_expected)) {
    print(group_count_check)
    stop(
      dataset_name,
      " does not have the expected diagnostic trajectory group counts.")
  }
  
  #3.5 Confirm completeness across all required model variables
  missingness_check <- tibble::tibble(
    variable = required_vars,
    n_missing = vapply(
      required_vars,
      function(variable_name) {
        variable_values <- data[[variable_name]]
        sum(
          is.na(variable_values) |
            (is.character(variable_values) &
                trimws(variable_values) == ""))
      },
      integer(1)))
  
  if (any(missingness_check$n_missing > 0)) {
    print(
      missingness_check %>%
        dplyr::filter(n_missing > 0))
    stop(
      dataset_name,
      " contains missing values in one or more required model variables.")
  }
  
  message(
    dataset_name,
    " validated successfully: ",
    observed_n_subjects,
    " participants and ",
    observed_n_rows,
    " rows.")
  invisible(group_count_check)
}


## Significant Connectivity Outcomes ##

#4. Read the FDR-significant results from the corrected 149-variable analysis
sig_df <- readxl::read_xlsx(sig_results_path)

#4.1 Confirm that the expected column identifying connectivity metrics exists
if (!("column_name" %in% names(sig_df))) {
  stop(
    "The significant results file does not contain a column_name variable.")
}

#4.2 Create the ordered vector of significant connectivity outcomes
sig_dvs <- sig_df$column_name %>%
  as.character() %>%
  trimws() %>%
  unique()

sig_dvs <- sig_dvs[
  !is.na(sig_dvs) &
    sig_dvs != ""]

#4.3 Confirm that exactly 20 significant connectivity metrics were retained
if (length(sig_dvs) != 20) {
  stop(
    "Expected 20 significant connectivity metrics, but found ",
    length(sig_dvs), ".")
}

#4.4 Write the ordered connectivity outcome manifest
readr::write_csv(
  tibble::tibble(
    analysis_order = seq_along(sig_dvs),
    column_name = sig_dvs),
  dv_manifest_out_path)


## Published Repeated Measures Cohort Scaffold ##

#5. Read the previous repeated measures dataset
rm_group_source <- read.csv(
  rm_group_source_path,
  stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)),
    group = trimws(as.character(group)),
    sex = trimws(as.character(sex)),
    site_name = trimws(as.character(site_name)),
    family_id = trimws(as.character(family_id))) %>%
  dplyr::filter(
    eventname %in% analysis_events)

#5.1 Confirm that the required scaffold columns are available
required_scaffold_vars <- c(
  "subjectkey",
  "eventname",
  "group",
  "interview_age",
  "sex",
  "site_name",
  "family_id",
  "rsfmri_c_ngd_meanmotion")

missing_scaffold_vars <- setdiff(
  required_scaffold_vars,
  names(rm_group_source))

if (length(missing_scaffold_vars) > 0) {
  stop(
    "The published repeated measures source is missing: ",
    paste(missing_scaffold_vars, collapse = ", "))
}

#5.2 Create the subject-event scaffold from the published analysis
published_rm_scaffold <- rm_group_source %>%
  dplyr::select(
    dplyr::all_of(required_scaffold_vars)) %>%
  dplyr::distinct()

#5.3 Confirm that the published scaffold has unique subject-event rows
published_scaffold_duplicates <- published_rm_scaffold %>%
  dplyr::count(subjectkey, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(published_scaffold_duplicates) > 0) {
  stop(
    "The published repeated measures scaffold contains duplicated ",
    "subject-event rows.")
}

#5.4 Create a subject-level lookup for group and family ID
rm_subject_lookup <- published_rm_scaffold %>%
  dplyr::select(
    subjectkey,
    group,
    family_id) %>%
  dplyr::distinct()

#5.5 Confirm that each participant maps to one group and family ID
duplicated_subject_lookup <- rm_subject_lookup %>%
  dplyr::count(subjectkey) %>%
  dplyr::filter(n > 1)

if (nrow(duplicated_subject_lookup) > 0) {
  stop(
    "Some participants map to multiple diagnostic groups or family IDs ",
    "in the published repeated measures source.")
}


## Primary QC-Complete Repeated Measures Dataset ##

#6. Read the corrected and QC-filtered 149-variable imaging dataset
qcd_rsfmri_data <- read.csv(
  qcd_rsfmri_path,
  stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)),
    site_name = trimws(as.character(site_name)),
    sex = trimws(as.character(sex))) %>%
  dplyr::filter(
    eventname %in% analysis_events)

#6.1 Confirm that all required primary imaging variables exist
required_qcd_vars <- c(
  "subjectkey",
  "eventname",
  "interview_age",
  "sex",
  "site_name",
  "imgincl_rsfmri_include",
  "rsfmri_c_ngd_meanmotion",
  sig_dvs)

missing_qcd_vars <- setdiff(
  required_qcd_vars,
  names(qcd_rsfmri_data))

if (length(missing_qcd_vars) > 0) {
  stop(
    "The corrected QC-filtered imaging data are missing: ",
    paste(missing_qcd_vars, collapse = ", "))
}

#6.2 Confirm that the corrected QC-filtered data have unique subject-event rows
qcd_key_duplicates <- qcd_rsfmri_data %>%
  dplyr::count(subjectkey, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(qcd_key_duplicates) > 0) {
  stop(
    "The corrected QC-filtered imaging data contain duplicated ",
    "subject-event rows.")
}

#6.3 Retain the repeated measures cohort and attach trajectory group and family ID
primary_repeated_measures_data <- qcd_rsfmri_data %>%
  dplyr::select(
    subjectkey,
    eventname,
    interview_age,
    sex,
    site_name,
    imgincl_rsfmri_include,
    rsfmri_c_ngd_meanmotion,
    dplyr::all_of(sig_dvs)) %>%
  dplyr::semi_join(
    rm_subject_lookup,
    by = "subjectkey") %>%
  dplyr::left_join(
    rm_subject_lookup,
    by = "subjectkey")

#6.4 Standardize empty model variables to missing values
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::mutate(
    sex = dplyr::na_if(trimws(as.character(sex)), ""),
    site_name =
      dplyr::na_if(trimws(as.character(site_name)), ""),
    family_id =
      dplyr::na_if(trimws(as.character(family_id)), ""))

#6.5 Convert imaging inclusion, age, motion, and connectivity metrics to numeric
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::mutate(
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include)),
    interview_age =
      suppressWarnings(as.numeric(interview_age)),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    dplyr::across(
      dplyr::all_of(sig_dvs),
      ~ suppressWarnings(as.numeric(.x))))

#6.6 Require rsfMRI inclusion at each retained assessment wave
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::filter(
    imgincl_rsfmri_include == 1)

#6.7 Require simultaneous completeness across every model covariate and all 20 connectivity outcomes
primary_required_vars <- c(
  "interview_age",
  "sex",
  "site_name",
  "family_id",
  "rsfmri_c_ngd_meanmotion",
  sig_dvs)

primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(primary_required_vars),
      ~ !is.na(.x)))

#6.8 Require one baseline and one two-year follow-up row after applying QC and simultaneous complete-case restrictions
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::group_by(subjectkey) %>%
  dplyr::filter(
    dplyr::n() == 2,
    dplyr::n_distinct(eventname) == 2,
    all(analysis_events %in% eventname)) %>%
  dplyr::ungroup()

#6.9 Create age in years
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
  dplyr::mutate(
    age_in_years = floor(interview_age / 12))

#6.10 Retain and order the primary modeling variables
primary_repeated_measures_data <-
  primary_repeated_measures_data %>%
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
    dplyr::all_of(sig_dvs))


## Published-Cohort Corrected Sensitivity Dataset ##

#7. Read and clean the official rsfMRI column correction mapping
map_raw <- readr::read_csv(
  column_mapping_path,
  show_col_types = FALSE)

mapping <- map_raw %>%
  dplyr::select(
    `name_nda (previously used)`,
    `name_nda (correct)`) %>%
  dplyr::filter(
    !is.na(`name_nda (previously used)`),
    !is.na(`name_nda (correct)`)) %>%
  dplyr::distinct() %>%
  dplyr::rename(
    old = `name_nda (previously used)`,
    correct = `name_nda (correct)`) %>%
  dplyr::mutate(
    old = trimws(old),
    correct = trimws(correct))

#7.1 Ensure that each old column name maps to one corrected destination
duplicate_old_names <- mapping %>%
  dplyr::count(old) %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_old_names) > 0) {
  stop(
    "The official correction mapping contains duplicated old column names.")
}

old_to_correct <- stats::setNames(
  mapping$correct,
  mapping$old)

#7.2 Read the unfiltered raw rsfMRI data
raw_rsfmri_corrected <- readr::read_csv(
  raw_rsfmri_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname))) %>%
  dplyr::filter(eventname %in% analysis_events)

#7.3 Rename the raw headers using the official correction mapping
corrected_names <- old_to_correct[names(raw_rsfmri_corrected)]

names(raw_rsfmri_corrected)[!is.na(corrected_names)] <- corrected_names[!is.na(corrected_names)]

#7.4 Confirm that header correction did not create duplicated column names
duplicated_corrected_headers <- unique(
  names(raw_rsfmri_corrected)[
    duplicated(names(raw_rsfmri_corrected))])

if (length(duplicated_corrected_headers) > 0) {
  stop(
    "Correcting the raw rsfMRI headers produced duplicated names: ",
    paste(duplicated_corrected_headers, collapse = ", "))
}

#7.5 Confirm that the corrected raw data have unique subject-event rows
raw_key_duplicates <- raw_rsfmri_corrected %>%
  dplyr::count(subjectkey, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(raw_key_duplicates) > 0) {
  stop(
    "The corrected raw rsfMRI data contain duplicated subject-event rows.")
}

#7.6 Confirm that all 20 significant outcomes exist in the corrected raw data
missing_raw_dvs <- setdiff(
  sig_dvs,
  names(raw_rsfmri_corrected))

if (length(missing_raw_dvs) > 0) {
  stop(
    "The corrected raw rsfMRI data are missing: ",
    paste(missing_raw_dvs, collapse = ", "))
}

#7.7 Create the corrected raw connectivity table
corrected_raw_fc <- raw_rsfmri_corrected %>%
  dplyr::select(
    subjectkey,
    eventname,
    dplyr::all_of(sig_dvs)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sig_dvs),
      ~ suppressWarnings(as.numeric(.x))))

#7.8 Join corrected connectivity values onto the published scaffold
sensitivity_repeated_measures_data <-
  published_rm_scaffold %>%
  dplyr::left_join(
    corrected_raw_fc,
    by = c("subjectkey", "eventname"))

#7.9 Standardize and convert the published model covariates
sensitivity_repeated_measures_data <-
  sensitivity_repeated_measures_data %>%
  dplyr::mutate(
    sex = dplyr::na_if(trimws(as.character(sex)), ""),
    site_name =
      dplyr::na_if(trimws(as.character(site_name)), ""),
    family_id =
      dplyr::na_if(trimws(as.character(family_id)), ""),
    interview_age =
      suppressWarnings(as.numeric(interview_age)),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    age_in_years = floor(interview_age / 12))

#7.10 Require simultaneous completeness while preserving the exact
# published cohort
sensitivity_required_vars <- c(
  "interview_age",
  "sex",
  "site_name",
  "family_id",
  "rsfmri_c_ngd_meanmotion",
  sig_dvs)

incomplete_sensitivity_rows <-
  sensitivity_repeated_measures_data %>%
  dplyr::filter(
    !dplyr::if_all(
      dplyr::all_of(sensitivity_required_vars),
      ~ !is.na(.x)))

if (nrow(incomplete_sensitivity_rows) > 0) {
  stop(
    "The published-cohort sensitivity dataset contains ",
    nrow(incomplete_sensitivity_rows),
    " incomplete rows after attaching the corrected connectivity data.")
}

#7.11 Retain and order the sensitivity modeling variables
sensitivity_repeated_measures_data <-
  sensitivity_repeated_measures_data %>%
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
    dplyr::all_of(sig_dvs))


## Sample and QC Audits ##

#8. Create a row-level audit of availability in the corrected QC-filtered data
qcd_key_presence <- qcd_rsfmri_data %>%
  dplyr::select(
    subjectkey,
    eventname,
    imgincl_rsfmri_include) %>%
  dplyr::mutate(
    qcd_row_present = TRUE)

repeated_measures_qc_row_audit <-
  published_rm_scaffold %>%
  dplyr::select(
    subjectkey,
    eventname,
    group) %>%
  dplyr::left_join(
    qcd_key_presence,
    by = c("subjectkey", "eventname")) %>%
  dplyr::mutate(
    qcd_row_present =
      tidyr::replace_na(qcd_row_present, FALSE),
    qc_status = dplyr::case_when(
      qcd_row_present ~
        "Present in corrected QC-filtered imaging data",
      !qcd_row_present ~
        "Absent from corrected QC-filtered imaging data",
      TRUE ~ NA_character_))

#8.1 Create a subject-level audit of baseline and follow-up row availability
repeated_measures_qc_subject_audit <-
  repeated_measures_qc_row_audit %>%
  dplyr::group_by(subjectkey, group) %>%
  dplyr::summarise(
    baseline_present_in_qcd = any(
      eventname == "baseline_year_1_arm_1" &
        qcd_row_present),
    followup_present_in_qcd = any(
      eventname == "2_year_follow_up_y_arm_1" &
        qcd_row_present),
    retained_in_primary =
      subjectkey %in%
      primary_repeated_measures_data$subjectkey,
    .groups = "drop")

#8.2 Create a sample-flow table for both analysis datasets
repeated_measures_sample_flow <-
  dplyr::bind_rows(
    sensitivity_repeated_measures_data %>%
      dplyr::distinct(subjectkey, group) %>%
      dplyr::count(group, name = "n_subjects") %>%
      dplyr::mutate(
        analysis_dataset =
          "Published cohort corrected sensitivity"),
    primary_repeated_measures_data %>%
      dplyr::distinct(subjectkey, group) %>%
      dplyr::count(group, name = "n_subjects") %>%
      dplyr::mutate(
        analysis_dataset =
          "Primary QC-complete cohort")) %>%
  dplyr::mutate(
    n_rows = n_subjects * 2L) %>%
  dplyr::select(analysis_dataset, group, n_subjects, n_rows)


## Final Validation ##

#9. Validate the primary QC-complete dataset
validate_repeated_measures_dataset(
  data = primary_repeated_measures_data,
  dataset_name = "Primary QC-complete repeated measures dataset",
  expected_n_subjects = 1710L,
  expected_n_rows = 3420L,
  expected_group_counts = expected_primary_group_counts,
  required_vars = c(
    "subjectkey",
    "eventname",
    "group",
    primary_required_vars))

#9.1 Validate the published-cohort sensitivity dataset
validate_repeated_measures_dataset(
  data = sensitivity_repeated_measures_data,
  dataset_name =
    "Published-cohort corrected repeated measures sensitivity dataset",
  expected_n_subjects = 2074L,
  expected_n_rows = 4148L,
  expected_group_counts = expected_sensitivity_group_counts,
  required_vars = c(
    "subjectkey",
    "eventname",
    "group",
    sensitivity_required_vars))


## Output ##

#10. Write the primary QC-complete repeated measures dataset
readr::write_csv(
  primary_repeated_measures_data,
  primary_out_path)

#10.1 Write the  published-cohort corrected sensitivity dataset
readr::write_csv(
  sensitivity_repeated_measures_data,
  sensitivity_out_path)

#10.2 Write the row-level QC audit
readr::write_csv(
  repeated_measures_qc_row_audit,
  qc_row_audit_out_path)

#10.3 Write the subject-level QC audit
readr::write_csv(
  repeated_measures_qc_subject_audit,
  qc_subject_audit_out_path)

#10.4 Write the repeated measures sample-flow table
readr::write_csv(
  repeated_measures_sample_flow,
  sample_flow_out_path)

#10.5 Print final output summary
message(
  "Primary output written to: ",
  primary_out_path)

message(
  "Sensitivity output written to: ",
  sensitivity_out_path)

message(
  "Primary sample: ",
  dplyr::n_distinct(primary_repeated_measures_data$subjectkey),
  " participants and ",
  nrow(primary_repeated_measures_data),
  " rows.")

message(
  "Published-cohort sensitivity sample: ",
  dplyr::n_distinct(sensitivity_repeated_measures_data$subjectkey),
  " participants and ",
  nrow(sensitivity_repeated_measures_data),
  " rows.")
