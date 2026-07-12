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

#1.1 Corrected 149-variable imaging subset retained only for validation.
# It is NOT used to define repeated-measures eligibility because its row
# universe is inherited from a restricted upstream subset.
qcd_rsfmri_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_149vars.csv"

#1.2 Unfiltered raw rsfMRI data
raw_rsfmri_path <- "./data_raw/ABCD_rsfMRI_Data.csv"

#1.3 Official mapping of previously used to corrected rsfMRI column names
column_mapping_path <- "./data_raw/rsfmri_gpnet_aseg_correction.csv"

#1.3.1 Direct rsfMRI QC source of truth
qc_source_path <- "./data_raw/abcd_imgincl01.csv"

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
    1661L,
    46L,
    48L))

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



## Corrected Raw rsfMRI and Direct QC Sources ##

#6. Read and clean the official rsfMRI column correction mapping
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

#6.1 Ensure that each old column name maps to one corrected destination
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

#6.2 Read the unfiltered raw rsfMRI data
raw_rsfmri_corrected <- readr::read_csv(
  raw_rsfmri_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname))) %>%
  dplyr::filter(eventname %in% analysis_events)

#6.3 Rename the raw headers using the official correction mapping.
# This happens only in memory on a freshly read raw file, so rerunning the
# script cannot rename an already corrected output a second time.
corrected_names <- old_to_correct[names(raw_rsfmri_corrected)]
names(raw_rsfmri_corrected)[!is.na(corrected_names)] <-
  corrected_names[!is.na(corrected_names)]

#6.4 Confirm that header correction did not create duplicated column names
duplicated_corrected_headers <- unique(
  names(raw_rsfmri_corrected)[
    duplicated(names(raw_rsfmri_corrected))])

if (length(duplicated_corrected_headers) > 0) {
  stop(
    "Correcting the raw rsfMRI headers produced duplicated names: ",
    paste(duplicated_corrected_headers, collapse = ", "))
}

#6.5 Confirm that the corrected raw data have unique subject-event rows
raw_key_duplicates <- raw_rsfmri_corrected %>%
  dplyr::count(subjectkey, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(raw_key_duplicates) > 0) {
  stop(
    "The corrected raw rsfMRI data contain duplicated subject-event rows.")
}

#6.6 Confirm that the raw rsfMRI table contains the subject-event keys
# and all 20 corrected connectivity outcomes. Non-connectivity model
# covariates are retained from the exact published repeated-measures scaffold.
required_raw_vars <- c(
  "subjectkey",
  "eventname",
  sig_dvs)

missing_raw_vars <- setdiff(
  required_raw_vars,
  names(raw_rsfmri_corrected))

if (length(missing_raw_vars) > 0) {
  stop(
    "The corrected raw rsfMRI data are missing: ",
    paste(missing_raw_vars, collapse = ", "))
}

#6.7 Create the corrected raw connectivity source. The raw rsfMRI table is the source of truth for the corrected FC values. All non-connectivity model covariates remain tied to the published subject-event scaffold so that cohort inclusion is the only intended change
corrected_raw_fc <- raw_rsfmri_corrected %>%
  dplyr::select(
    subjectkey,
    eventname,
    dplyr::all_of(sig_dvs)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sig_dvs),
      ~ suppressWarnings(as.numeric(.x))))

#6.8 Read the direct ABCD rsfMRI inclusion variable
direct_qc <- readr::read_csv(
  qc_source_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  dplyr::transmute(
    subjectkey = trimws(as.character(src_subject_id)),
    eventname = trimws(as.character(eventname)),
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include))) %>%
  dplyr::filter(eventname %in% analysis_events)

direct_qc_duplicates <- direct_qc %>%
  dplyr::count(subjectkey, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(direct_qc_duplicates) > 0) {
  stop(
    "The direct rsfMRI QC source contains duplicated subject-event rows.")
}


## Primary QC-Complete Repeated Measures Dataset ##

#7. Begin with the exact published subject-event scaffold, then attach:
#   - row-specific corrected raw imaging/covariate data; and
#   - the direct rsfMRI QC flag.
# Presence in subset_qcd_imaging_data_149vars.csv is deliberately not an
# eligibility condition.
primary_row_universe <- published_rm_scaffold %>%
  dplyr::mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)),
    group = trimws(as.character(group)),
    interview_age =
      suppressWarnings(as.numeric(interview_age)),
    sex =
      dplyr::na_if(trimws(as.character(sex)), ""),
    site_name =
      dplyr::na_if(trimws(as.character(site_name)), ""),
    family_id =
      dplyr::na_if(trimws(as.character(family_id)), ""),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion))) %>%
  dplyr::left_join(
    corrected_raw_fc,
    by = c("subjectkey", "eventname")) %>%
  dplyr::left_join(
    direct_qc,
    by = c("subjectkey", "eventname"))

if (nrow(primary_row_universe) != nrow(published_rm_scaffold)) {
  stop(
    "Joining raw imaging and direct QC changed the number of published ",
    "subject-event scaffold rows.")
}

#7.1 Standardize the directly attached QC variable
primary_row_universe <- primary_row_universe %>%
  dplyr::mutate(
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include)))

#7.2 Define simultaneous model completeness
primary_required_vars <- c(
  "interview_age",
  "sex",
  "site_name",
  "family_id",
  "rsfmri_c_ngd_meanmotion",
  sig_dvs)

primary_row_universe <- primary_row_universe %>%
  dplyr::mutate(
    row_model_complete =
      dplyr::if_all(
        dplyr::all_of(primary_required_vars),
        ~ !is.na(.x)),
    row_qc_pass =
      !is.na(imgincl_rsfmri_include) &
      imgincl_rsfmri_include == 1)

#7.3 Apply the direct row-level QC and simultaneous complete-case rules
primary_repeated_measures_data <- primary_row_universe %>%
  dplyr::filter(
    row_qc_pass,
    row_model_complete) %>%
  dplyr::group_by(subjectkey) %>%
  dplyr::filter(
    dplyr::n() == 2,
    dplyr::n_distinct(eventname) == 2,
    all(analysis_events %in% eventname)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    age_in_years = floor(interview_age / 12)) %>%
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

#8. The corrected raw connectivity table was created in Section 6.7

#8.1 Join corrected connectivity values onto the exact published scaffold
sensitivity_repeated_measures_data <-
  published_rm_scaffold %>%
  dplyr::left_join(
    corrected_raw_fc,
    by = c("subjectkey", "eventname"))

#8.2 Standardize and convert the published model covariates
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

#8.3 Require simultaneous completeness while preserving the exact
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

#8.4 Retain and order the sensitivity modeling variables
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


## Validation Against Existing Pipeline Artifacts ##

#9. Confirm that every primary row belongs to the published sensitivity scaffold
unexpected_primary_keys <- primary_repeated_measures_data %>%
  dplyr::anti_join(
    sensitivity_repeated_measures_data,
    by = c("subjectkey", "eventname"))

if (nrow(unexpected_primary_keys) > 0) {
  stop(
    "The rebuilt primary cohort contains subject-event rows outside the ",
    "published repeated-measures scaffold.")
}

#9.1 Confirm directly that every retained primary row has QC == 1
primary_qc_check <- primary_repeated_measures_data %>%
  dplyr::select(subjectkey, eventname) %>%
  dplyr::left_join(
    direct_qc,
    by = c("subjectkey", "eventname"))

if (
  nrow(primary_qc_check) != nrow(primary_repeated_measures_data) ||
    any(is.na(primary_qc_check$imgincl_rsfmri_include)) ||
    any(primary_qc_check$imgincl_rsfmri_include != 1)) {
  stop(
    "At least one retained primary row does not have a direct rsfMRI QC flag ",
    "equal to 1.")
}

#9.2 Validate corrected raw connectivity against the existing corrected
# 149-variable subset wherever their subject-event rows overlap
if (file.exists(qcd_rsfmri_path)) {
  qcd_validation_source <- readr::read_csv(
    qcd_rsfmri_path,
    show_col_types = FALSE) %>%
    dplyr::mutate(
      subjectkey = trimws(as.character(subjectkey)),
      eventname = trimws(as.character(eventname)),
      interview_age =
        suppressWarnings(as.numeric(interview_age)),
      sex =
        dplyr::na_if(trimws(as.character(sex)), ""),
      site_name =
        dplyr::na_if(trimws(as.character(site_name)), ""),
      rsfmri_c_ngd_meanmotion =
        suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
      dplyr::across(
        dplyr::all_of(sig_dvs),
        ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::select(
      subjectkey,
      eventname,
      interview_age,
      sex,
      site_name,
      rsfmri_c_ngd_meanmotion,
      dplyr::all_of(sig_dvs))

  qcd_validation_duplicates <- qcd_validation_source %>%
    dplyr::count(subjectkey, eventname) %>%
    dplyr::filter(n > 1)

  if (nrow(qcd_validation_duplicates) > 0) {
    stop(
      "The existing corrected 149-variable subset contains duplicated ",
      "subject-event rows.")
  }

  overlap_validation <- corrected_raw_fc %>%
    dplyr::inner_join(
      qcd_validation_source,
      by = c("subjectkey", "eventname"),
      suffix = c("_raw", "_qcd149"))

  tolerance <- 1e-10

  overlap_mismatches <- vapply(
    sig_dvs,
    function(dv) {
      raw_values <- overlap_validation[[paste0(dv, "_raw")]]
      qcd_values <- overlap_validation[[paste0(dv, "_qcd149")]]

      sum(
        !(
          (is.na(raw_values) & is.na(qcd_values)) |
            (!is.na(raw_values) &
               !is.na(qcd_values) &
               abs(raw_values - qcd_values) <= tolerance)))
    },
    integer(1))

  if (any(overlap_mismatches > 0)) {
    stop(
      "Corrected raw connectivity values do not match the existing corrected ",
      "149-variable subset for: ",
      paste(names(overlap_mismatches)[overlap_mismatches > 0],
        collapse = ", "))
  }

  eligible_row_keys <- primary_row_universe %>%
    dplyr::filter(
      row_qc_pass,
      row_model_complete) %>%
    dplyr::select(subjectkey, eventname)

  eligible_rows_missing_from_qcd149 <- eligible_row_keys %>%
    dplyr::anti_join(
      qcd_validation_source,
      by = c("subjectkey", "eventname"))

  message(
    "Eligible scaffold rows absent from the corrected 149-variable subset: ",
    nrow(eligible_rows_missing_from_qcd149),
    ". These rows are recovered from corrected raw data for the primary cohort.")
}

#9.2.1 Validate that scaffold model covariates agree with the existing
# corrected imaging subset wherever both sources contain the row
scaffold_qcd_covariate_overlap <- published_rm_scaffold %>%
  dplyr::transmute(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)),
    interview_age_scaffold =
      suppressWarnings(as.numeric(interview_age)),
    sex_scaffold =
      dplyr::na_if(trimws(as.character(sex)), ""),
    site_name_scaffold =
      dplyr::na_if(trimws(as.character(site_name)), ""),
    meanmotion_scaffold =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion))) %>%
  dplyr::inner_join(
    qcd_validation_source %>%
      dplyr::transmute(
        subjectkey,
        eventname,
        interview_age_qcd = interview_age,
        sex_qcd = sex,
        site_name_qcd = site_name,
        meanmotion_qcd = rsfmri_c_ngd_meanmotion),
    by = c("subjectkey", "eventname"))

numeric_mismatch_count <- function(x, y, tolerance = 1e-10) {
  sum(
    !(
      (is.na(x) & is.na(y)) |
        (!is.na(x) &
           !is.na(y) &
           abs(x - y) <= tolerance)))
}

character_mismatch_count <- function(x, y) {
  x <- trimws(as.character(x))
  y <- trimws(as.character(y))
  
  sum(
    !(
      (is.na(x) & is.na(y)) |
        (!is.na(x) &
           !is.na(y) &
           x == y)))
}

scaffold_qcd_covariate_mismatches <- c(
  interview_age = numeric_mismatch_count(
    scaffold_qcd_covariate_overlap$interview_age_scaffold,
    scaffold_qcd_covariate_overlap$interview_age_qcd),
  sex = character_mismatch_count(
    scaffold_qcd_covariate_overlap$sex_scaffold,
    scaffold_qcd_covariate_overlap$sex_qcd),
  site_name = character_mismatch_count(
    scaffold_qcd_covariate_overlap$site_name_scaffold,
    scaffold_qcd_covariate_overlap$site_name_qcd),
  rsfmri_c_ngd_meanmotion = numeric_mismatch_count(
    scaffold_qcd_covariate_overlap$meanmotion_scaffold,
    scaffold_qcd_covariate_overlap$meanmotion_qcd))

if (any(scaffold_qcd_covariate_mismatches > 0)) {
  stop(
    "Published scaffold covariates do not match the corrected imaging ",
    "subset for: ",
    paste(
      names(scaffold_qcd_covariate_mismatches)[
        scaffold_qcd_covariate_mismatches > 0],
      collapse = ", "),
    ". No output files were overwritten.")
}

message(
  "Published scaffold covariates matched the corrected imaging subset ",
  "across ",
  nrow(scaffold_qcd_covariate_overlap),
  " overlapping subject-event rows.")

#9.3 If a prior primary output exists, ensure that rebuilding only adds
# eligible rows and does not alter values for rows already present. This check
# is idempotent: after the corrected output has been written, the existing and
# rebuilt datasets should match completely.
if (file.exists(primary_out_path)) {
  existing_primary <- readr::read_csv(
    primary_out_path,
    show_col_types = FALSE) %>%
    dplyr::mutate(
      subjectkey = trimws(as.character(subjectkey)),
      eventname = trimws(as.character(eventname)),
      group = trimws(as.character(group)),
      sex = trimws(as.character(sex)),
      site_name = trimws(as.character(site_name)),
      family_id = trimws(as.character(family_id)))

  existing_rows_removed <- existing_primary %>%
    dplyr::select(subjectkey, eventname) %>%
    dplyr::anti_join(
      primary_repeated_measures_data,
      by = c("subjectkey", "eventname"))

  if (nrow(existing_rows_removed) > 0) {
    stop(
      "The rebuilt primary cohort would remove ",
      nrow(existing_rows_removed),
      " rows from the existing primary output. No files were overwritten.")
  }

  comparison_columns <- c(
    "group",
    "interview_age",
    "age_in_years",
    "sex",
    "site_name",
    "family_id",
    "rsfmri_c_ngd_meanmotion",
    sig_dvs)

  missing_comparison_columns <- setdiff(
    comparison_columns,
    names(existing_primary))

  if (length(missing_comparison_columns) > 0) {
    stop(
      "The existing primary output is missing columns needed for the ",
      "retained-row parity check: ",
      paste(missing_comparison_columns, collapse = ", "))
  }

  shared_rows <- existing_primary %>%
    dplyr::select(
      subjectkey,
      eventname,
      dplyr::all_of(comparison_columns)) %>%
    dplyr::inner_join(
      primary_repeated_measures_data %>%
        dplyr::select(
          subjectkey,
          eventname,
          dplyr::all_of(comparison_columns)),
      by = c("subjectkey", "eventname"),
      suffix = c("_existing", "_rebuilt"))

  if (nrow(shared_rows) != nrow(existing_primary)) {
    stop(
      "The retained-row parity check did not recover every existing primary ",
      "row in the rebuilt cohort.")
  }

  numeric_comparison_columns <- c(
    "interview_age",
    "age_in_years",
    "rsfmri_c_ngd_meanmotion",
    sig_dvs)

  character_comparison_columns <- setdiff(
    comparison_columns,
    numeric_comparison_columns)

  retained_row_mismatch_counts <- c(
    stats::setNames(
      vapply(
        numeric_comparison_columns,
        function(variable_name) {
          existing_values <- suppressWarnings(as.numeric(
            shared_rows[[paste0(variable_name, "_existing")]]))
          rebuilt_values <- suppressWarnings(as.numeric(
            shared_rows[[paste0(variable_name, "_rebuilt")]]))

          sum(
            !(
              (is.na(existing_values) & is.na(rebuilt_values)) |
                (!is.na(existing_values) &
                   !is.na(rebuilt_values) &
                   abs(existing_values - rebuilt_values) <= 1e-10)))
        },
        integer(1)),
      numeric_comparison_columns),
    stats::setNames(
      vapply(
        character_comparison_columns,
        function(variable_name) {
          existing_values <- trimws(as.character(
            shared_rows[[paste0(variable_name, "_existing")]]))
          rebuilt_values <- trimws(as.character(
            shared_rows[[paste0(variable_name, "_rebuilt")]]))

          sum(
            !(
              (is.na(existing_values) & is.na(rebuilt_values)) |
                (!is.na(existing_values) &
                   !is.na(rebuilt_values) &
                   existing_values == rebuilt_values)))
        },
        integer(1)),
      character_comparison_columns))

  if (any(retained_row_mismatch_counts > 0)) {
    stop(
      "Rebuilding the primary cohort changed retained-row values for: ",
      paste(
        names(retained_row_mismatch_counts)[
          retained_row_mismatch_counts > 0],
        collapse = ", "),
      ". No files were overwritten.")
  }

  message(
    "Existing primary rows retained with identical model values: ",
    nrow(existing_primary),
    ".")
}


## Aggregate, ID-Free QC Audits ##

#10. Row-level aggregate QC audit. No subject IDs are written.
repeated_measures_qc_row_audit <- primary_row_universe %>%
  dplyr::mutate(
    qc_status = dplyr::case_when(
      is.na(imgincl_rsfmri_include) ~ "qc_flag_missing",
      imgincl_rsfmri_include == 1 ~ "qc_flag_equal_1",
      TRUE ~ "qc_flag_not_1"),
    model_data_status = dplyr::if_else(
      row_model_complete,
      "model_data_complete",
      "model_data_incomplete")) %>%
  dplyr::count(
    group,
    eventname,
    qc_status,
    model_data_status,
    name = "n_rows") %>%
  dplyr::arrange(
    group,
    eventname,
    qc_status,
    model_data_status)

#10.1 Subject-level aggregate QC audit. No subject IDs are written.
subject_qc_status <- primary_row_universe %>%
  dplyr::group_by(subjectkey, group) %>%
  dplyr::summarise(
    n_qc_passing_rows = sum(row_qc_pass),
    n_complete_rows = sum(row_model_complete),
    n_eligible_rows = sum(row_qc_pass & row_model_complete),
    retained_in_primary =
      subjectkey[[1]] %in%
      primary_repeated_measures_data$subjectkey,
    .groups = "drop") %>%
  dplyr::mutate(
    qc_pattern = dplyr::case_when(
      n_qc_passing_rows == 2 ~ "two_qc_passing_rows",
      n_qc_passing_rows == 1 ~ "one_qc_passing_row",
      n_qc_passing_rows == 0 ~ "no_qc_passing_rows",
      TRUE ~ "unexpected_qc_pattern"))

repeated_measures_qc_subject_audit <- subject_qc_status %>%
  dplyr::count(
    group,
    qc_pattern,
    n_complete_rows,
    n_eligible_rows,
    retained_in_primary,
    name = "n_subjects") %>%
  dplyr::arrange(
    group,
    qc_pattern,
    retained_in_primary)

#10.2 Create a sample-flow table for both analysis datasets
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
  dplyr::select(
    analysis_dataset,
    group,
    n_subjects,
    n_rows)


## Final Validation ##

#11. Validate the independently derived primary QC-complete dataset.
# The expected counts are regression checks applied only after eligibility
# has been derived from the scaffold, raw data, direct QC, and completeness.
validate_repeated_measures_dataset(
  data = primary_repeated_measures_data,
  dataset_name = "Primary QC-complete repeated measures dataset",
  expected_n_subjects = 1767L,
  expected_n_rows = 3534L,
  expected_group_counts = expected_primary_group_counts,
  required_vars = c(
    "subjectkey",
    "eventname",
    "group",
    primary_required_vars))

#11.1 Validate the published-cohort sensitivity dataset
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

#12. All validation occurs before any existing output is overwritten.
readr::write_csv(
  primary_repeated_measures_data,
  primary_out_path)

readr::write_csv(
  sensitivity_repeated_measures_data,
  sensitivity_out_path)

# These existing audit filenames now contain aggregate ID-free summaries.
readr::write_csv(
  repeated_measures_qc_row_audit,
  qc_row_audit_out_path)

readr::write_csv(
  repeated_measures_qc_subject_audit,
  qc_subject_audit_out_path)

readr::write_csv(
  repeated_measures_sample_flow,
  sample_flow_out_path)

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
