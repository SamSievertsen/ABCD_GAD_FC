## Setup ##

# Load packages for reading, wrangling, validation, and writing data
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

options(digits = 8, scipen = 999)


## Paths ##

# ABCD release 5.1 parent-report KSADS-COMP data
parent_ksads_path <- "./data_raw/mh_p_ksads_ss.csv"

# Parent-report KSADS variables used to define the HC pool
ksads_variable_manifest_path <- "./data_raw/supplementary_ksads_variable_desc.csv"

# Corrected, QC-filtered 149-variable imaging source
corrected_imaging_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_149vars.csv"

# Direct rs-fMRI QC source of truth
qc_source_path <- "./data_raw/abcd_imgincl01.csv"

# Ordered manifest of the 20 corrected outcomes carried forward
metric_manifest_path <-"./results/results_149_vars/repeated_measures_dv_manifest_149vars.csv"

# Output directory
output_dir <- "./data_processed/supplement_149vars"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

healthy_control_pool_out_path <- file.path(
  output_dir,
  "supplementary_healthy_control_pool.csv")

healthy_control_pool_summary_out_path <- file.path(
  output_dir,
  "supplementary_healthy_control_pool_summary.csv")


## Constants ##

analysis_events <- c(
  "baseline_year_1_arm_1",
  "2_year_follow_up_y_arm_1")


## Helper Functions ##

clean_character <- function(x) {
  dplyr::na_if(trimws(as.character(x)), "")
}

assert_columns <- function(data, required, dataset_name) {
  missing_columns <- setdiff(required, names(data))

  if (length(missing_columns) > 0) {
    stop(
      dataset_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      ".")
  }
}

assert_unique_keys <- function(data, dataset_name) {
  duplicate_keys <- data %>%
    dplyr::count(subjectkey, eventname, name = "n_rows") %>%
    dplyr::filter(n_rows > 1)

  if (nrow(duplicate_keys) > 0) {
    stop(
      dataset_name,
      " contains ",
      nrow(duplicate_keys),
      " duplicated subjectkey-eventname keys.")
  }
}

standardize_subject_id <- function(
    data,
    dataset_name,
    id_candidates = c("subjectkey", "src_subject_id")) {

  matching_ids <- intersect(id_candidates, names(data))

  if (length(matching_ids) == 0) {
    stop(
      dataset_name,
      " contains none of the supported participant-ID columns: ",
      paste(id_candidates, collapse = ", "),
      ".")
  }

  id_column <- matching_ids[[1]]

  data$subjectkey <- clean_character(data[[id_column]])

  if (any(is.na(data$subjectkey))) {
    stop(
      dataset_name,
      " contains missing participant identifiers after standardization.")
  }

  message(
    dataset_name,
    ": using ",
    id_column,
    " as the participant identifier.")

  data
}

clean_ksads_value <- function(x) {
  x <- clean_character(x)
  x[x %in% c("555", "999")] <- NA_character_
  x[x == "888"] <- "0"
  suppressWarnings(as.numeric(x))
}

write_csv_safely <- function(data, path) {
  output_directory <- dirname(path)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  temporary_path <- tempfile(
    pattern = paste0(basename(path), ".tmp_"),
    tmpdir = output_directory,
    fileext = ".csv")

  readr::write_csv(data, temporary_path)

  copied <- file.copy(
    from = temporary_path,
    to = path,
    overwrite = TRUE)

  unlink(temporary_path)

  if (!isTRUE(copied)) {
    stop("Failed to replace output file: ", path)
  }

  message("Output written: ", path)
}


## Connectivity Outcome Manifest ##

metric_manifest <- readr::read_csv(
  metric_manifest_path,
  show_col_types = FALSE) %>%
  dplyr::transmute(
    analysis_order = suppressWarnings(as.integer(analysis_order)),
    column_name = trimws(as.character(column_name))) %>%
  dplyr::arrange(analysis_order)

assert_columns(
  metric_manifest,
  c("analysis_order", "column_name"),
  "Connectivity outcome manifest")

if (
  nrow(metric_manifest) != 20L ||
    dplyr::n_distinct(metric_manifest$analysis_order) != 20L ||
    dplyr::n_distinct(metric_manifest$column_name) != 20L ||
    any(is.na(metric_manifest$analysis_order)) ||
    any(is.na(metric_manifest$column_name))) {
  stop("The connectivity outcome manifest must contain 20 unique outcomes.")
}

if (!identical(metric_manifest$analysis_order, seq_len(20L))) {
  stop("The connectivity manifest analysis_order must be exactly 1 through 20.")
}

connectivity_vars <- metric_manifest$column_name


## HC KSADS Variables ##

ksads_variable_manifest <- readr::read_csv(
  ksads_variable_manifest_path,
  show_col_types = FALSE)

assert_columns(
  ksads_variable_manifest,
  "variable_name",
  "HC KSADS variable manifest")

ksads_diagnosis_vars <- ksads_variable_manifest %>%
  dplyr::transmute(
    variable_name = trimws(as.character(variable_name))) %>%
  dplyr::filter(
    !is.na(variable_name),
    variable_name != "",
    stringr::str_starts(variable_name, "ksads_")) %>%
  dplyr::distinct() %>%
  dplyr::pull(variable_name)

if (length(ksads_diagnosis_vars) == 0L) {
  stop("No KSADS diagnosis variables were identified in the manifest.")
}


## Parent-Report KSADS Data ##

parent_ksads <- readr::read_csv(
  parent_ksads_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Parent-report KSADS source",
    id_candidates = c("src_subject_id", "subjectkey"))

assert_columns(
  parent_ksads,
  c("subjectkey", "eventname", ksads_diagnosis_vars),
  "Parent-report KSADS source")

parent_ksads <- parent_ksads %>%
  dplyr::mutate(
    eventname = clean_character(eventname)) %>%
  dplyr::filter(eventname %in% analysis_events) %>%
  dplyr::select(
    subjectkey,
    eventname,
    dplyr::all_of(ksads_diagnosis_vars)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(ksads_diagnosis_vars),
      clean_ksads_value))

assert_unique_keys(
  parent_ksads,
  "Filtered parent-report KSADS source")

unexpected_ksads_values <- parent_ksads %>%
  dplyr::select(dplyr::all_of(ksads_diagnosis_vars)) %>%
  unlist(use.names = FALSE) %>%
  unique() %>%
  stats::na.omit() %>%
  setdiff(c(0, 1))

if (length(unexpected_ksads_values) > 0) {
  stop(
    "Unexpected cleaned KSADS diagnosis values were found: ",
    paste(unexpected_ksads_values, collapse = ", "),
    ".")
}

diagnosis_matrix <- as.matrix(
  parent_ksads %>%
    dplyr::select(dplyr::all_of(ksads_diagnosis_vars)))

parent_ksads <- parent_ksads %>%
  dplyr::mutate(
    row_has_any_observed_diagnosis_item =
      rowSums(!is.na(diagnosis_matrix)) > 0,
    row_has_any_positive_diagnosis =
      rowSums(diagnosis_matrix == 1, na.rm = TRUE) > 0)

healthy_control_subjects <- parent_ksads %>%
  dplyr::group_by(subjectkey) %>%
  dplyr::summarise(
    has_any_observed_diagnosis_item =
      any(row_has_any_observed_diagnosis_item),
    has_any_positive_diagnosis =
      any(row_has_any_positive_diagnosis),
    .groups = "drop") %>%
  dplyr::filter(
    has_any_observed_diagnosis_item,
    !has_any_positive_diagnosis) %>%
  dplyr::select(subjectkey)

if (nrow(healthy_control_subjects) == 0L) {
  stop("No HC-eligible participants were identified.")
}


## Corrected Imaging Eligibility ##

corrected_imaging <- readr::read_csv(
  corrected_imaging_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Corrected 149-variable imaging source")

assert_columns(
  corrected_imaging,
  c(
    "subjectkey",
    "eventname",
    "rsfmri_c_ngd_meanmotion",
    connectivity_vars),
  "Corrected 149-variable imaging source")

corrected_imaging <- corrected_imaging %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    dplyr::across(
      dplyr::all_of(connectivity_vars),
      ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::filter(eventname %in% analysis_events)

assert_unique_keys(
  corrected_imaging,
  "Corrected 149-variable imaging source")

complete_imaging_keys <- corrected_imaging %>%
  dplyr::filter(
    !is.na(rsfmri_c_ngd_meanmotion),
    dplyr::if_all(
      dplyr::all_of(connectivity_vars),
      ~ !is.na(.x))) %>%
  dplyr::select(
    subjectkey,
    eventname)


## Direct rs-fMRI QC Validation ##

direct_qc <- readr::read_csv(
  qc_source_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  standardize_subject_id(
    dataset_name = "Direct rs-fMRI QC source",
    id_candidates = c("src_subject_id")) %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include))) %>%
  dplyr::filter(eventname %in% analysis_events) %>%
  dplyr::select(
    subjectkey,
    eventname,
    imgincl_rsfmri_include)

assert_columns(
  direct_qc,
  c(
    "subjectkey",
    "eventname",
    "imgincl_rsfmri_include"),
  "Direct rs-fMRI QC source")

assert_unique_keys(
  direct_qc,
  "Direct rs-fMRI QC source")

eligible_imaging_keys <- complete_imaging_keys %>%
  dplyr::left_join(
    direct_qc,
    by = c("subjectkey", "eventname"))

if (any(is.na(eligible_imaging_keys$imgincl_rsfmri_include))) {
  stop(
    "At least one corrected imaging row did not match the direct QC source.")
}

if (any(eligible_imaging_keys$imgincl_rsfmri_include != 1)) {
  stop(
    "The corrected QC-filtered imaging source contains one or more rows with ",
    "direct imgincl_rsfmri_include != 1.")
}

eligible_imaging_keys <- eligible_imaging_keys %>%
  dplyr::select(
    subjectkey,
    eventname)


## Final HC Pool ##

healthy_control_pool <- eligible_imaging_keys %>%
  dplyr::semi_join(
    healthy_control_subjects,
    by = "subjectkey") %>%
  dplyr::mutate(
    group = "HC",
    comorbidity_group = "HC",
    broad_clinical_group = "HC") %>%
  dplyr::arrange(
    subjectkey,
    factor(eventname, levels = analysis_events))

assert_unique_keys(
  healthy_control_pool,
  "Final HC pool")

if (nrow(healthy_control_pool) == 0L) {
  stop("The final HC imaging pool is empty.")
}

healthy_control_pool_summary <- dplyr::bind_rows(
  healthy_control_pool %>%
    dplyr::summarise(
      summary_level = "overall",
      eventname = "all",
      n_subject_event_rows = dplyr::n(),
      n_unique_subjects = dplyr::n_distinct(subjectkey)),
  healthy_control_pool %>%
    dplyr::group_by(eventname) %>%
    dplyr::summarise(
      summary_level = "event",
      n_subject_event_rows = dplyr::n(),
      n_unique_subjects = dplyr::n_distinct(subjectkey),
      .groups = "drop"))

message(
  "HC pool validated: ",
  nrow(healthy_control_pool),
  " subject-event rows from ",
  dplyr::n_distinct(healthy_control_pool$subjectkey),
  " unique participants.")

print(healthy_control_pool_summary)


## Output ##

write_csv_safely(
  healthy_control_pool,
  healthy_control_pool_out_path)

write_csv_safely(
  healthy_control_pool_summary,
  healthy_control_pool_summary_out_path)
