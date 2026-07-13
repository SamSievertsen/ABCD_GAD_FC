## Setup ##

# Load packages for reading, wrangling, validation, and writing data
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

options(digits = 8, scipen = 999)


## Paths ##

# Current cross-sectional sample used for the corrected 149-metric
site_visit_scaffold_path <- "./data_processed/main_analysis/site_visit_analysis_data.csv"

# Raw rsfMRI source and official DAIRC column-name correction mapping
raw_rsfmri_path <- "./data_raw/ABCD_rsfMRI_Data.csv"
column_mapping_path <- "./data_raw/rsfmri_gpnet_aseg_correction.csv"

# Direct rsfMRI QC source of truth
qc_source_path <- "./data_raw/abcd_imgincl01.csv"

# Corrected 149-variable source used only to validate corrected FC values for
# the selected scaffold rows
corrected_149_validation_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_149vars.csv"

# Ordered manifest of the 20 corrected cross-sectional outcomes carried forward
metric_manifest_path <- "./results/results_149_vars/repeated_measures_dv_manifest_149vars.csv"

# Symptom sources
cbcl_path <- "./data_raw/abcd_cbcls01.csv"
bpm_path <- "./data_raw/abcd_yssbpm01.csv"

# Existing dimensional output paths to replace only after all checks pass
cbcl_out_path <- "./data_processed/main_analysis/dimensional_analysis_imaging_cbcl_data.csv"
dimensional_out_path <- "./data_processed/main_analysis/dimensional_analysis_data.csv"


## Constants ##

analysis_events <- c(
  "baseline_year_1_arm_1",
  "2_year_follow_up_y_arm_1")

analysis_groups <- c("control", "GAD")

cbcl_vars <- c(
  "cbcl_scr_syn_anxdep_t",
  "cbcl_scr_syn_internal_t",
  "cbcl_scr_syn_external_t",
  "cbcl_scr_dsm5_anxdisord_t")

bpm_vars <- c(
  "bpm_y_scr_internal_t",
  "bpm_y_scr_external_t")

bpm_missing_item_vars <- c(
  "bpm_y_scr_internal_nm",
  "bpm_y_scr_external_nm")

# Current selected cross-sectional scaffold
expected_scaffold_n <- 3322L
expected_scaffold_group_counts <- tibble::tribble(
  ~analysis_group, ~expected_n,
  "control", 3158L,
  "GAD",      164L)

# Intended/published complete-case CBCL dimensional sample
expected_cbcl_n <- 3146L
expected_cbcl_group_event_counts <- tibble::tribble(
  ~analysis_group, ~eventname, ~expected_n,
  "control", "baseline_year_1_arm_1", 2003L,
  "control", "2_year_follow_up_y_arm_1", 986L,
  "GAD", "baseline_year_1_arm_1", 98L,
  "GAD", "2_year_follow_up_y_arm_1", 59L)


## Helper Functions ##

clean_character <- function(x) {
  dplyr::na_if(trimws(as.character(x)), "")
}

numeric_equal <- function(x, y, tolerance = 1e-10) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))

  (is.na(x) & is.na(y)) |
    (!is.na(x) & !is.na(y) & abs(x - y) <= tolerance)
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

assert_complete <- function(data, variables, dataset_name) {
  missing_counts <- vapply(
    variables,
    function(variable) {
      values <- data[[variable]]
      sum(
        is.na(values) |
          (is.character(values) & trimws(values) == ""))
    },
    integer(1))

  if (any(missing_counts > 0)) {
    print(missing_counts[missing_counts > 0])
    stop(
      dataset_name,
      " contains missing values in required variables.")
  }
}

assert_exact_group_counts <- function(
    data,
    expected_counts,
    grouping_variables,
    dataset_name) {

  observed_counts <- data %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(
      observed_n = dplyr::n(),
      .groups = "drop")

  count_check <- expected_counts %>%
    dplyr::full_join(
      observed_counts,
      by = grouping_variables) %>%
    dplyr::mutate(
      expected_n = tidyr::replace_na(expected_n, 0L),
      observed_n = tidyr::replace_na(observed_n, 0L),
      matches = expected_n == observed_n)

  if (any(!count_check$matches)) {
    print(count_check)
    stop(dataset_name, " does not have the expected counts.")
  }

  invisible(count_check)
}

# Use one explicitly selected participant identifier for each NDA table.
standardize_subject_id <- function(data, dataset_name, id_column) {
  if (!(id_column %in% names(data))) {
    stop(
      dataset_name,
      " is missing the required participant-ID column: ",
      id_column,
      ".")
  }

  selected_id <- clean_character(data[[id_column]])

  if (any(is.na(selected_id))) {
    stop(
      dataset_name,
      " contains ",
      sum(is.na(selected_id)),
      " rows with a missing ",
      id_column,
      " after removal of the description row.")
  }

  data$subjectkey <- selected_id

  message(
    dataset_name,
    ": using ",
    id_column,
    " as the participant identifier.")

  data
}

# Write to a temporary file in the destination directory, then replace the
# target only after the write completes successfully.
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
  show_col_types = FALSE)

assert_columns(
  metric_manifest,
  c("analysis_order", "column_name"),
  "Connectivity outcome manifest")

metric_manifest <- metric_manifest %>%
  dplyr::transmute(
    analysis_order = suppressWarnings(as.integer(analysis_order)),
    column_name = clean_character(column_name)) %>%
  dplyr::arrange(analysis_order)

if (
  nrow(metric_manifest) != 20L ||
    dplyr::n_distinct(metric_manifest$analysis_order) != 20L ||
    dplyr::n_distinct(metric_manifest$column_name) != 20L ||
    any(is.na(metric_manifest$analysis_order)) ||
    any(is.na(metric_manifest$column_name))) {
  stop(
    "The connectivity outcome manifest must contain exactly 20 unique, ",
    "nonmissing analysis orders and column names.")
}

if (!identical(metric_manifest$analysis_order, seq_len(20L))) {
  stop("The connectivity manifest analysis_order must be exactly 1 through 20.")
}

connectivity_vars <- metric_manifest$column_name


## Current Cross-Sectional Scaffold ##

site_visit_scaffold <- readr::read_csv(
  site_visit_scaffold_path,
  show_col_types = FALSE)

# Normalize the family-ID field before selecting the scaffold variables.
if ("family_id" %in% names(site_visit_scaffold)) {
  site_visit_scaffold <- site_visit_scaffold %>%
    dplyr::mutate(family_id = clean_character(family_id))
} else if ("rel_family_id" %in% names(site_visit_scaffold)) {
  site_visit_scaffold <- site_visit_scaffold %>%
    dplyr::rename(family_id = rel_family_id) %>%
    dplyr::mutate(family_id = clean_character(family_id))
} else {
  stop("site_visit_analysis_data.csv has neither family_id nor rel_family_id.")
}

required_scaffold_vars <- c(
  "subjectkey",
  "family_id",
  "eventname",
  "analysis_group",
  "group",
  "interview_age",
  "age_in_years",
  "sex",
  "site_name",
  "rsfmri_c_ngd_meanmotion")

assert_columns(
  site_visit_scaffold,
  required_scaffold_vars,
  "Current site-visit scaffold")

site_visit_scaffold <- site_visit_scaffold %>%
  dplyr::transmute(
    subjectkey = clean_character(subjectkey),
    family_id = clean_character(family_id),
    eventname = clean_character(eventname),
    analysis_group = clean_character(analysis_group),
    group = clean_character(group),
    interview_age = suppressWarnings(as.numeric(interview_age)),
    age_in_years = suppressWarnings(as.numeric(age_in_years)),
    sex = clean_character(sex),
    site_name = clean_character(site_name),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion))) %>%
  dplyr::filter(
    eventname %in% analysis_events,
    analysis_group %in% analysis_groups)

assert_unique_keys(site_visit_scaffold, "Current site-visit scaffold")
assert_complete(
  site_visit_scaffold,
  required_scaffold_vars,
  "Current site-visit scaffold")

if (nrow(site_visit_scaffold) != expected_scaffold_n) {
  stop(
    "The current site-visit scaffold contains ",
    nrow(site_visit_scaffold),
    " rows; expected ",
    expected_scaffold_n,
    ".")
}

assert_exact_group_counts(
  data = site_visit_scaffold,
  expected_counts = expected_scaffold_group_counts,
  grouping_variables = "analysis_group",
  dataset_name = "Current site-visit scaffold")

if (dplyr::n_distinct(site_visit_scaffold$subjectkey) != nrow(site_visit_scaffold)) {
  stop(
    "The current site-visit scaffold must contain exactly one selected ",
    "assessment per participant.")
}

message(
  "Current site-visit scaffold validated: ",
  expected_scaffold_n,
  " unique participants/subject-event rows.")


## Correct Raw rsfMRI Headers and Extract the 20 Outcomes ##

map_raw <- readr::read_csv(
  column_mapping_path,
  show_col_types = FALSE)

assert_columns(
  map_raw,
  c("name_nda (previously used)", "name_nda (correct)"),
  "Official DAIRC mapping")

mapping <- map_raw %>%
  dplyr::transmute(
    old = clean_character(`name_nda (previously used)`),
    correct = clean_character(`name_nda (correct)`)) %>%
  dplyr::filter(!is.na(old), !is.na(correct)) %>%
  dplyr::distinct()

mapping_conflicts <- mapping %>%
  dplyr::count(old, name = "n_destinations") %>%
  dplyr::filter(n_destinations > 1)

if (nrow(mapping_conflicts) > 0) {
  stop("The official correction mapping contains conflicting old column names.")
}

old_to_correct <- stats::setNames(mapping$correct, mapping$old)

raw_rsfmri_corrected <- readr::read_csv(
  raw_rsfmri_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate(
    subjectkey = clean_character(subjectkey),
    eventname = clean_character(eventname)) %>%
  dplyr::filter(eventname %in% analysis_events)

corrected_names <- old_to_correct[names(raw_rsfmri_corrected)]
name_positions_to_replace <- which(!is.na(corrected_names))
names(raw_rsfmri_corrected)[name_positions_to_replace] <-
  unname(corrected_names[name_positions_to_replace])

duplicated_corrected_headers <- unique(
  names(raw_rsfmri_corrected)[duplicated(names(raw_rsfmri_corrected))])

if (length(duplicated_corrected_headers) > 0) {
  stop(
    "Correcting raw rsfMRI headers produced duplicated names: ",
    paste(duplicated_corrected_headers, collapse = ", "),
    ".")
}

assert_columns(
  raw_rsfmri_corrected,
  c("subjectkey", "eventname", connectivity_vars),
  "Corrected raw rsfMRI data")
assert_unique_keys(raw_rsfmri_corrected, "Corrected raw rsfMRI data")

corrected_raw_fc <- raw_rsfmri_corrected %>%
  dplyr::select(
    subjectkey,
    eventname,
    dplyr::all_of(connectivity_vars)) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(connectivity_vars),
      ~ suppressWarnings(as.numeric(.x))))


## Direct rsfMRI QC ##

direct_qc <- readr::read_csv(
  qc_source_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  standardize_subject_id(
    dataset_name = "Direct rsfMRI QC source",
    id_column = "src_subject_id")

assert_columns(
  direct_qc,
  c("subjectkey", "eventname", "imgincl_rsfmri_include"),
  "Direct rsfMRI QC source")

direct_qc <- direct_qc %>%
  dplyr::transmute(
    subjectkey = clean_character(subjectkey),
    eventname = clean_character(eventname),
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include))) %>%
  dplyr::filter(eventname %in% analysis_events)

assert_unique_keys(direct_qc, "Direct rsfMRI QC source")


## Build and Validate the Corrected Imaging Scaffold ##

corrected_imaging_scaffold <- site_visit_scaffold %>%
  dplyr::left_join(
    corrected_raw_fc,
    by = c("subjectkey", "eventname")) %>%
  dplyr::left_join(
    direct_qc,
    by = c("subjectkey", "eventname"))

if (nrow(corrected_imaging_scaffold) != expected_scaffold_n) {
  stop("The imaging/QC joins unexpectedly changed the scaffold row count.")
}

assert_unique_keys(
  corrected_imaging_scaffold,
  "Corrected imaging scaffold")

missing_fc_counts <- vapply(
  connectivity_vars,
  function(variable) sum(is.na(corrected_imaging_scaffold[[variable]])),
  integer(1))

if (any(missing_fc_counts > 0)) {
  print(missing_fc_counts[missing_fc_counts > 0])
  stop(
    "One or more current site-visit observations lack corrected values for ",
    "at least one of the 20 connectivity outcomes.")
}

invalid_qc_rows <- corrected_imaging_scaffold %>%
  dplyr::filter(
    is.na(imgincl_rsfmri_include) |
      imgincl_rsfmri_include != 1)

if (nrow(invalid_qc_rows) > 0) {
  print(
    invalid_qc_rows %>%
      dplyr::count(
        analysis_group,
        eventname,
        imgincl_rsfmri_include,
        name = "n_rows"))
  stop(
    "At least one current site-visit observation does not have direct ",
    "imgincl_rsfmri_include == 1.")
}

message(
  "All ",
  expected_scaffold_n,
  " selected observations have complete corrected FC values and direct QC == 1.")


## Validate Selected Corrected Values Against the Existing 149-Variable Source ##

if (file.exists(corrected_149_validation_path)) {
  corrected_149_validation <- readr::read_csv(
    corrected_149_validation_path,
    show_col_types = FALSE)

  assert_columns(
    corrected_149_validation,
    c("subjectkey", "eventname", connectivity_vars),
    "Corrected 149-variable validation source")

  corrected_149_validation <- corrected_149_validation %>%
    dplyr::transmute(
      subjectkey = clean_character(subjectkey),
      eventname = clean_character(eventname),
      dplyr::across(
        dplyr::all_of(connectivity_vars),
        ~ suppressWarnings(as.numeric(.x))))

  assert_unique_keys(
    corrected_149_validation,
    "Corrected 149-variable validation source")

  selected_validation_overlap <- corrected_imaging_scaffold %>%
    dplyr::select(
      subjectkey,
      eventname,
      dplyr::all_of(connectivity_vars)) %>%
    dplyr::inner_join(
      corrected_149_validation,
      by = c("subjectkey", "eventname"),
      suffix = c("_raw", "_qcd149"))

  if (nrow(selected_validation_overlap) != expected_scaffold_n) {
    stop(
      "Only ",
      nrow(selected_validation_overlap),
      " of ",
      expected_scaffold_n,
      " selected rows were available in the corrected 149-variable validation source.")
  }

  validation_mismatch_counts <- vapply(
    connectivity_vars,
    function(metric) {
      sum(
        !numeric_equal(
          selected_validation_overlap[[paste0(metric, "_raw")]],
          selected_validation_overlap[[paste0(metric, "_qcd149")]]))
    },
    integer(1))

  if (any(validation_mismatch_counts > 0)) {
    print(validation_mismatch_counts[validation_mismatch_counts > 0])
    stop(
      "Corrected raw FC values disagree with the existing 149-variable ",
      "source for one or more selected outcomes.")
  }

  message(
    "Selected corrected raw FC values matched the 149-variable source across ",
    expected_scaffold_n,
    " rows and all 20 outcomes.")
}


## Clean and Join CBCL Data ##

# Remove the single NDA description row exactly once.
cbcl_raw <- readr::read_csv(
  cbcl_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  standardize_subject_id(
    dataset_name = "CBCL source",
    id_column = "subjectkey")

assert_columns(
  cbcl_raw,
  c("subjectkey", "eventname", cbcl_vars),
  "CBCL source")

cbcl_raw <- cbcl_raw %>%
  dplyr::transmute(
    subjectkey = clean_character(subjectkey),
    eventname = clean_character(eventname),
    dplyr::across(
      dplyr::all_of(cbcl_vars),
      ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::filter(eventname %in% analysis_events)

assert_unique_keys(cbcl_raw, "CBCL source")

cbcl_data <- corrected_imaging_scaffold %>%
  dplyr::left_join(
    cbcl_raw,
    by = c("subjectkey", "eventname")) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(cbcl_vars),
      ~ !is.na(.x))) %>%
  dplyr::select(
    subjectkey,
    family_id,
    eventname,
    analysis_group,
    group,
    interview_age,
    age_in_years,
    sex,
    site_name,
    rsfmri_c_ngd_meanmotion,
    dplyr::all_of(connectivity_vars),
    dplyr::all_of(cbcl_vars))

assert_unique_keys(cbcl_data, "Primary CBCL dimensional dataset")

if (nrow(cbcl_data) != expected_cbcl_n) {
  stop(
    "The primary CBCL complete-case sample contains ",
    nrow(cbcl_data),
    " rows; expected ",
    expected_cbcl_n,
    ".")
}

if (dplyr::n_distinct(cbcl_data$subjectkey) != expected_cbcl_n) {
  stop("The primary CBCL sample must contain one observation per participant.")
}

assert_exact_group_counts(
  data = cbcl_data,
  expected_counts = expected_cbcl_group_event_counts,
  grouping_variables = c("analysis_group", "eventname"),
  dataset_name = "Primary CBCL complete-case sample")

assert_complete(
  cbcl_data,
  c(
    required_scaffold_vars,
    connectivity_vars,
    cbcl_vars),
  "Primary CBCL dimensional dataset")

message(
  "Primary CBCL dimensional sample validated: ",
  expected_cbcl_n,
  " participants (2,989 controls; 157 GAD).")


## Clean BPM Data and Preserve the Published Nested Sample Definition ##

bpm_raw <- readr::read_csv(
  bpm_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  standardize_subject_id(
    dataset_name = "BPM source",
    id_column = "subjectkey")

assert_columns(
  bpm_raw,
  c("subjectkey", "eventname", bpm_vars, bpm_missing_item_vars),
  "BPM source")

bpm_complete_source <- bpm_raw %>%
  dplyr::transmute(
    subjectkey = clean_character(subjectkey),
    eventname = clean_character(eventname),
    dplyr::across(
      dplyr::all_of(c(bpm_vars, bpm_missing_item_vars)),
      ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::filter(
    eventname %in% analysis_events,
    !is.na(bpm_y_scr_internal_t),
    !is.na(bpm_y_scr_external_t),
    bpm_y_scr_internal_nm == 0,
    bpm_y_scr_external_nm == 0) %>%
  dplyr::select(
    subjectkey,
    eventname,
    dplyr::all_of(bpm_vars))

assert_unique_keys(bpm_complete_source, "Complete BPM source")

# Preserve the published approach: BPM analyses are nested within the complete
# four-score CBCL sample. Participants without complete BPM data remain in the
# combined file with missing BPM values and are excluded only from BPM models
dimensional_analysis_data <- cbcl_data %>%
  dplyr::left_join(
    bpm_complete_source,
    by = c("subjectkey", "eventname"))

assert_unique_keys(
  dimensional_analysis_data,
  "Primary combined dimensional dataset")

if (nrow(dimensional_analysis_data) != expected_cbcl_n) {
  stop("Joining BPM unexpectedly changed the primary dimensional row count.")
}

bpm_complete_data <- dimensional_analysis_data %>%
  dplyr::filter(
    !is.na(bpm_y_scr_internal_t),
    !is.na(bpm_y_scr_external_t))

if (nrow(bpm_complete_data) == 0L) {
  stop("No participants met the published BPM completeness definition.")
}

if (
  dplyr::n_distinct(bpm_complete_data$analysis_group) != 2L ||
    dplyr::n_distinct(bpm_complete_data$eventname) < 1L) {
  stop("The BPM-complete sample does not contain both analysis groups.")
}

message(
  "Primary BPM-complete sample nested within CBCL: ",
  nrow(bpm_complete_data),
  " participants.")
print(
  bpm_complete_data %>%
    dplyr::count(analysis_group, eventname, name = "n"))

if (nrow(bpm_complete_data) != 974L) {
  message(
    "Informational: the current primary BPM-complete sample contains ",
    nrow(bpm_complete_data),
    " participants rather than the legacy generated-output count of 974. ",
    "The current explicit participant-event scaffold and published BPM ",
    "completeness rules remain authoritative.")
}


## Final Output Structure Checks ##

expected_cbcl_output_columns <- c(
  "subjectkey",
  "family_id",
  "eventname",
  "analysis_group",
  "group",
  "interview_age",
  "age_in_years",
  "sex",
  "site_name",
  "rsfmri_c_ngd_meanmotion",
  connectivity_vars,
  cbcl_vars)

expected_dimensional_output_columns <- c(
  expected_cbcl_output_columns,
  bpm_vars)

if (!identical(names(cbcl_data), expected_cbcl_output_columns)) {
  stop("The CBCL output columns are not in the expected order.")
}

if (!identical(
  names(dimensional_analysis_data),
  expected_dimensional_output_columns)) {
  stop("The combined dimensional output columns are not in the expected order.")
}


## Output ##

# All data construction and validation steps have passed. Replace only the two
# existing primary dimensional data files.
write_csv_safely(cbcl_data, cbcl_out_path)
write_csv_safely(dimensional_analysis_data, dimensional_out_path)

message("Primary dimensional data wrangling completed successfully.")
