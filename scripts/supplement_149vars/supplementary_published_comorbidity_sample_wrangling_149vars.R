## Setup ##

# Load packages for reading, wrangling, validation, and output
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

options(digits = 8, scipen = 999)


## Paths ##

# Final participant, assessment, diagnostic-group, and resampled-HC scaffold used for the published common-comorbidity sensitivity analysis
published_scaffold_path <- paste0(
  "./data_processed/supplement/",
  "comorbidity_hc_resampled_merged_groups.csv")

# Historical resampled HC file used to validate the HC portion of the scaffold
published_hc_path <- paste0(
  "./data_processed/supplement/",
  "supplementary_healthy_control_group_resampled.csv")

# Full corrected rs-fMRI source containing all available subject-event rows
corrected_imaging_all_rows_path <- paste0(
  "./data_processed/main_analysis/",
  "corrected_rsfmri_data_all_rows.csv")

# Primary 149-variable subset, used only to validate the full corrected source
corrected_imaging_validation_path <- paste0(
  "./data_processed/main_analysis/",
  "subset_qcd_imaging_data_149vars.csv")

# Direct rs-fMRI QC and family ID sources
qc_source_path <- "./data_raw/abcd_imgincl01.csv"
family_data_path <- "./data_raw/acspsw03.txt"

# Ordered manifest of the 20 corrected outcomes carried forward
metric_manifest_path <- paste0(
  "./results/results_149_vars/",
  "repeated_measures_dv_manifest_149vars.csv")

# Output directory
output_dir <- "./data_processed/supplement_149vars"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data_out_path <- file.path(
  output_dir,
  "comorbidity_connectivity_analysis_data_published_sample_149vars.csv")

sample_validation_out_path <- file.path(
  output_dir,
  "comorbidity_published_sample_validation_summary.csv")

group_event_counts_out_path <- file.path(
  output_dir,
  "comorbidity_published_sample_group_event_counts.csv")

site_counts_out_path <- file.path(
  output_dir,
  "comorbidity_published_sample_site_counts.csv")

imaging_coverage_out_path <- file.path(
  output_dir,
  "comorbidity_published_sample_imaging_coverage.csv")

source_comparison_out_path <- file.path(
  output_dir,
  "comorbidity_corrected_imaging_source_comparison.csv")


## Constants ##

analysis_events <- c(
  "baseline_year_1_arm_1",
  "2_year_follow_up_y_arm_1")

comorbidity_group_levels <- c(
  "HC",
  "GAD_Only",
  "GAD_with_comorbidities",
  "MDD_Only",
  "Separation_Anxiety_Only",
  "Social_Anxiety_Only")

# Locks this analysis to the  audited historical files
expected_total_n <- 3688L

expected_group_counts <- tibble::tribble(
  ~comorbidity_group, ~expected_n,
  "HC", 3364L,
  "GAD_Only", 123L,
  "GAD_with_comorbidities", 38L,
  "MDD_Only", 19L,
  "Separation_Anxiety_Only", 33L,
  "Social_Anxiety_Only", 111L)

expected_group_event_counts <- tibble::tribble(
  ~comorbidity_group, ~eventname, ~expected_n,
  "HC", "baseline_year_1_arm_1", 2202L,
  "HC", "2_year_follow_up_y_arm_1", 1162L,
  "GAD_Only", "baseline_year_1_arm_1", 74L,
  "GAD_Only", "2_year_follow_up_y_arm_1", 49L,
  "GAD_with_comorbidities", "baseline_year_1_arm_1", 24L,
  "GAD_with_comorbidities", "2_year_follow_up_y_arm_1", 14L,
  "MDD_Only", "baseline_year_1_arm_1", 10L,
  "MDD_Only", "2_year_follow_up_y_arm_1", 9L,
  "Separation_Anxiety_Only", "baseline_year_1_arm_1", 25L,
  "Separation_Anxiety_Only", "2_year_follow_up_y_arm_1", 8L,
  "Social_Anxiety_Only", "baseline_year_1_arm_1", 73L,
  "Social_Anxiety_Only", "2_year_follow_up_y_arm_1", 38L)

expected_event_counts <- tibble::tribble(
  ~eventname, ~expected_n,
  "baseline_year_1_arm_1", 2408L,
  "2_year_follow_up_y_arm_1", 1280L)

numeric_comparison_tolerance <- 1e-12


## Helper Functions ##

clean_character <- function(x) {
  dplyr::na_if(trimws(as.character(x)), "")
}

assert_columns <- function(data, required, dataset_name) {
  missing_columns <- setdiff(required, names(data))

  if (length(missing_columns) > 0L) {
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
    dplyr::filter(n_rows > 1L)

  if (nrow(duplicate_keys) > 0L) {
    stop(
      dataset_name,
      " contains duplicated subjectkey-eventname keys.")
  }
}

standardize_subject_id <- function(
    data,
    dataset_name,
    id_candidates = c("subjectkey", "src_subject_id"),
    compare_existing_subjectkey = TRUE) {
  
  matching_ids <- intersect(
    id_candidates,
    names(data))
  
  if (length(matching_ids) == 0L) {
    stop(
      dataset_name,
      " contains none of the supported participant-ID columns: ",
      paste(id_candidates, collapse = ", "),
      ".")
  }
  
  id_column <- matching_ids[[1]]
  
  standardized_id <- clean_character(
    data[[id_column]])
  
  if (
    isTRUE(compare_existing_subjectkey) &&
    "subjectkey" %in% names(data) &&
    id_column != "subjectkey") {
    
    existing_subjectkey <- clean_character(
      data$subjectkey)
    
    both_observed <-
      !is.na(existing_subjectkey) &
      !is.na(standardized_id)
    
    if (
      any(
        existing_subjectkey[both_observed] !=
        standardized_id[both_observed])) {
      stop(
        dataset_name,
        " contains disagreeing participant-ID columns.")
    }
  }
  
  data$subjectkey <- standardized_id
  
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

assert_exact_counts <- function(
    data,
    expected_counts,
    grouping_variables,
    dataset_name) {

  observed_counts <- data %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(observed_n = dplyr::n(), .groups = "drop")

  count_check <- expected_counts %>%
    dplyr::full_join(observed_counts, by = grouping_variables) %>%
    dplyr::mutate(
      expected_n = tidyr::replace_na(expected_n, 0L),
      observed_n = tidyr::replace_na(observed_n, 0L),
      matches = expected_n == observed_n)

  if (any(!count_check$matches)) {
    print(count_check)
    stop(dataset_name, " does not have the exact expected counts.")
  }

  invisible(count_check)
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


## Published Sample Scaffold ##

published_scaffold <- readr::read_csv(
  published_scaffold_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Published comorbidity scaffold",
    id_candidates = c("src_subject_id", "subjectkey")) %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    interview_age = suppressWarnings(as.numeric(interview_age)),
    sex = clean_character(sex),
    site_name = clean_character(site_name),
    comorbidity_group = clean_character(comorbidity_group),
    scaffold_row_order = dplyr::row_number()) %>%
  dplyr::filter(eventname %in% analysis_events)

assert_columns(
  published_scaffold,
  c(
    "subjectkey",
    "eventname",
    "interview_age",
    "sex",
    "comorbidity_group"),
  "Published comorbidity scaffold")

assert_unique_keys(published_scaffold, "Published comorbidity scaffold")

if (
  nrow(published_scaffold) != expected_total_n ||
    dplyr::n_distinct(published_scaffold$subjectkey) != expected_total_n) {
  stop(
    "The published comorbidity scaffold must contain exactly ",
    expected_total_n,
    " unique participants.")
}

unexpected_groups <- setdiff(
  unique(published_scaffold$comorbidity_group),
  comorbidity_group_levels)

if (length(unexpected_groups) > 0L) {
  stop(
    "Unexpected comorbidity groups in the published scaffold: ",
    paste(unexpected_groups, collapse = ", "),
    ".")
}

required_covariates <- c(
  "eventname",
  "interview_age",
  "sex",
  "site_name",
  "comorbidity_group")

covariate_missingness <- vapply(
  required_covariates,
  function(variable) {
    values <- published_scaffold[[variable]]

    sum(
      is.na(values) |
        (is.character(values) & trimws(values) == ""))
  },
  integer(1))

if (any(covariate_missingness > 0L)) {
  print(covariate_missingness)
  stop("The published scaffold contains incomplete required covariates.")
}

assert_exact_counts(
  published_scaffold,
  expected_group_counts,
  "comorbidity_group",
  "Published group distribution")

assert_exact_counts(
  published_scaffold,
  expected_group_event_counts,
  c("comorbidity_group", "eventname"),
  "Published group-event distribution")

assert_exact_counts(
  published_scaffold,
  expected_event_counts,
  "eventname",
  "Published event distribution")


## Historical Resampled HC Validation ##

published_hc <- readr::read_csv(
  published_hc_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Published resampled HC file",
    id_candidates = c("src_subject_id", "subjectkey")) %>%
  dplyr::mutate(eventname = clean_character(eventname)) %>%
  dplyr::select(subjectkey, eventname) %>%
  dplyr::distinct()

published_scaffold_hc <- published_scaffold %>%
  dplyr::filter(comorbidity_group == "HC") %>%
  dplyr::select(subjectkey, eventname)

if (
  nrow(published_hc) != 3364L ||
    nrow(
      dplyr::anti_join(
        published_hc,
        published_scaffold_hc,
        by = c("subjectkey", "eventname"))) > 0L ||
    nrow(
      dplyr::anti_join(
        published_scaffold_hc,
        published_hc,
        by = c("subjectkey", "eventname"))) > 0L) {
  stop(
    "The HC rows in the published scaffold do not exactly match the ",
    "historical resampled HC file.")
}


## Full Corrected Imaging Source ##

corrected_imaging_all_rows <- readr::read_csv(
  corrected_imaging_all_rows_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Full corrected rs-fMRI source")

assert_columns(
  corrected_imaging_all_rows,
  c(
    "subjectkey",
    "eventname",
    "rsfmri_c_ngd_meanmotion",
    connectivity_vars),
  "Full corrected rs-fMRI source")

corrected_imaging_all_rows <- corrected_imaging_all_rows %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    dplyr::across(
      dplyr::all_of(connectivity_vars),
      ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::filter(eventname %in% analysis_events) %>%
  dplyr::select(
    subjectkey,
    eventname,
    rsfmri_c_ngd_meanmotion,
    dplyr::all_of(connectivity_vars)) %>%
  dplyr::mutate(
    corrected_imaging_row_present = TRUE)

assert_unique_keys(
  corrected_imaging_all_rows,
  "Full corrected rs-fMRI source")


## Validate the Full Corrected Source ##

source_comparison <- tibble::tibble(
  variable = character(),
  n_overlapping_nonmissing = integer(),
  n_mismatched = integer(),
  maximum_absolute_difference = double())

if (file.exists(corrected_imaging_validation_path)) {

  corrected_imaging_validation <- readr::read_csv(
    corrected_imaging_validation_path,
    show_col_types = FALSE) %>%
    standardize_subject_id(
      dataset_name = "Primary corrected 149-variable validation source") %>%
    dplyr::mutate(
      eventname = clean_character(eventname),
      rsfmri_c_ngd_meanmotion =
        suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
      dplyr::across(
        dplyr::all_of(connectivity_vars),
        ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::filter(eventname %in% analysis_events) %>%
    dplyr::select(
      subjectkey,
      eventname,
      rsfmri_c_ngd_meanmotion,
      dplyr::all_of(connectivity_vars))

  assert_unique_keys(
    corrected_imaging_validation,
    "Primary corrected 149-variable validation source")

  overlapping_sources <- corrected_imaging_all_rows %>%
    dplyr::inner_join(
      corrected_imaging_validation,
      by = c("subjectkey", "eventname"),
      suffix = c("_full", "_validation"))

  if (nrow(overlapping_sources) == 0L) {
    stop(
      "The full corrected source and validation subset have no overlapping rows.")
  }

  comparison_vars <- c(
    "rsfmri_c_ngd_meanmotion",
    connectivity_vars)

  source_comparison <- dplyr::bind_rows(
    lapply(
      comparison_vars,
      function(variable) {
        full_values <-
          overlapping_sources[[paste0(variable, "_full")]]

        validation_values <-
          overlapping_sources[[paste0(variable, "_validation")]]

        observed_pair <-
          !is.na(full_values) &
          !is.na(validation_values)

        differences <- abs(
          full_values[observed_pair] -
            validation_values[observed_pair])

        tibble::tibble(
          variable = variable,
          n_overlapping_nonmissing = sum(observed_pair),
          n_mismatched =
            sum(differences > numeric_comparison_tolerance),
          maximum_absolute_difference =
            if (length(differences) == 0L) {
              NA_real_
            } else {
              max(differences)
            })
      }))

  if (any(source_comparison$n_mismatched > 0L)) {
    print(
      source_comparison %>%
        dplyr::filter(n_mismatched > 0L))

    stop(
      "The full corrected imaging source does not match the validated ",
      "149-variable subset on overlapping rows.")
  }

} else {

  source_comparison <- tibble::tibble(
    variable = "Validation subset unavailable",
    n_overlapping_nonmissing = NA_integer_,
    n_mismatched = NA_integer_,
    maximum_absolute_difference = NA_real_)

  warning(
    "The primary corrected 149-variable validation source was not found.")
}


## Direct rs-fMRI QC ##

direct_qc <- readr::read_csv(
  qc_source_path,
  show_col_types = FALSE) %>%
  dplyr::slice(-1) %>%
  standardize_subject_id(
    dataset_name = "Direct rs-fMRI QC source",
    id_candidates = "src_subject_id",
    compare_existing_subjectkey = FALSE) %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include))) %>%
  dplyr::filter(eventname %in% analysis_events) %>%
  dplyr::select(
    subjectkey,
    eventname,
    imgincl_rsfmri_include)

assert_unique_keys(direct_qc, "Direct rs-fMRI QC source")


## Family ID Data ##

family_data <- readr::read_tsv(
  family_data_path,
  show_col_types = FALSE) %>%
  standardize_subject_id(
    dataset_name = "Family ID source") %>%
  dplyr::mutate(
    eventname = clean_character(eventname),
    family_id = clean_character(rel_family_id)) %>%
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::select(subjectkey, family_id) %>%
  dplyr::distinct()

family_duplicates <- family_data %>%
  dplyr::count(subjectkey, name = "n_rows") %>%
  dplyr::filter(n_rows > 1L)

if (nrow(family_duplicates) > 0L) {
  stop(
    "The family ID source contains multiple baseline family IDs for ",
    nrow(family_duplicates),
    " participants.")
}


## Attach Corrected Imaging Without Changing the Published Sample ##

analysis_data <- published_scaffold %>%
  dplyr::left_join(
    corrected_imaging_all_rows,
    by = c("subjectkey", "eventname")) %>%
  dplyr::left_join(
    direct_qc,
    by = c("subjectkey", "eventname")) %>%
  dplyr::left_join(
    family_data,
    by = "subjectkey") %>%
  dplyr::mutate(
    age_in_years = interview_age / 12)

if (
  nrow(analysis_data) != expected_total_n ||
    dplyr::n_distinct(analysis_data$subjectkey) != expected_total_n) {
  stop("Joining corrected data changed the number of rows or participants.")
}

if (
  !identical(analysis_data$subjectkey, published_scaffold$subjectkey) ||
    !identical(analysis_data$eventname, published_scaffold$eventname) ||
    !identical(
      analysis_data$comorbidity_group,
      published_scaffold$comorbidity_group) ||
    !identical(
      analysis_data$scaffold_row_order,
      published_scaffold$scaffold_row_order)) {
  stop(
    "Joining corrected data changed the order, event, or diagnostic-group ",
    "assignment of the published scaffold.")
}

complete_all_20 <- apply(
  !is.na(analysis_data[connectivity_vars]),
  1,
  all)

analysis_missingness <- tibble::tibble(
  criterion = c(
    "No full corrected imaging match",
    "Missing mean motion",
    "Missing at least one corrected outcome",
    "No direct QC match",
    "Direct rs-fMRI QC is not 1",
    "Missing family ID"),
  n_rows = c(
    sum(!analysis_data$corrected_imaging_row_present),
    sum(is.na(analysis_data$rsfmri_c_ngd_meanmotion)),
    sum(!complete_all_20),
    sum(is.na(analysis_data$imgincl_rsfmri_include)),
    sum(
      is.na(analysis_data$imgincl_rsfmri_include) |
        analysis_data$imgincl_rsfmri_include != 1),
    sum(is.na(analysis_data$family_id))))

imaging_coverage <- analysis_data %>%
  dplyr::mutate(
    complete_all_20_outcomes = complete_all_20,
    direct_rsfmri_qc_pass = imgincl_rsfmri_include == 1) %>%
  dplyr::group_by(comorbidity_group, eventname) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_without_full_corrected_imaging_match =
      sum(!corrected_imaging_row_present),
    n_without_mean_motion =
      sum(is.na(rsfmri_c_ngd_meanmotion)),
    n_without_complete_all_20_outcomes =
      sum(!complete_all_20_outcomes),
    n_without_direct_rsfmri_qc_pass =
      sum(
        is.na(direct_rsfmri_qc_pass) |
          !direct_rsfmri_qc_pass),
    n_without_family_id =
      sum(is.na(family_id)),
    .groups = "drop")

group_event_counts <- analysis_data %>%
  dplyr::count(
    comorbidity_group,
    eventname,
    name = "n_participants") %>%
  dplyr::arrange(
    factor(comorbidity_group, levels = comorbidity_group_levels),
    factor(eventname, levels = analysis_events))

site_counts <- analysis_data %>%
  dplyr::count(site_name, name = "n_participants") %>%
  dplyr::arrange(site_name)

sample_validation <- tibble::tibble(
  n_rows = nrow(analysis_data),
  n_unique_participants =
    dplyr::n_distinct(analysis_data$subjectkey),
  n_unique_subject_event_keys =
    dplyr::n_distinct(
      paste(
        analysis_data$subjectkey,
        analysis_data$eventname,
        sep = "||")),
  n_hc =
    sum(analysis_data$comorbidity_group == "HC"),
  n_clinical =
    sum(analysis_data$comorbidity_group != "HC"),
  n_baseline =
    sum(analysis_data$eventname == "baseline_year_1_arm_1"),
  n_followup =
    sum(analysis_data$eventname == "2_year_follow_up_y_arm_1"),
  n_sites =
    dplyr::n_distinct(analysis_data$site_name),
  n_without_full_corrected_imaging_match =
    analysis_missingness$n_rows[
      analysis_missingness$criterion ==
        "No full corrected imaging match"],
  n_without_complete_all_20_outcomes =
    analysis_missingness$n_rows[
      analysis_missingness$criterion ==
        "Missing at least one corrected outcome"],
  n_without_direct_rsfmri_qc_pass =
    analysis_missingness$n_rows[
      analysis_missingness$criterion ==
        "Direct rs-fMRI QC is not 1"],
  n_without_family_id =
    analysis_missingness$n_rows[
      analysis_missingness$criterion ==
        "Missing family ID"])

print(sample_validation)
print(analysis_missingness)
print(group_event_counts)
print(site_counts)


## Write Aggregate Audit Outputs Before the Final Eligibility Check ##

write_csv_safely(
  sample_validation,
  sample_validation_out_path)

write_csv_safely(
  group_event_counts,
  group_event_counts_out_path)

write_csv_safely(
  site_counts,
  site_counts_out_path)

write_csv_safely(
  imaging_coverage,
  imaging_coverage_out_path)

write_csv_safely(
  source_comparison,
  source_comparison_out_path)


## Final Eligibility Check and Analysis Dataset Output ##

if (any(analysis_missingness$n_rows > 0L)) {
  stop(
    "The exact published scaffold could not be fully populated with corrected ",
    "imaging, direct QC, and family data. Aggregate audit files were written ",
    "to ",
    output_dir,
    ". No participant was dropped and no analysis dataset was written.")
}

analysis_data_for_output <- analysis_data %>%
  dplyr::arrange(scaffold_row_order) %>%
  dplyr::select(
    subjectkey,
    family_id,
    eventname,
    interview_age,
    age_in_years,
    sex,
    site_name,
    group,
    comorbidity_group,
    broad_clinical_group,
    dplyr::everything(),
    -scaffold_row_order,
    -corrected_imaging_row_present)

write_csv_safely(
  analysis_data_for_output,
  analysis_data_out_path)

message(
  "Exact published common-comorbidity sample retained: ",
  expected_total_n,
  " participants. Corrected imaging values were attached without resampling, ",
  "reassignment, event substitution, or participant exclusion.")
