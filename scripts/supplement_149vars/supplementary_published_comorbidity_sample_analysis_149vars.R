## Setup ##

# Load packages for mixed modeling, inference, multiplicity correction,
# validation, and output
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)

options(digits = 8, scipen = 999, contrasts = c("contr.treatment", "contr.poly"))

emmeans::emm_options(
  lmer.df = "satterthwaite",
  lmerTest.limit = 5000,
  disable.lmerTest = FALSE,
  disable.pbkrtest = TRUE)


## Paths ##

analysis_data_path <- paste0(
  "./data_processed/supplement_149vars/",
  "comorbidity_connectivity_analysis_data_published_sample_149vars.csv")

published_scaffold_path <- paste0(
  "./data_processed/supplement/",
  "comorbidity_hc_resampled_merged_groups.csv")

metric_manifest_path <- paste0(
  "./results/results_149_vars/",
  "repeated_measures_dv_manifest_149vars.csv")

results_dir <- paste0(
  "./results/results_149_vars/supplement/",
  "comorbidity_published_sample")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

omnibus_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_omnibus_results.csv")

estimated_marginal_means_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_estimated_marginal_means.csv")

pairwise_posthocs_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_pairwise_posthocs.csv")

model_diagnostics_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_model_diagnostics.csv")

sample_characteristics_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_sample_characteristics.csv")

group_event_counts_out_path <- file.path(
  results_dir,
  "comorbidity_connectivity_group_event_counts.csv")


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

expected_published_scaffold_md5 <-
  "2f83ba4c6079abd5f4eafa3b0fe6bca9"

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


## Connectivity Outcome Manifest and Labels ##

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
  stop("The connectivity manifest must contain exactly 20 unique outcomes.")
}

if (!identical(metric_manifest$analysis_order, seq_len(20L))) {
  stop("The connectivity manifest analysis_order must be exactly 1 through 20.")
}

metric_labels <- tibble::tribble(
  ~column_name, ~DV_Summary,
  "rsfmri_cor_ngd_au_scs_thplh", "AN - Left Thalamus",
  "rsfmri_cor_ngd_vta_scs_plrh", "VAN - Right Pallidum",
  "rsfmri_cor_ngd_smh_scs_plrh", "SMHN - Right Pallidum",
  "rsfmri_cor_ngd_au_scs_thprh", "AN - Right Thalamus",
  "rsfmri_cor_ngd_au_scs_cdelh", "AN - Left Caudate",
  "rsfmri_cor_ngd_au_scs_ptlh", "AN - Left Putamen",
  "rsfmri_c_ngd_ad_ngd_sa", "AN - SN",
  "rsfmri_cor_ngd_au_scs_pllh", "AN - Left Pallidum",
  "rsfmri_cor_ngd_smh_scs_pllh", "SMHN - Left Pallidum",
  "rsfmri_cor_ngd_au_scs_plrh", "AN - Right Pallidum",
  "rsfmri_c_ngd_vta_ngd_vta", "Within-VAN",
  "rsfmri_cor_ngd_smh_scs_thplh", "SMHN - Left Thalamus",
  "rsfmri_cor_ngd_smh_scs_thprh", "SMHN - Right Thalamus",
  "rsfmri_cor_ngd_vta_scs_ptrh", "VAN - Right Putamen",
  "rsfmri_cor_ngd_vta_scs_thplh", "VAN - Left Thalamus",
  "rsfmri_cor_ngd_vta_scs_ptlh", "VAN - Left Putamen",
  "rsfmri_cor_ngd_smh_scs_cderh", "SMHN - Right Caudate",
  "rsfmri_cor_ngd_au_scs_cderh", "AN - Right Caudate",
  "rsfmri_cor_ngd_au_scs_ptrh", "AN - Right Putamen",
  "rsfmri_cor_ngd_smh_scs_cdelh", "SMHN - Left Caudate")

metric_manifest <- metric_manifest %>%
  dplyr::left_join(metric_labels, by = "column_name")

if (any(is.na(metric_manifest$DV_Summary))) {
  stop(
    "Missing human-readable labels for: ",
    paste(
      metric_manifest$column_name[is.na(metric_manifest$DV_Summary)],
      collapse = ", "),
    ".")
}

connectivity_vars <- metric_manifest$column_name


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
    stop(dataset_name, " contains duplicated subjectkey-eventname rows.")
  }
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

normalize_term <- function(x) {
  x %>%
    stringr::str_replace_all("`", "") %>%
    stringr::str_replace_all("\\s", "")
}

get_numeric_column <- function(data, candidates) {
  matching_columns <- intersect(candidates, names(data))

  if (length(matching_columns) == 0L) {
    return(rep(NA_real_, nrow(data)))
  }

  suppressWarnings(as.numeric(data[[matching_columns[[1]]]]))
}

fit_lmer_safe <- function(formula, data) {
  model_warnings <- character()
  model_error <- NA_character_

  model <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        formula = formula,
        data = data,
        REML = TRUE,
        na.action = na.omit),
      warning = function(warning_condition) {
        model_warnings <<- c(
          model_warnings,
          conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }),
    error = function(error_condition) {
      model_error <<- conditionMessage(error_condition)
      NULL
    })

  list(
    model = model,
    warnings = if (length(model_warnings) == 0L) {
      NA_character_
    } else {
      paste(unique(model_warnings), collapse = " | ")
    },
    error = model_error)
}

extract_anova_term <- function(model, term, type = "II") {
  empty_result <- tibble::tibble(
    f_value = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    p_value = NA_real_)

  if (is.null(model)) {
    return(empty_result)
  }

  anova_object <- tryCatch(
    car::Anova(
      model,
      type = type,
      test.statistic = "F"),
    error = function(error_condition) NULL)

  if (is.null(anova_object)) {
    return(empty_result)
  }

  anova_table <- tibble::as_tibble(
    as.data.frame(anova_object),
    rownames = "term") %>%
    dplyr::mutate(normalized_term = normalize_term(term))

  matched_row <- which(
    anova_table$normalized_term == normalize_term(term))

  if (length(matched_row) != 1L) {
    return(empty_result)
  }

  row <- anova_table[matched_row, , drop = FALSE]

  tibble::tibble(
    f_value = get_numeric_column(
      row,
      c("F", "F value", "F.value"))[[1]],
    numerator_df = get_numeric_column(
      row,
      c("Df", "num Df", "NumDF"))[[1]],
    denominator_df = get_numeric_column(
      row,
      c("Df.res", "den Df", "DenDF"))[[1]],
    p_value = get_numeric_column(
      row,
      c("Pr(>F)", "Pr..F."))[[1]])
}

extract_coefficient <- function(model, term) {
  empty_result <- tibble::tibble(
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    t_value = NA_real_,
    p_value = NA_real_)

  if (is.null(model)) {
    return(empty_result)
  }

  coefficient_table <- as.data.frame(summary(model)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  matched_row <- which(coefficient_table$term == term)

  if (length(matched_row) != 1L) {
    return(empty_result)
  }

  row <- coefficient_table[matched_row, , drop = FALSE]

  tibble::tibble(
    estimate = get_numeric_column(row, "Estimate")[[1]],
    std_error = get_numeric_column(
      row,
      c("Std. Error", "Std..Error"))[[1]],
    df = get_numeric_column(row, "df")[[1]],
    t_value = get_numeric_column(
      row,
      c("t value", "t.value"))[[1]],
    p_value = get_numeric_column(
      row,
      c("Pr(>|t|)", "Pr...t.."))[[1]])
}

extract_diagnostics <- function(fit_result, data) {
  model <- fit_result$model

  if (is.null(model)) {
    return(
      tibble::tibble(
        n_observations = NA_integer_,
        n_subjects = dplyr::n_distinct(data$subjectkey),
        singular = NA,
        optimizer_code = NA_character_,
        convergence_messages = NA_character_,
        warnings = fit_result$warnings,
        error = fit_result$error,
        AIC = NA_real_,
        BIC = NA_real_,
        logLik = NA_real_))
  }

  convergence_messages <- model@optinfo$conv$lme4$messages
  convergence_messages <- if (is.null(convergence_messages)) {
    NA_character_
  } else {
    paste(convergence_messages, collapse = " | ")
  }

  optimizer_code <- model@optinfo$conv$opt
  optimizer_code <- if (is.null(optimizer_code)) {
    NA_character_
  } else {
    paste(optimizer_code, collapse = " | ")
  }

  tibble::tibble(
    n_observations = stats::nobs(model),
    n_subjects = dplyr::n_distinct(data$subjectkey),
    singular = lme4::isSingular(model, tol = 1e-4),
    optimizer_code = optimizer_code,
    convergence_messages = convergence_messages,
    warnings = fit_result$warnings,
    error = fit_result$error,
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)))
}

prefix_columns <- function(data, prefix) {
  names(data) <- paste0(prefix, names(data))
  data
}

has_text <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
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
    stop("Failed to replace result file: ", path)
  }

  message("Result written: ", path)
}


## Read and Validate the Exact Published Analysis Sample ##

observed_scaffold_md5 <- unname(tools::md5sum(published_scaffold_path))

if (!identical(observed_scaffold_md5, expected_published_scaffold_md5)) {
  stop("The historical published scaffold no longer matches the audited file.")
}

analysis_data <- readr::read_csv(
  analysis_data_path,
  show_col_types = FALSE)

required_variables <- c(
  "subjectkey",
  "family_id",
  "eventname",
  "comorbidity_group",
  "interview_age",
  "age_in_years",
  "sex",
  "site_name",
  "imgincl_rsfmri_include",
  "rsfmri_c_ngd_meanmotion",
  connectivity_vars)

assert_columns(
  analysis_data,
  required_variables,
  "Published-sample comorbidity analysis dataset")

analysis_data <- analysis_data %>%
  dplyr::mutate(
    subjectkey = clean_character(subjectkey),
    family_id = clean_character(family_id),
    eventname = clean_character(eventname),
    comorbidity_group = clean_character(comorbidity_group),
    interview_age = suppressWarnings(as.numeric(interview_age)),
    age_in_years = suppressWarnings(as.numeric(age_in_years)),
    sex = clean_character(sex),
    site_name = clean_character(site_name),
    imgincl_rsfmri_include =
      suppressWarnings(as.numeric(imgincl_rsfmri_include)),
    rsfmri_c_ngd_meanmotion =
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    dplyr::across(
      dplyr::all_of(connectivity_vars),
      ~ suppressWarnings(as.numeric(.x))))

assert_unique_keys(
  analysis_data,
  "Published-sample comorbidity analysis dataset")

if (
  nrow(analysis_data) != expected_total_n ||
    dplyr::n_distinct(analysis_data$subjectkey) != expected_total_n) {
  stop(
    "The analysis dataset must contain exactly ",
    expected_total_n,
    " unique participants.")
}

assert_exact_counts(
  analysis_data,
  expected_group_counts,
  "comorbidity_group",
  "Comorbidity analysis group distribution")

assert_exact_counts(
  analysis_data,
  expected_group_event_counts,
  c("comorbidity_group", "eventname"),
  "Comorbidity analysis group-event distribution")

missing_counts <- vapply(
  required_variables,
  function(variable) {
    values <- analysis_data[[variable]]

    sum(
      is.na(values) |
        (is.character(values) & trimws(values) == ""))
  },
  integer(1))

if (any(missing_counts > 0L)) {
  print(missing_counts[missing_counts > 0L])
  stop("The exact published analysis sample contains incomplete required data.")
}

if (any(analysis_data$imgincl_rsfmri_include != 1)) {
  stop(
    "At least one published-sample row does not have direct ",
    "imgincl_rsfmri_include == 1.")
}

analysis_data <- analysis_data %>%
  dplyr::mutate(
    subjectkey = factor(subjectkey),
    family_id = factor(family_id),
    eventname = factor(eventname, levels = analysis_events),
    comorbidity_group = factor(
      comorbidity_group,
      levels = comorbidity_group_levels),
    sex = factor(sex),
    site_name = factor(site_name))

message(
  "Exact published common-comorbidity analysis sample validated: ",
  nrow(analysis_data),
  " participants.")

print(
  analysis_data %>%
    dplyr::count(
      comorbidity_group,
      eventname,
      name = "n_participants"))


## Model Function ##

fit_comorbidity_outcome <- function(
    outcome,
    analysis_order,
    outcome_label) {

  model_formula <- stats::as.formula(
    paste0(
      outcome,
      " ~ comorbidity_group + rsfmri_c_ngd_meanmotion + ",
      "sex + eventname + (1 | site_name) + (1 | family_id)"))

  fit_result <- fit_lmer_safe(model_formula, analysis_data)
  model <- fit_result$model

  group_test <- extract_anova_term(
    model,
    "comorbidity_group",
    type = "II")

  motion_result <- extract_coefficient(
    model,
    "rsfmri_c_ngd_meanmotion")

  sex_test <- extract_anova_term(model, "sex", type = "II")
  event_test <- extract_anova_term(model, "eventname", type = "II")
  diagnostics <- extract_diagnostics(fit_result, analysis_data)

  base_columns <- tibble::tibble(
    analysis_order = analysis_order,
    column_name = outcome,
    DV_Summary = outcome_label,
    model_formula = paste(deparse(model_formula), collapse = " "))

  if (is.null(model)) {
    return(
      list(
        omnibus = dplyr::bind_cols(
          base_columns,
          prefix_columns(group_test, "comorbidity_group_"),
          prefix_columns(motion_result, "motion_"),
          prefix_columns(sex_test, "sex_"),
          prefix_columns(event_test, "event_")),
        emmeans = NULL,
        posthocs = NULL,
        diagnostics = dplyr::bind_cols(base_columns, diagnostics)))
  }

  emm <- emmeans::emmeans(model, ~ comorbidity_group)

  emm_summary <- summary(
    emm,
    infer = c(TRUE, FALSE)) %>%
    as.data.frame() %>%
    dplyr::transmute(
      comorbidity_group = as.character(comorbidity_group),
      emmean = suppressWarnings(as.numeric(emmean)),
      std_error = suppressWarnings(as.numeric(SE)),
      df = suppressWarnings(as.numeric(df)),
      lower_ci = suppressWarnings(as.numeric(lower.CL)),
      upper_ci = suppressWarnings(as.numeric(upper.CL)))

  pairwise_summary <- emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = "none") %>%
    summary(
      infer = c(TRUE, TRUE),
      adjust = "none") %>%
    as.data.frame() %>%
    dplyr::transmute(
      contrast = as.character(contrast),
      estimate = suppressWarnings(as.numeric(estimate)),
      std_error = suppressWarnings(as.numeric(SE)),
      df = suppressWarnings(as.numeric(df)),
      lower_ci = suppressWarnings(as.numeric(lower.CL)),
      upper_ci = suppressWarnings(as.numeric(upper.CL)),
      t_value = suppressWarnings(as.numeric(t.ratio)),
      p_value = suppressWarnings(as.numeric(p.value))) %>%
    dplyr::mutate(
      p_value_fdr_within_outcome =
        stats::p.adjust(p_value, method = "fdr"))

  list(
    omnibus = dplyr::bind_cols(
      base_columns,
      prefix_columns(group_test, "comorbidity_group_"),
      prefix_columns(motion_result, "motion_"),
      prefix_columns(sex_test, "sex_"),
      prefix_columns(event_test, "event_")),
    emmeans = dplyr::bind_cols(
      base_columns[rep(1L, nrow(emm_summary)), ],
      emm_summary),
    posthocs = dplyr::bind_cols(
      base_columns[rep(1L, nrow(pairwise_summary)), ],
      pairwise_summary),
    diagnostics = dplyr::bind_cols(base_columns, diagnostics))
}


## Run the 20 Corrected Outcomes ##

model_grid <- metric_manifest %>%
  dplyr::select(
    analysis_order,
    column_name,
    DV_Summary)

model_results <- purrr::pmap(
  model_grid,
  function(
      analysis_order,
      column_name,
      DV_Summary) {

    message(
      "Modeling ",
      DV_Summary,
      " by common-comorbidity group.")

    fit_comorbidity_outcome(
      outcome = column_name,
      analysis_order = analysis_order,
      outcome_label = DV_Summary)
  })

omnibus_results <- purrr::map_dfr(
  model_results,
  "omnibus") %>%
  dplyr::arrange(analysis_order) %>%
  dplyr::mutate(
    comorbidity_group_p_value_fdr_across_20 =
      stats::p.adjust(
        comorbidity_group_p_value,
        method = "fdr"),
    comorbidity_group_nominal_significant =
      comorbidity_group_p_value < 0.05,
    comorbidity_group_fdr_significant =
      comorbidity_group_p_value_fdr_across_20 < 0.05)

estimated_marginal_means <- purrr::map_dfr(
  model_results,
  "emmeans") %>%
  dplyr::arrange(
    analysis_order,
    factor(
      comorbidity_group,
      levels = comorbidity_group_levels))

pairwise_posthocs <- purrr::map_dfr(
  model_results,
  "posthocs") %>%
  dplyr::left_join(
    omnibus_results %>%
      dplyr::select(
        column_name,
        comorbidity_group_p_value,
        comorbidity_group_p_value_fdr_across_20,
        comorbidity_group_nominal_significant,
        comorbidity_group_fdr_significant),
    by = "column_name") %>%
  dplyr::arrange(analysis_order, p_value)

model_diagnostics <- purrr::map_dfr(
  model_results,
  "diagnostics") %>%
  dplyr::arrange(analysis_order)


## Validate Model Outputs ##

if (nrow(omnibus_results) != 20L) {
  stop("The omnibus result table must contain exactly 20 rows.")
}

expected_emmean_rows <-
  20L * length(comorbidity_group_levels)

if (nrow(estimated_marginal_means) != expected_emmean_rows) {
  stop(
    "The estimated-marginal-means table contains ",
    nrow(estimated_marginal_means),
    " rows; expected ",
    expected_emmean_rows,
    ".")
}

expected_pairwise_rows <-
  20L * choose(length(comorbidity_group_levels), 2L)

if (nrow(pairwise_posthocs) != expected_pairwise_rows) {
  stop(
    "The pairwise post-hoc table contains ",
    nrow(pairwise_posthocs),
    " rows; expected ",
    expected_pairwise_rows,
    ".")
}

error_rows <- model_diagnostics %>%
  dplyr::filter(has_text(error))

convergence_rows <- model_diagnostics %>%
  dplyr::filter(has_text(convergence_messages))

optimizer_failure_rows <- model_diagnostics %>%
  dplyr::filter(
    has_text(optimizer_code) &
      trimws(as.character(optimizer_code)) != "0")

if (nrow(error_rows) > 0L) {
  print(error_rows)
  stop("One or more comorbidity models failed. No result files were written.")
}

if (nrow(convergence_rows) > 0L) {
  print(convergence_rows)
  stop(
    "One or more comorbidity models produced convergence messages. ",
    "No result files were written.")
}

if (nrow(optimizer_failure_rows) > 0L) {
  print(optimizer_failure_rows)
  stop(
    "One or more comorbidity models returned a nonzero optimizer code. ",
    "No result files were written.")
}

required_p_values <- c(
  "comorbidity_group_p_value",
  "motion_p_value",
  "sex_p_value",
  "event_p_value")

missing_inference <- vapply(
  required_p_values,
  function(variable) {
    sum(!is.finite(omnibus_results[[variable]]))
  },
  integer(1))

if (any(missing_inference > 0L)) {
  print(missing_inference[missing_inference > 0L])
  stop("Required model inference is missing. No result files were written.")
}

if (any(model_diagnostics$n_observations != expected_total_n)) {
  stop(
    "At least one model used fewer than the exact ",
    expected_total_n,
    " published-sample observations.")
}

message(
  "All 20 common-comorbidity models passed model-error, convergence, ",
  "optimizer, observation-count, and inferential checks.")

print(
  tibble::tibble(
    n_models = nrow(model_diagnostics),
    n_singular = sum(model_diagnostics$singular %in% TRUE),
    n_nominal_group_effects =
      sum(omnibus_results$comorbidity_group_nominal_significant),
    n_fdr_group_effects =
      sum(omnibus_results$comorbidity_group_fdr_significant)))


## Core Sample Characteristics ##

sample_characteristics <- analysis_data %>%
  dplyr::mutate(
    scanner_model = dplyr::case_when(
      stringr::str_starts(as.character(site_name), "S") ~ "Siemens",
      stringr::str_starts(as.character(site_name), "P") ~
        "Philips Healthcare",
      stringr::str_starts(as.character(site_name), "G") ~
        "GE Healthcare",
      TRUE ~ "Unknown")) %>%
  dplyr::group_by(comorbidity_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_age_years = mean(age_in_years),
    sd_age_years = stats::sd(age_in_years),
    n_female = sum(sex == "F"),
    percent_female = 100 * n_female / n,
    mean_framewise_displacement = mean(rsfmri_c_ngd_meanmotion),
    sd_framewise_displacement = stats::sd(rsfmri_c_ngd_meanmotion),
    n_baseline = sum(eventname == "baseline_year_1_arm_1"),
    n_followup = sum(eventname == "2_year_follow_up_y_arm_1"),
    n_siemens = sum(scanner_model == "Siemens"),
    n_philips = sum(scanner_model == "Philips Healthcare"),
    n_ge = sum(scanner_model == "GE Healthcare"),
    n_unknown_scanner = sum(scanner_model == "Unknown"),
    .groups = "drop")

group_event_counts <- analysis_data %>%
  dplyr::count(
    comorbidity_group,
    eventname,
    name = "n_participants")


## Output ##

write_csv_safely(
  omnibus_results,
  omnibus_out_path)

write_csv_safely(
  estimated_marginal_means,
  estimated_marginal_means_out_path)

write_csv_safely(
  pairwise_posthocs,
  pairwise_posthocs_out_path)

write_csv_safely(
  model_diagnostics,
  model_diagnostics_out_path)

write_csv_safely(
  sample_characteristics,
  sample_characteristics_out_path)

write_csv_safely(
  group_event_counts,
  group_event_counts_out_path)

message(
  "Updated common-comorbidity analyses completed successfully using the ",
  "exact published 3,688-participant sample. Results are in: ",
  results_dir)

message(
  "Primary multiplicity correction: BH-FDR across the 20 omnibus ",
  "comorbidity-group tests. Pairwise post-hoc p values are BH-FDR corrected ",
  "within each outcome across the 15 group contrasts, matching the published ",
  "post-hoc strategy.")
