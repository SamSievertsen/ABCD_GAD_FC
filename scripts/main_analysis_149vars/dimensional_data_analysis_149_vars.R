## Setup ##

# Load packages for mixed modeling, inference, multiplicity correction, and output
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)

options(
  digits = 8,
  scipen = 999,
  contrasts = c("contr.treatment", "contr.poly"))

emmeans::emm_options(
  lmer.df = "satterthwaite",
  lmerTest.limit = 5000,
  disable.lmerTest = FALSE,
  disable.pbkrtest = TRUE)


## Paths ##

dimensional_data_path <- "./data_processed/main_analysis/dimensional_analysis_data.csv"

metric_manifest_path <- "./results/results_149_vars/repeated_measures_dv_manifest_149vars.csv"

results_dir <- "./results/results_149_vars/dimensional"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cbcl_merged_out_path <-
  file.path(results_dir, "cbcl_merged_intx_analysis_results.csv")
cbcl_gad_out_path <-
  file.path(results_dir, "cbcl_gad_group_analysis_results.csv")
cbcl_hc_out_path <-
  file.path(results_dir, "cbcl_hc_group_analysis_results.csv")
bpm_merged_out_path <-
  file.path(results_dir, "bpm_merged_group_results.csv")
bpm_gad_out_path <-
  file.path(results_dir, "bpm_gad_group_results.csv")
bpm_hc_out_path <-
  file.path(results_dir, "bpm_control_group_results.csv")


## Constants ##

analysis_dataset <- "primary_current_complete_case"
dataset_label <- "Primary current complete-case dimensional cohort"

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

expected_cbcl_n <- 3146L
expected_cbcl_group_counts <- tibble::tribble(
  ~analysis_group, ~expected_n,
  "control", 2989L,
  "GAD", 157L)

expected_cbcl_group_event_counts <- tibble::tribble(
  ~analysis_group, ~eventname, ~expected_n,
  "control", "baseline_year_1_arm_1", 2003L,
  "control", "2_year_follow_up_y_arm_1", 986L,
  "GAD", "baseline_year_1_arm_1", 98L,
  "GAD", "2_year_follow_up_y_arm_1", 59L)


## Connectivity Outcome Manifest and Labels ##

metric_manifest <- readr::read_csv(
  metric_manifest_path,
  show_col_types = FALSE)

if (!all(c("analysis_order", "column_name") %in% names(metric_manifest))) {
  stop("The connectivity manifest lacks analysis_order or column_name.")
}

metric_manifest <- metric_manifest %>%
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
  stop("The outcome manifest must contain 20 unique outcomes.")
}

if (!identical(metric_manifest$analysis_order, seq_len(20L))) {
  stop("The outcome manifest analysis_order must be exactly 1 through 20.")
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
      " contains duplicated subjectkey-eventname rows.")
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
    stop(dataset_name, " contains incomplete required model data.")
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

get_numeric_column <- function(data, candidates) {
  matching_columns <- intersect(candidates, names(data))

  if (length(matching_columns) == 0) {
    return(rep(NA_real_, nrow(data)))
  }

  suppressWarnings(as.numeric(data[[matching_columns[[1]]]]))
}

normalize_term <- function(x) {
  x %>%
    stringr::str_replace_all("`", "") %>%
    stringr::str_replace_all("\\s", "")
}

safe_formula <- function(outcome, rhs) {
  stats::as.formula(paste(outcome, "~", rhs))
}

has_text <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

# Fit a mixed model while retaining warnings and errors for auditability.
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
    warnings = if (length(model_warnings) == 0) {
      NA_character_
    } else {
      paste(unique(model_warnings), collapse = " | ")
    },
    error = model_error)
}

extract_diagnostics <- function(fit_result, data) {
  model <- fit_result$model

  if (is.null(model)) {
    return(tibble::tibble(
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

extract_coefficient <- function(model, term) {
  empty_result <- tibble::tibble(
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    t_value = NA_real_,
    p_value = NA_real_,
    lower_ci = NA_real_,
    upper_ci = NA_real_)

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

  estimate <- get_numeric_column(row, "Estimate")[[1]]
  std_error <- get_numeric_column(row, c("Std. Error", "Std..Error"))[[1]]
  df <- get_numeric_column(row, "df")[[1]]
  t_value <- get_numeric_column(row, c("t value", "t.value"))[[1]]
  p_value <- get_numeric_column(row, c("Pr(>|t|)", "Pr...t.."))[[1]]

  critical_value <- if (is.finite(df)) {
    stats::qt(0.975, df = df)
  } else {
    stats::qnorm(0.975)
  }

  tibble::tibble(
    estimate = estimate,
    std_error = std_error,
    df = df,
    t_value = t_value,
    p_value = p_value,
    lower_ci = estimate - critical_value * std_error,
    upper_ci = estimate + critical_value * std_error)
}

extract_anova_term <- function(model, term_candidates, type) {
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

  normalized_candidates <- normalize_term(term_candidates)
  matched_row <- which(
    anova_table$normalized_term %in% normalized_candidates)

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

# Extract Control and GAD symptom slopes and their GAD-minus-Control difference
# from the published interaction model
extract_group_trends <- function(model, symptom) {
  empty_slope <- tibble::tibble(
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    lower_ci = NA_real_,
    upper_ci = NA_real_,
    t_value = NA_real_,
    p_value = NA_real_)

  empty_result <- list(
    control = empty_slope,
    gad = empty_slope,
    difference = empty_slope)

  if (is.null(model)) {
    return(empty_result)
  }

  tryCatch({
    trends <- emmeans::emtrends(
      model,
      ~ analysis_group,
      var = symptom)

    trend_summary <- as.data.frame(
      summary(
        trends,
        infer = c(TRUE, TRUE)))

    trend_column <- names(trend_summary)[
      stringr::str_detect(names(trend_summary), "\\.trend$")]

    if (length(trend_column) != 1L) {
      stop("Could not identify the emtrends slope column.")
    }

    trend_summary <- trend_summary %>%
      dplyr::transmute(
        analysis_group = as.character(analysis_group),
        estimate = suppressWarnings(as.numeric(.data[[trend_column]])),
        std_error = suppressWarnings(as.numeric(SE)),
        df = suppressWarnings(as.numeric(df)),
        lower_ci = suppressWarnings(as.numeric(lower.CL)),
        upper_ci = suppressWarnings(as.numeric(upper.CL)),
        t_value = suppressWarnings(as.numeric(t.ratio)),
        p_value = suppressWarnings(as.numeric(p.value)))

    control_row <- trend_summary %>%
      dplyr::filter(analysis_group == "control") %>%
      dplyr::select(-analysis_group)

    gad_row <- trend_summary %>%
      dplyr::filter(analysis_group == "GAD") %>%
      dplyr::select(-analysis_group)

    if (nrow(control_row) != 1L || nrow(gad_row) != 1L) {
      stop("Could not extract one Control and one GAD trend.")
    }

    difference_row <- emmeans::contrast(
      trends,
      method = list("GAD - control" = c(-1, 1)),
      adjust = "none") %>%
      summary(
        infer = c(TRUE, TRUE),
        adjust = "none") %>%
      as.data.frame() %>%
      dplyr::transmute(
        estimate = suppressWarnings(as.numeric(estimate)),
        std_error = suppressWarnings(as.numeric(SE)),
        df = suppressWarnings(as.numeric(df)),
        lower_ci = suppressWarnings(as.numeric(lower.CL)),
        upper_ci = suppressWarnings(as.numeric(upper.CL)),
        t_value = suppressWarnings(as.numeric(t.ratio)),
        p_value = suppressWarnings(as.numeric(p.value)))

    if (nrow(difference_row) != 1L) {
      stop("Could not extract one GAD-minus-Control slope contrast.")
    }

    list(
      control = control_row,
      gad = gad_row,
      difference = difference_row)
  }, error = function(error_condition) empty_result)
}

prefix_columns <- function(data, prefix) {
  names(data) <- paste0(prefix, names(data))
  data
}

add_fdr_columns <- function(data, p_columns) {
  for (p_column in p_columns) {
    if (!(p_column %in% names(data))) {
      stop("Cannot FDR-correct missing p-value column: ", p_column)
    }

    across_family_name <- paste0(p_column, "_fdr_across_family")
    within_symptom_name <- paste0(p_column, "_fdr_within_symptom")

    data[[across_family_name]] <- stats::p.adjust(
      data[[p_column]],
      method = "fdr")

    data <- data %>%
      dplyr::group_by(symptom_variable) %>%
      dplyr::mutate(
        !!within_symptom_name := stats::p.adjust(
          .data[[p_column]],
          method = "fdr")) %>%
      dplyr::ungroup()
  }

  data
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


## Read and Validate the Primary Dimensional Dataset ##

dimensional_data <- readr::read_csv(
  dimensional_data_path,
  show_col_types = FALSE)

required_variables <- c(
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
  cbcl_vars,
  bpm_vars)

assert_columns(
  dimensional_data,
  required_variables,
  "Primary dimensional dataset")

dimensional_data <- dimensional_data %>%
  dplyr::mutate(
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
      suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
    dplyr::across(
      dplyr::all_of(c(connectivity_vars, cbcl_vars, bpm_vars)),
      ~ suppressWarnings(as.numeric(.x))))

assert_unique_keys(dimensional_data, "Primary dimensional dataset")

unexpected_groups <- setdiff(
  unique(dimensional_data$analysis_group),
  analysis_groups)

if (length(unexpected_groups) > 0) {
  stop(
    "Unexpected analysis groups: ",
    paste(unexpected_groups, collapse = ", "),
    ".")
}

# The wrangling script already creates the simultaneous four-score CBCL
cbcl_analysis_data <- dimensional_data %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(c(
        connectivity_vars,
        cbcl_vars,
        "subjectkey",
        "family_id",
        "eventname",
        "analysis_group",
        "group",
        "interview_age",
        "age_in_years",
        "sex",
        "site_name",
        "rsfmri_c_ngd_meanmotion")),
      ~ !is.na(.x) & trimws(as.character(.x)) != ""))

if (nrow(cbcl_analysis_data) != expected_cbcl_n) {
  stop(
    "The CBCL analysis sample contains ",
    nrow(cbcl_analysis_data),
    " rows; expected ",
    expected_cbcl_n,
    ".")
}

if (dplyr::n_distinct(cbcl_analysis_data$subjectkey) != expected_cbcl_n) {
  stop("The CBCL analysis sample must contain one observation per participant.")
}

assert_exact_counts(
  data = cbcl_analysis_data,
  expected_counts = expected_cbcl_group_counts,
  grouping_variables = "analysis_group",
  dataset_name = "Primary CBCL analysis sample")

assert_exact_counts(
  data = cbcl_analysis_data,
  expected_counts = expected_cbcl_group_event_counts,
  grouping_variables = c("analysis_group", "eventname"),
  dataset_name = "Primary CBCL analysis sample")

assert_complete(
  cbcl_analysis_data,
  c(
    connectivity_vars,
    cbcl_vars,
    "subjectkey",
    "family_id",
    "eventname",
    "analysis_group",
    "sex",
    "site_name",
    "rsfmri_c_ngd_meanmotion"),
  "Primary CBCL analysis sample")

# Preserve the published BPM approach: the BPM sample is the complete BPM
# subset nested within the simultaneous four-score CBCL sample
bpm_analysis_data <- cbcl_analysis_data %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(bpm_vars),
      ~ !is.na(.x)))

if (nrow(bpm_analysis_data) == 0L) {
  stop("The primary BPM analysis sample is empty.")
}

if (dplyr::n_distinct(bpm_analysis_data$analysis_group) != 2L) {
  stop("The primary BPM analysis sample does not contain both groups.")
}

# Match the published factor construction. Control is explicitly retained as
# the reference diagnostic group
prepare_factors <- function(data) {
  data %>%
    dplyr::mutate(
      subjectkey = factor(subjectkey),
      family_id = factor(family_id),
      eventname = factor(eventname),
      analysis_group = factor(
        analysis_group,
        levels = analysis_groups),
      group = factor(group),
      sex = factor(sex),
      site_name = factor(site_name))
}

cbcl_analysis_data <- prepare_factors(cbcl_analysis_data)
bpm_analysis_data <- prepare_factors(bpm_analysis_data)

cbcl_gad_n <- sum(cbcl_analysis_data$analysis_group == "GAD")
cbcl_control_n <- sum(cbcl_analysis_data$analysis_group == "control")
bpm_gad_n <- sum(bpm_analysis_data$analysis_group == "GAD")
bpm_control_n <- sum(bpm_analysis_data$analysis_group == "control")

message(
  "CBCL analysis sample validated: ",
  nrow(cbcl_analysis_data),
  " participants (",
  cbcl_control_n,
  " Controls; ",
  cbcl_gad_n,
  " GAD).")

message(
  "BPM analysis sample validated: ",
  nrow(bpm_analysis_data),
  " participants (",
  bpm_control_n,
  " Controls; ",
  bpm_gad_n,
  " GAD).")

print(
  bpm_analysis_data %>%
    dplyr::count(analysis_group, eventname, name = "n"))


## Model Functions ##

# Fit one symptom-outcome pair using:
# 1. the published merged interaction specification;
# 2. a no-interaction pooled model with otherwise identical specification;
# 3. the published GAD-only specification; and
# 4. the published Control-only specification
fit_dimensional_pair <- function(
    data,
    symptom,
    outcome,
    analysis_order,
    outcome_label,
    family_name) {

  is_cbcl <- identical(family_name, "CBCL")

  if (is_cbcl) {
    interaction_rhs <- paste0(
      symptom,
      " * analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + ",
      "(1 | site_name) + (1 | family_id)")

    pooled_rhs <- paste0(
      symptom,
      " + analysis_group + rsfmri_c_ngd_meanmotion + sex + eventname + ",
      "(1 | site_name) + (1 | family_id)")

    subgroup_rhs <- paste0(
      symptom,
      " + rsfmri_c_ngd_meanmotion + sex + eventname + ",
      "(1 | site_name) + (1 | family_id)")
  } else {
    
    # Published BPM specifications: no event fixed effect in the merged
    # model and a site-only random intercept in the subgroup models
    interaction_rhs <- paste0(
      symptom,
      " * analysis_group + rsfmri_c_ngd_meanmotion + sex + ",
      "(1 | site_name) + (1 | family_id)")

    pooled_rhs <- paste0(
      symptom,
      " + analysis_group + rsfmri_c_ngd_meanmotion + sex + ",
      "(1 | site_name) + (1 | family_id)")

    subgroup_rhs <- paste0(
      symptom,
      " + rsfmri_c_ngd_meanmotion + sex + (1 | site_name)")
  }

  interaction_formula <- safe_formula(outcome, interaction_rhs)
  pooled_formula <- safe_formula(outcome, pooled_rhs)
  subgroup_formula <- safe_formula(outcome, subgroup_rhs)

  interaction_fit <- fit_lmer_safe(interaction_formula, data)
  pooled_fit <- fit_lmer_safe(pooled_formula, data)

  gad_data <- data %>%
    dplyr::filter(analysis_group == "GAD") %>%
    droplevels()

  control_data <- data %>%
    dplyr::filter(analysis_group == "control") %>%
    droplevels()

  gad_fit <- fit_lmer_safe(subgroup_formula, gad_data)
  control_fit <- fit_lmer_safe(subgroup_formula, control_data)

  interaction_term_candidates <- c(
    paste0(symptom, ":analysis_group"),
    paste0("analysis_group:", symptom))

  interaction_test <- extract_anova_term(
    interaction_fit$model,
    interaction_term_candidates,
    type = "III")

  diagnosis_group_test <- extract_anova_term(
    interaction_fit$model,
    "analysis_group",
    type = "III")

  interaction_motion <- extract_coefficient(
    interaction_fit$model,
    "rsfmri_c_ngd_meanmotion")

  interaction_sex <- extract_anova_term(
    interaction_fit$model,
    "sex",
    type = "III")

  interaction_event <- if (is_cbcl) {
    extract_anova_term(
      interaction_fit$model,
      "eventname",
      type = "III")
  } else {
    tibble::tibble(
      f_value = NA_real_,
      numerator_df = NA_real_,
      denominator_df = NA_real_,
      p_value = NA_real_)
  }

  group_trends <- extract_group_trends(
    interaction_fit$model,
    symptom)

  pooled_slope <- extract_coefficient(
    pooled_fit$model,
    symptom)

  gad_slope <- extract_coefficient(gad_fit$model, symptom)
  gad_motion <- extract_coefficient(
    gad_fit$model,
    "rsfmri_c_ngd_meanmotion")
  gad_sex <- extract_anova_term(gad_fit$model, "sex", type = "II")
  gad_event <- if (is_cbcl) {
    extract_anova_term(gad_fit$model, "eventname", type = "II")
  } else {
    tibble::tibble(
      f_value = NA_real_,
      numerator_df = NA_real_,
      denominator_df = NA_real_,
      p_value = NA_real_)
  }

  control_slope <- extract_coefficient(control_fit$model, symptom)
  control_motion <- extract_coefficient(
    control_fit$model,
    "rsfmri_c_ngd_meanmotion")
  control_sex <- extract_anova_term(
    control_fit$model,
    "sex",
    type = "II")
  control_event <- if (is_cbcl) {
    extract_anova_term(control_fit$model, "eventname", type = "II")
  } else {
    tibble::tibble(
      f_value = NA_real_,
      numerator_df = NA_real_,
      denominator_df = NA_real_,
      p_value = NA_real_)
  }

  base_columns <- tibble::tibble(
    analysis_dataset = analysis_dataset,
    dataset_label = dataset_label,
    symptom_family = family_name,
    symptom_variable = symptom,
    analysis_order = analysis_order,
    column_name = outcome,
    DV_Summary = outcome_label,
    published_interaction_formula =
      paste(deparse(interaction_formula), collapse = " "),
    added_pooled_formula = paste(deparse(pooled_formula), collapse = " "),
    published_subgroup_formula =
      paste(deparse(subgroup_formula), collapse = " "))

  merged_result <- dplyr::bind_cols(
    base_columns,
    prefix_columns(interaction_test, "interaction_"),
    prefix_columns(
      group_trends$control,
      "control_slope_from_interaction_"),
    prefix_columns(
      group_trends$gad,
      "gad_slope_from_interaction_"),
    prefix_columns(
      group_trends$difference,
      "gad_minus_control_slope_"),
    prefix_columns(pooled_slope, "pooled_slope_"),
    prefix_columns(diagnosis_group_test, "diagnosis_group_"),
    prefix_columns(interaction_motion, "interaction_motion_"),
    prefix_columns(interaction_sex, "interaction_sex_"),
    prefix_columns(interaction_event, "interaction_event_"),
    prefix_columns(
      extract_diagnostics(interaction_fit, data),
      "interaction_model_"),
    prefix_columns(
      extract_diagnostics(pooled_fit, data),
      "pooled_model_"))

  gad_result <- dplyr::bind_cols(
    base_columns,
    tibble::tibble(analysis_group = "GAD"),
    prefix_columns(gad_slope, "symptom_slope_"),
    prefix_columns(gad_motion, "motion_"),
    prefix_columns(gad_sex, "sex_"),
    prefix_columns(gad_event, "event_"),
    prefix_columns(
      extract_diagnostics(gad_fit, gad_data),
      "model_"))

  control_result <- dplyr::bind_cols(
    base_columns,
    tibble::tibble(analysis_group = "control"),
    prefix_columns(control_slope, "symptom_slope_"),
    prefix_columns(control_motion, "motion_"),
    prefix_columns(control_sex, "sex_"),
    prefix_columns(control_event, "event_"),
    prefix_columns(
      extract_diagnostics(control_fit, control_data),
      "model_"))

  list(
    merged = merged_result,
    gad = gad_result,
    control = control_result)
}

run_dimensional_family <- function(
    data,
    symptom_variables,
    family_name) {

  model_grid <- tidyr::crossing(
    symptom_variable = symptom_variables,
    metric_manifest %>%
      dplyr::select(
        analysis_order,
        column_name,
        DV_Summary)) %>%
    dplyr::arrange(
      match(symptom_variable, symptom_variables),
      analysis_order)

  model_results <- purrr::pmap(
    model_grid,
    function(
        symptom_variable,
        analysis_order,
        column_name,
        DV_Summary) {
      message(
        family_name,
        ": modeling ",
        DV_Summary,
        " ~ ",
        symptom_variable)

      fit_dimensional_pair(
        data = data,
        symptom = symptom_variable,
        outcome = column_name,
        analysis_order = analysis_order,
        outcome_label = DV_Summary,
        family_name = family_name)
    })

  merged_results <- purrr::map_dfr(model_results, "merged") %>%
    dplyr::arrange(
      match(symptom_variable, symptom_variables),
      analysis_order) %>%
    add_fdr_columns(
      c(
        "interaction_p_value",
        "control_slope_from_interaction_p_value",
        "gad_slope_from_interaction_p_value",
        "gad_minus_control_slope_p_value",
        "pooled_slope_p_value"))

  gad_results <- purrr::map_dfr(model_results, "gad") %>%
    dplyr::arrange(
      match(symptom_variable, symptom_variables),
      analysis_order) %>%
    add_fdr_columns("symptom_slope_p_value")

  control_results <- purrr::map_dfr(model_results, "control") %>%
    dplyr::arrange(
      match(symptom_variable, symptom_variables),
      analysis_order) %>%
    add_fdr_columns("symptom_slope_p_value")

  list(
    merged = merged_results,
    gad = gad_results,
    control = control_results)
}

validate_result_family <- function(
    results,
    family_name,
    expected_rows,
    expected_merged_n,
    expected_gad_n,
    expected_control_n) {

  if (
    nrow(results$merged) != expected_rows ||
      nrow(results$gad) != expected_rows ||
      nrow(results$control) != expected_rows) {
    stop(
      family_name,
      " did not return the expected ",
      expected_rows,
      " rows per result family.")
  }

  merged_error_rows <- results$merged %>%
    dplyr::filter(
      has_text(interaction_model_error) |
        has_text(pooled_model_error))

  gad_error_rows <- results$gad %>%
    dplyr::filter(has_text(model_error))

  control_error_rows <- results$control %>%
    dplyr::filter(has_text(model_error))

  if (
    nrow(merged_error_rows) > 0 ||
      nrow(gad_error_rows) > 0 ||
      nrow(control_error_rows) > 0) {
    print(
      tibble::tibble(
        result_family = c("merged", "GAD", "control"),
        n_model_error_rows = c(
          nrow(merged_error_rows),
          nrow(gad_error_rows),
          nrow(control_error_rows))))
    stop(family_name, " contains one or more model-fitting errors.")
  }

  merged_convergence_rows <- results$merged %>%
    dplyr::filter(
      has_text(interaction_model_convergence_messages) |
        has_text(pooled_model_convergence_messages))

  gad_convergence_rows <- results$gad %>%
    dplyr::filter(has_text(model_convergence_messages))

  control_convergence_rows <- results$control %>%
    dplyr::filter(has_text(model_convergence_messages))

  if (
    nrow(merged_convergence_rows) > 0 ||
      nrow(gad_convergence_rows) > 0 ||
      nrow(control_convergence_rows) > 0) {
    print(
      tibble::tibble(
        result_family = c("merged", "GAD", "control"),
        n_convergence_rows = c(
          nrow(merged_convergence_rows),
          nrow(gad_convergence_rows),
          nrow(control_convergence_rows))))
    stop(
      family_name,
      " contains one or more optimizer convergence messages. No outputs were written.")
  }

  optimizer_failed <- function(x) {
    has_text(x) & trimws(as.character(x)) != "0"
  }

  merged_optimizer_rows <- results$merged %>%
    dplyr::filter(
      optimizer_failed(interaction_model_optimizer_code) |
        optimizer_failed(pooled_model_optimizer_code))

  gad_optimizer_rows <- results$gad %>%
    dplyr::filter(optimizer_failed(model_optimizer_code))

  control_optimizer_rows <- results$control %>%
    dplyr::filter(optimizer_failed(model_optimizer_code))

  if (
    nrow(merged_optimizer_rows) > 0 ||
      nrow(gad_optimizer_rows) > 0 ||
      nrow(control_optimizer_rows) > 0) {
    print(
      tibble::tibble(
        result_family = c("merged", "GAD", "control"),
        n_nonzero_optimizer_code_rows = c(
          nrow(merged_optimizer_rows),
          nrow(gad_optimizer_rows),
          nrow(control_optimizer_rows))))
    stop(
      family_name,
      " contains one or more nonzero optimizer return codes. No outputs were written.")
  }

  required_merged_inference <- c(
    "interaction_p_value",
    "control_slope_from_interaction_p_value",
    "gad_slope_from_interaction_p_value",
    "gad_minus_control_slope_p_value",
    "pooled_slope_p_value")

  merged_missing_inference <- vapply(
    required_merged_inference,
    function(variable) sum(!is.finite(results$merged[[variable]])),
    integer(1))

  gad_missing_inference <- sum(
    !is.finite(results$gad$symptom_slope_p_value))

  control_missing_inference <- sum(
    !is.finite(results$control$symptom_slope_p_value))

  if (
    any(merged_missing_inference > 0) ||
      gad_missing_inference > 0 ||
      control_missing_inference > 0) {
    print(merged_missing_inference[merged_missing_inference > 0])
    print(
      c(
        gad_missing_inference = gad_missing_inference,
        control_missing_inference = control_missing_inference))
    stop(
      family_name,
      " contains missing required inferential results. No outputs were written.")
  }

  if (
    any(results$merged$interaction_model_n_observations != expected_merged_n) ||
      any(results$merged$pooled_model_n_observations != expected_merged_n) ||
      any(results$gad$model_n_observations != expected_gad_n) ||
      any(results$control$model_n_observations != expected_control_n)) {
    stop(
      family_name,
      " model observation counts do not match the validated analysis samples.")
  }

  singular_summary <- tibble::tibble(
    model_family = c(
      "merged interaction",
      "merged pooled",
      "GAD subgroup",
      "Control subgroup"),
    n_singular = c(
      sum(results$merged$interaction_model_singular %in% TRUE),
      sum(results$merged$pooled_model_singular %in% TRUE),
      sum(results$gad$model_singular %in% TRUE),
      sum(results$control$model_singular %in% TRUE)))

  message(family_name, " model validation passed.")
  print(singular_summary)

  invisible(singular_summary)
}


## Run CBCL Models ##

cbcl_results <- run_dimensional_family(
  data = cbcl_analysis_data,
  symptom_variables = cbcl_vars,
  family_name = "CBCL")

# Four symptoms x 20 outcomes = 80 tests per inferential family.
validate_result_family(
  results = cbcl_results,
  family_name = "CBCL",
  expected_rows = 80L,
  expected_merged_n = nrow(cbcl_analysis_data),
  expected_gad_n = cbcl_gad_n,
  expected_control_n = cbcl_control_n)


## Run BPM Models ##

bpm_results <- run_dimensional_family(
  data = bpm_analysis_data,
  symptom_variables = bpm_vars,
  family_name = "BPM")

# Two symptoms x 20 outcomes = 40 tests per inferential family.
validate_result_family(
  results = bpm_results,
  family_name = "BPM",
  expected_rows = 40L,
  expected_merged_n = nrow(bpm_analysis_data),
  expected_gad_n = bpm_gad_n,
  expected_control_n = bpm_control_n)


## Output ##

# Write none of the six result files until both complete model families have
# finished and passed all row-count, model-error, convergence, and inference
# checks. FDR corrections are performed independently within each primary
# inferential family
write_csv_safely(cbcl_results$merged, cbcl_merged_out_path)
write_csv_safely(cbcl_results$gad, cbcl_gad_out_path)
write_csv_safely(cbcl_results$control, cbcl_hc_out_path)
write_csv_safely(bpm_results$merged, bpm_merged_out_path)
write_csv_safely(bpm_results$gad, bpm_gad_out_path)
write_csv_safely(bpm_results$control, bpm_hc_out_path)


## Completion Summary ##

message(
  "Primary dimensional analyses completed successfully. Results are in: ",
  results_dir)

message(
  "Published model specifications were replicated across the 20 corrected ",
  "outcomes. Added pooled models test the actual group-adjusted whole-sample ",
  "symptom association without a symptom-by-group interaction.")

message(
  "Principal BH-FDR families: 80 tests for each CBCL inferential family and ",
  "40 tests for each BPM inferential family. Secondary within-symptom FDR ",
  "columns correct across the 20 outcomes for each symptom.")
