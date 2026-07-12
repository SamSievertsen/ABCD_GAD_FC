## Setup ##

# Load packages for reading, wrangling, modeling, post hoc testing, and plotting
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(writexl)

# Set numeric display options
options(digits = 6, scipen = 999)

emmeans::emm_options(
  lmer.df = "satterthwaite",
  lmerTest.limit = 5000,
  disable.lmerTest = FALSE,
  disable.pbkrtest = TRUE)

## Paths ##

#1.1 Primary QC-complete repeated measures dataset
primary_data_path <- "./data_processed/main_analysis/repeated_measures_grouped_imaging_data_149vars.csv"

#1.2 Exact published-cohort sensitivity dataset with corrected connectivity values
sensitivity_data_path <- paste0(
  "./data_processed/main_analysis/",
  "repeated_measures_grouped_imaging_data_149vars_",
  "published_cohort_sensitivity.csv")

#1.3 Ordered manifest of the 20 connectivity outcomes
metric_manifest_path <- "./results/results_149_vars/repeated_measures_dv_manifest_149vars.csv"

#1.4 Race and ethnicity data
ethnicity_data_path <- "./data_raw/acspsw03.txt"

#1.5 CBCL data
cbcl_data_path <- "./data_raw/abcd_cbcls01.csv"

#1.6 Root output directory
results_root <- "./results/results_149_vars/repeated_measures"

#1.7 Dataset-specific output directories
primary_results_dir <- file.path(results_root, "primary")
sensitivity_results_dir <- file.path(
  results_root,
  "published_cohort_sensitivity")

comparison_results_dir <- file.path(results_root, "comparison")

#1.8 Create output directories
purrr::walk(
  c(
    results_root,
    primary_results_dir,
    sensitivity_results_dir,
    comparison_results_dir,
    file.path(primary_results_dir, "plots"),
    file.path(sensitivity_results_dir, "plots")),
  ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))


## Constants ##

#2.1 Define the expected assessment waves and their order
analysis_events <- c(
  "baseline_year_1_arm_1",
  "2_year_follow_up_y_arm_1")

#2.2 Define the diagnostic trajectory group order
analysis_groups <- c(
  "Control",
  "Continuous GAD",
  "GAD Converter",
  "GAD Remitter")

#2.3 Define expected subject counts for the primary cohort
expected_primary_group_counts <- tibble::tibble(
  group = analysis_groups,
  expected_n = c(1661L, 12L, 46L, 48L))

#2.4 Define expected subject counts for the published-cohort sensitivity sample
expected_sensitivity_group_counts <- tibble::tibble(
  group = analysis_groups,
  expected_n = c(1955L, 12L, 55L, 52L))

#2.5 Define consistent plot colors
plot_colors <- c(
  "HC" = "#218380",
  "Continuous GAD" = "#D81159",
  "GAD Converter" = "#FDB833",
  "GAD Remitter" = "#5390D9")


## Connectivity Outcome Manifest ##

#3. Read the ordered manifest created by the repeated measures wrangling script
metric_manifest <- readr::read_csv(
  metric_manifest_path,
  show_col_types = FALSE) %>%
  dplyr::mutate(
    analysis_order = as.integer(analysis_order),
    column_name = trimws(as.character(column_name))) %>%
  dplyr::arrange(analysis_order)

#3.1 Confirm that exactly 20 unique outcomes are present
if (
  nrow(metric_manifest) != 20 ||
    dplyr::n_distinct(metric_manifest$column_name) != 20) {
  stop("The repeated measures outcome manifest must contain exactly 20 unique metrics.")
}

#3.2 Add human-readable labels for all 20 connectivity outcomes
metric_labels <- tibble::tribble(
  ~column_name, ~metric_label,
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

if (any(is.na(metric_manifest$metric_label))) {
  stop(
    "Human-readable labels are missing for: ",
    paste(
      metric_manifest$column_name[is.na(metric_manifest$metric_label)],
      collapse = ", "))
}

sig_dvs <- metric_manifest$column_name


## Helper Functions ##

#4. Return the first matching column from a dataframe or an NA vector
extract_numeric_column <- function(data, candidate_names) {
  matched_name <- intersect(candidate_names, names(data))

  if (length(matched_name) == 0) {
    return(rep(NA_real_, nrow(data)))
  }

  suppressWarnings(as.numeric(data[[matched_name[[1]]]]))
}

#4.1 Create a safe filename from a metric label
sanitize_filename <- function(value) {
  value %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

#4.2 Validate and format a repeated measures analysis dataset
prepare_analysis_data <- function(
    data_path,
    dataset_name,
    expected_n_subjects,
    expected_n_rows,
    expected_group_counts
) {

  analysis_data <- read.csv(
    data_path,
    stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      subjectkey = trimws(as.character(subjectkey)),
      eventname = trimws(as.character(eventname)),
      group = trimws(as.character(group)),
      sex = trimws(as.character(sex)),
      site_name = trimws(as.character(site_name)),
      family_id = trimws(as.character(family_id)))

  required_vars <- c(
    "subjectkey",
    "eventname",
    "group",
    "interview_age",
    "age_in_years",
    "sex",
    "site_name",
    "family_id",
    "rsfmri_c_ngd_meanmotion",
    sig_dvs)

  missing_vars <- setdiff(required_vars, names(analysis_data))

  if (length(missing_vars) > 0) {
    stop(
      dataset_name,
      " is missing required variables: ",
      paste(missing_vars, collapse = ", "))
  }

  duplicate_keys <- analysis_data %>%
    dplyr::count(subjectkey, eventname) %>%
    dplyr::filter(n > 1)

  if (nrow(duplicate_keys) > 0) {
    stop(dataset_name, " contains duplicated subject-event rows.")
  }

  analysis_data <- analysis_data %>%
    dplyr::mutate(
      interview_age = suppressWarnings(as.numeric(interview_age)),
      age_in_years = suppressWarnings(as.numeric(age_in_years)),
      rsfmri_c_ngd_meanmotion =
        suppressWarnings(as.numeric(rsfmri_c_ngd_meanmotion)),
      dplyr::across(
        dplyr::all_of(sig_dvs),
        ~ suppressWarnings(as.numeric(.x))))

  complete_case_vars <- c(
    "subjectkey",
    "eventname",
    "group",
    "interview_age",
    "age_in_years",
    "sex",
    "site_name",
    "family_id",
    "rsfmri_c_ngd_meanmotion",
    sig_dvs)

  incomplete_rows <- analysis_data %>%
    dplyr::filter(
      !dplyr::if_all(
        dplyr::all_of(complete_case_vars),
        ~ !is.na(.x) & trimws(as.character(.x)) != ""))

  if (nrow(incomplete_rows) > 0) {
    stop(
      dataset_name,
      " contains ",
      nrow(incomplete_rows),
      " incomplete rows across required model variables.")
  }

  subject_wave_check <- analysis_data %>%
    dplyr::group_by(subjectkey) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_events = dplyr::n_distinct(eventname),
      has_baseline = any(eventname == analysis_events[[1]]),
      has_followup = any(eventname == analysis_events[[2]]),
      .groups = "drop")

  invalid_waves <- subject_wave_check %>%
    dplyr::filter(
      n_rows != 2 |
        n_events != 2 |
        !has_baseline |
        !has_followup)

  if (nrow(invalid_waves) > 0) {
    stop(
      dataset_name,
      " contains participants without exactly one baseline and one follow-up row.")
  }

  observed_n_subjects <- dplyr::n_distinct(analysis_data$subjectkey)
  observed_n_rows <- nrow(analysis_data)

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

  observed_group_counts <- analysis_data %>%
    dplyr::distinct(subjectkey, group) %>%
    dplyr::count(group, name = "observed_n")

  group_count_check <- expected_group_counts %>%
    dplyr::full_join(observed_group_counts, by = "group") %>%
    dplyr::mutate(
      expected_n = tidyr::replace_na(expected_n, 0L),
      observed_n = tidyr::replace_na(observed_n, 0L),
      matches_expected = expected_n == observed_n)

  if (any(!group_count_check$matches_expected)) {
    print(group_count_check)
    stop(dataset_name, " does not have the expected group counts.")
  }

  unexpected_groups <- setdiff(unique(analysis_data$group), analysis_groups)
  unexpected_events <- setdiff(unique(analysis_data$eventname), analysis_events)

  if (length(unexpected_groups) > 0) {
    stop(
      dataset_name,
      " contains unexpected groups: ",
      paste(unexpected_groups, collapse = ", "))
  }

  if (length(unexpected_events) > 0) {
    stop(
      dataset_name,
      " contains unexpected events: ",
      paste(unexpected_events, collapse = ", "))
  }

  analysis_data <- analysis_data %>%
    dplyr::mutate(
      eventname = factor(eventname, levels = analysis_events),
      group = factor(group, levels = analysis_groups),
      sex = factor(sex),
      site_name = factor(site_name),
      family_id = factor(family_id),
      subjectkey = factor(subjectkey))

  message(
    dataset_name,
    " validated: ",
    observed_n_subjects,
    " participants and ",
    observed_n_rows,
    " rows.")

  analysis_data
}

#4.3 Fit one repeated measures mixed model and extract model results
fit_repeated_measures_model <- function(
    analysis_data,
    dataset_key,
    dataset_label,
    dv,
    metric_label,
    analysis_order) {

  model_warnings <- character()
  model_error <- NA_character_

  model_formula <- stats::as.formula(
    paste0(
      dv,
      " ~ group * eventname + ",
      "rsfmri_c_ngd_meanmotion + sex + ",
      "(1 | site_name) + (1 | family_id) + (1 | subjectkey)"))

  fitted_model <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        formula = model_formula,
        data = analysis_data,
        na.action = na.omit,
        REML = TRUE),
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

  empty_interaction <- tibble::tibble(
    dataset = dataset_key,
    dataset_label = dataset_label,
    analysis_order = analysis_order,
    DV = dv,
    DV_Summary = metric_label,
    term = "group:eventname",
    f_value = NA_real_,
    numerator_df = NA_real_,
    denominator_df = NA_real_,
    p_value = NA_real_)

  if (is.null(fitted_model)) {
    diagnostic_row <- tibble::tibble(
      dataset = dataset_key,
      dataset_label = dataset_label,
      analysis_order = analysis_order,
      DV = dv,
      DV_Summary = metric_label,
      n_observations = NA_integer_,
      n_subjects = dplyr::n_distinct(analysis_data$subjectkey),
      singular = NA,
      convergence_messages = NA_character_,
      warnings = paste(unique(model_warnings), collapse = " | "),
      error = model_error,
      AIC = NA_real_,
      BIC = NA_real_,
      logLik = NA_real_,
      marginal_R2 = NA_real_,
      conditional_R2 = NA_real_)

    return(list(
      model = NULL,
      interaction = empty_interaction,
      anova = tibble::tibble(),
      fixed_effects = tibble::tibble(),
      diagnostics = diagnostic_row))
  }

  anova_object <- tryCatch(
    car::Anova(
      fitted_model,
      type = "III",
      test.statistic = "F"),
    error = function(error_condition) {
      model_error <<- paste0(
        "Type III ANCOVA failed: ",
        conditionMessage(error_condition))
      NULL
    })

  if (is.null(anova_object)) {
    anova_results <- tibble::tibble()
    interaction_result <- empty_interaction
  } else {
    anova_raw <- tibble::as_tibble(
      as.data.frame(anova_object),
      rownames = "term")

    anova_results <- tibble::tibble(
      dataset = dataset_key,
      dataset_label = dataset_label,
      analysis_order = analysis_order,
      DV = dv,
      DV_Summary = metric_label,
      term = anova_raw$term,
      sum_sq = extract_numeric_column(
        anova_raw,
        c("Sum Sq", "Sum.Sq")),
      numerator_df = extract_numeric_column(
        anova_raw,
        c("Df", "num Df", "NumDF")),
      denominator_df = extract_numeric_column(
        anova_raw,
        c("Df.res", "den Df", "DenDF")),
      f_value = extract_numeric_column(
        anova_raw,
        c("F", "F value", "F.value")),
      p_value = extract_numeric_column(
        anova_raw,
        c("Pr(>F)", "Pr..F.")))

    interaction_result <- anova_results %>%
      dplyr::filter(
        stringr::str_replace_all(term, "\\s", "") %in%
          c("group:eventname", "eventname:group")) %>%
      dplyr::select(
        dataset,
        dataset_label,
        analysis_order,
        DV,
        DV_Summary,
        term,
        f_value,
        numerator_df,
        denominator_df,
        p_value)

    if (nrow(interaction_result) != 1) {
      interaction_result <- empty_interaction
    }
  }

  coefficient_matrix <- summary(fitted_model)$coefficients

  fixed_effects <- tibble::as_tibble(
    as.data.frame(coefficient_matrix),
    rownames = "term")

  fixed_effects <- tibble::tibble(
    dataset = dataset_key,
    dataset_label = dataset_label,
    analysis_order = analysis_order,
    DV = dv,
    DV_Summary = metric_label,
    term = fixed_effects$term,
    estimate = extract_numeric_column(
      fixed_effects,
      c("Estimate")),
    std_error = extract_numeric_column(
      fixed_effects,
      c("Std. Error", "Std..Error")),
    df = extract_numeric_column(
      fixed_effects,
      c("df")),
    t_value = extract_numeric_column(
      fixed_effects,
      c("t value", "t.value")),
    p_value = extract_numeric_column(
      fixed_effects,
      c("Pr(>|t|)", "Pr...t..")))

  convergence_messages <- fitted_model@optinfo$conv$lme4$messages

  if (is.null(convergence_messages)) {
    convergence_messages <- NA_character_
  } else {
    convergence_messages <- paste(
      convergence_messages,
      collapse = " | ")
  }

  marginal_r2 <- NA_real_
  conditional_r2 <- NA_real_

  if (requireNamespace("performance", quietly = TRUE)) {
    r2_result <- tryCatch(
      performance::r2_nakagawa(fitted_model),
      error = function(error_condition) NULL)

    if (!is.null(r2_result)) {
      r2_result <- as.data.frame(r2_result)

      if ("R2_marginal" %in% names(r2_result)) {
        marginal_r2 <- as.numeric(r2_result$R2_marginal[[1]])
      }

      if ("R2_conditional" %in% names(r2_result)) {
        conditional_r2 <- as.numeric(r2_result$R2_conditional[[1]])
      }
    }
  }

  diagnostic_row <- tibble::tibble(
    dataset = dataset_key,
    dataset_label = dataset_label,
    analysis_order = analysis_order,
    DV = dv,
    DV_Summary = metric_label,
    n_observations = stats::nobs(fitted_model),
    n_subjects = dplyr::n_distinct(analysis_data$subjectkey),
    singular = lme4::isSingular(fitted_model, tol = 1e-4),
    convergence_messages = convergence_messages,
    warnings = ifelse(
      length(model_warnings) == 0,
      NA_character_,
      paste(unique(model_warnings), collapse = " | ")),
    error = model_error,
    AIC = stats::AIC(fitted_model),
    BIC = stats::BIC(fitted_model),
    logLik = as.numeric(stats::logLik(fitted_model)),
    marginal_R2 = marginal_r2,
    conditional_R2 = conditional_r2)

  list(
    model = fitted_model,
    interaction = interaction_result,
    anova = anova_results,
    fixed_effects = fixed_effects,
    diagnostics = diagnostic_row)
}

#4.4 Generate EMMs, post hoc tests, and a plot for one nominally significant outcome
run_nominal_posthocs <- function(
    fitted_model,
    analysis_data,
    dataset_key,
    dataset_label,
    dv,
    metric_label,
    analysis_order,
    plot_directory) {

  posthoc_error <- NA_character_

  posthoc_output <- tryCatch({

    emm_grid <- emmeans::emmeans(
      fitted_model,
      ~ eventname | group)
    
    emm_results <- as.data.frame(
      summary(
        emm_grid,
        infer = c(TRUE, TRUE))) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        dataset = dataset_key,
        dataset_label = dataset_label,
        analysis_order = analysis_order,
        DV = dv,
        DV_Summary = metric_label,
        .before = 1)

    #4.4.1 Estimate follow-up minus baseline within each diagnostic trajectory group
    within_group_change <- emmeans::contrast(
      emm_grid,
      method = "revpairwise",
      by = "group",
      adjust = "none")

    within_group_results <- as.data.frame(
      summary(
        within_group_change,
        infer = c(TRUE, TRUE),
        adjust = "none")) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        contrast_family = "Within-group follow-up minus baseline",
        dataset = dataset_key,
        dataset_label = dataset_label,
        analysis_order = analysis_order,
        DV = dv,
        DV_Summary = metric_label,
        .before = 1)

    #4.4.2 Compare longitudinal change across every pair of diagnostic groups
    between_group_change <- emmeans::contrast(
      emm_grid,
      interaction = c("revpairwise", "pairwise"),
      by = NULL,
      adjust = "none")

    between_group_results <- as.data.frame(
      summary(
        between_group_change,
        infer = c(TRUE, TRUE),
        adjust = "none")) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        contrast_family = "Between-group difference in longitudinal change",
        dataset = dataset_key,
        dataset_label = dataset_label,
        analysis_order = analysis_order,
        DV = dv,
        DV_Summary = metric_label,
        .before = 1)

    #4.4.3 Store raw p-values and apply FDR across all post hoc tests for this outcome
    posthoc_results <- dplyr::bind_rows(
      within_group_results,
      between_group_results) %>%
      dplyr::mutate(
        p_fdr_within_outcome = stats::p.adjust(
          p.value,
          method = "fdr")) %>%
      dplyr::group_by(contrast_family) %>%
      dplyr::mutate(
        p_fdr_within_family = stats::p.adjust(
          p.value,
          method = "fdr")) %>%
      dplyr::ungroup()

    #4.4.4 Prepare EMM plot data
    plot_data <- emm_results %>%
      dplyr::mutate(
        group_display = dplyr::recode(
          as.character(group),
          "Control" = "HC"),
        group_display = factor(
          group_display,
          levels = c(
            "HC",
            "Continuous GAD",
            "GAD Converter",
            "GAD Remitter")),
        time_numeric = dplyr::case_when(
          as.character(eventname) == analysis_events[[1]] ~ 1,
          as.character(eventname) == analysis_events[[2]] ~ 2,
          TRUE ~ NA_real_))

    subject_counts <- analysis_data %>%
      dplyr::distinct(subjectkey, group) %>%
      dplyr::count(group, name = "n") %>%
      dplyr::mutate(
        group_display = dplyr::recode(
          as.character(group),
          "Control" = "HC"),
        legend_label = paste0(group_display, " (n = ", n, ")"))

    legend_labels <- stats::setNames(
      subject_counts$legend_label,
      subject_counts$group_display)

    outcome_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = time_numeric,
        y = emmean,
        group = group_display,
        color = group_display,
        fill = group_display)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = lower.CL,
          ymax = upper.CL),
        alpha = 0.07,
        linewidth = 0) +
      ggplot2::geom_line(
        linewidth = 0.7,
        alpha = 0.7) +
      ggplot2::geom_point(
        size = 2.2,
        shape = 21) +
      ggplot2::scale_color_manual(
        name = "Diagnostic Group",
        values = plot_colors,
        labels = legend_labels) +
      ggplot2::scale_fill_manual(
        name = "Diagnostic Group",
        values = plot_colors,
        labels = legend_labels) +
      ggplot2::scale_x_continuous(
        breaks = c(1, 2),
        labels = c("Baseline", "2 Year Follow-up"),
        expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
      ggplot2::labs(
        x = "Timepoint",
        y = paste0(metric_label, " Connectivity"),
        title = metric_label,
        subtitle = dataset_label) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold"),
        plot.margin = grid::unit(c(1, 1, 1, 1), "lines"))

    safe_metric_name <- sanitize_filename(metric_label)

    png_path <- file.path(
      plot_directory,
      paste0(safe_metric_name, "_Repeated_Measures_Plot.png"))

    pdf_path <- file.path(
      plot_directory,
      paste0(safe_metric_name, "_Repeated_Measures_Plot.pdf"))

    ggplot2::ggsave(
      filename = png_path,
      plot = outcome_plot,
      dpi = 720,
      width = 8,
      height = 6,
      bg = "white")

    ggplot2::ggsave(
      filename = pdf_path,
      plot = outcome_plot,
      width = 8,
      height = 6,
      bg = "white",
      device = "pdf")

    list(
      emmeans = emm_results,
      posthocs = posthoc_results,
      plot_png = png_path,
      plot_pdf = pdf_path,
      error = NA_character_)
  }, error = function(error_condition) {
    posthoc_error <<- conditionMessage(error_condition)

    list(
      emmeans = tibble::tibble(),
      posthocs = tibble::tibble(),
      plot_png = NA_character_,
      plot_pdf = NA_character_,
      error = posthoc_error)
  })

  posthoc_output
}

#4.5 Run the full 20-outcome repeated measures analysis for one dataset
run_dataset_analysis <- function(
    analysis_data,
    dataset_key,
    dataset_label,
    results_directory) {

  message("Running repeated measures models for: ", dataset_label)

  model_results <- purrr::pmap(
    metric_manifest,
    function(analysis_order, column_name, metric_label) {
      fit_repeated_measures_model(
        analysis_data = analysis_data,
        dataset_key = dataset_key,
        dataset_label = dataset_label,
        dv = column_name,
        metric_label = metric_label,
        analysis_order = analysis_order)
    })

  names(model_results) <- metric_manifest$column_name

  interaction_results <- purrr::map_dfr(
    model_results,
    "interaction") %>%
    dplyr::arrange(analysis_order) %>%
    dplyr::mutate(
      p_fdr_20_outcomes = stats::p.adjust(
        p_value,
        method = "fdr"),
      nominal_interaction = !is.na(p_value) & p_value < 0.05,
      fdr_significant_interaction =
        !is.na(p_fdr_20_outcomes) & p_fdr_20_outcomes < 0.05)

  anova_results <- purrr::map_dfr(
    model_results,
    "anova") %>%
    dplyr::arrange(analysis_order)

  fixed_effect_results <- purrr::map_dfr(
    model_results,
    "fixed_effects") %>%
    dplyr::arrange(analysis_order)

  diagnostic_results <- purrr::map_dfr(
    model_results,
    "diagnostics") %>%
    dplyr::arrange(analysis_order)

  fitted_models <- purrr::map(
    model_results,
    "model")

  nominal_outcomes <- interaction_results %>%
    dplyr::filter(nominal_interaction) %>%
    dplyr::arrange(p_value)

  posthoc_results_list <- purrr::map(
    nominal_outcomes$DV,
    function(current_dv) {
      current_row <- nominal_outcomes %>%
        dplyr::filter(DV == current_dv) %>%
        dplyr::slice(1)

      run_nominal_posthocs(
        fitted_model = fitted_models[[current_dv]],
        analysis_data = analysis_data,
        dataset_key = dataset_key,
        dataset_label = dataset_label,
        dv = current_dv,
        metric_label = current_row$DV_Summary,
        analysis_order = current_row$analysis_order,
        plot_directory = file.path(results_directory, "plots"))
    })

  names(posthoc_results_list) <- nominal_outcomes$DV

  emmeans_results <- purrr::map_dfr(
    posthoc_results_list,
    "emmeans")

  posthoc_results <- purrr::map_dfr(
    posthoc_results_list,
    "posthocs")

  posthoc_error_results <- tibble::tibble(
    DV = names(posthoc_results_list),
    error = purrr::map_chr(
      posthoc_results_list,
      ~ ifelse(is.na(.x$error), NA_character_, .x$error))) %>%
    dplyr::left_join(
      metric_manifest %>%
        dplyr::select(
          DV = column_name,
          analysis_order,
          DV_Summary = metric_label),
      by = "DV") %>%
    dplyr::mutate(
      dataset = dataset_key,
      dataset_label = dataset_label,
      .before = 1) %>%
    dplyr::arrange(analysis_order)

  readr::write_csv(
    interaction_results,
    file.path(results_directory, "interaction_results_20_outcomes.csv"))

  readr::write_csv(
    anova_results,
    file.path(results_directory, "anova_all_terms.csv"))

  readr::write_csv(
    fixed_effect_results,
    file.path(results_directory, "fixed_effects.csv"))

  readr::write_csv(
    diagnostic_results,
    file.path(results_directory, "model_diagnostics.csv"))

  readr::write_csv(
    nominal_outcomes,
    file.path(results_directory, "nominally_significant_interactions.csv"))

  readr::write_csv(
    emmeans_results,
    file.path(results_directory, "emmeans_nominal_outcomes.csv"))

  readr::write_csv(
    posthoc_results,
    file.path(results_directory, "posthocs_nominal_outcomes.csv"))

  readr::write_csv(
    posthoc_error_results,
    file.path(results_directory, "posthoc_errors.csv"))

  writexl::write_xlsx(
    list(
      interaction = interaction_results,
      nominal_interactions = nominal_outcomes,
      anova_terms = anova_results,
      fixed_effects = fixed_effect_results,
      diagnostics = diagnostic_results,
      emmeans_nominal = emmeans_results,
      posthocs_nominal = posthoc_results,
      posthoc_errors = posthoc_error_results),
    path = file.path(
      results_directory,
      "repeated_measures_model_results.xlsx"))

  saveRDS(
    fitted_models,
    file = file.path(
      results_directory,
      "fitted_repeated_measures_models.rds"))

  message(
    dataset_label,
    ": ",
    sum(interaction_results$nominal_interaction, na.rm = TRUE),
    " nominally significant interactions and ",
    sum(
      interaction_results$fdr_significant_interaction,
      na.rm = TRUE),
    " FDR-significant interactions.")

  list(
    interaction = interaction_results,
    nominal_interactions = nominal_outcomes,
    anova = anova_results,
    fixed_effects = fixed_effect_results,
    diagnostics = diagnostic_results,
    emmeans = emmeans_results,
    posthocs = posthoc_results,
    posthoc_errors = posthoc_error_results,
    models = fitted_models)
}


## Demographic and Clinical Data Preparation ##

#5. Read and recode baseline race and ethnicity data
ethnicity_data <- read.delim(
  ethnicity_data_path,
  stringsAsFactors = FALSE) %>%
  dplyr::filter(eventname == "baseline_year_1_arm_1") %>%
  dplyr::transmute(
    subjectkey = trimws(as.character(subjectkey)),
    race_ethnicity = dplyr::case_when(
      as.character(race_ethnicity) == "1" ~ "White",
      as.character(race_ethnicity) == "2" ~ "Black",
      as.character(race_ethnicity) == "3" ~ "Hispanic",
      as.character(race_ethnicity) == "4" ~ "Asian",
      as.character(race_ethnicity) == "5" ~ "Other",
      is.na(race_ethnicity) |
        trimws(as.character(race_ethnicity)) == "" ~ "Unsure",
      TRUE ~ "Unsure"))

ethnicity_duplicates <- ethnicity_data %>%
  dplyr::count(subjectkey) %>%
  dplyr::filter(n > 1)

if (nrow(ethnicity_duplicates) > 0) {
  stop("The baseline race and ethnicity data contain duplicated subjectkeys.")
}

#5.1 Read and clean CBCL data at baseline and two-year follow-up
cbcl_data <- read.csv(
  cbcl_data_path,
  stringsAsFactors = FALSE) %>%
  dplyr::slice(-1) %>%
  dplyr::mutate(
    src_subject_id = trimws(as.character(src_subject_id)),
    eventname = trimws(as.character(eventname)),
    dplyr::across(
      c(
        cbcl_scr_syn_internal_t,
        cbcl_scr_syn_external_t,
        cbcl_scr_syn_anxdep_t,
        cbcl_scr_dsm5_anxdisord_t,
        cbcl_scr_dsm5_anxdisord_nm,
        cbcl_scr_syn_anxdep_nm,
        cbcl_scr_syn_internal_nm,
        cbcl_scr_syn_external_nm),
      ~ suppressWarnings(as.numeric(.x)))) %>%
  dplyr::filter(
    eventname %in% analysis_events,
    !is.na(cbcl_scr_dsm5_anxdisord_t),
    !is.na(cbcl_scr_syn_anxdep_t),
    !is.na(cbcl_scr_syn_internal_t),
    !is.na(cbcl_scr_syn_external_t),
    cbcl_scr_dsm5_anxdisord_nm == 0,
    cbcl_scr_syn_anxdep_nm == 0,
    cbcl_scr_syn_internal_nm == 0,
    cbcl_scr_syn_external_nm == 0) %>%
  dplyr::select(
    src_subject_id,
    eventname,
    cbcl_scr_syn_internal_t,
    cbcl_scr_syn_external_t,
    cbcl_scr_syn_anxdep_t,
    cbcl_scr_dsm5_anxdisord_t)

cbcl_duplicates <- cbcl_data %>%
  dplyr::count(src_subject_id, eventname) %>%
  dplyr::filter(n > 1)

if (nrow(cbcl_duplicates) > 0) {
  stop("The cleaned CBCL data contain duplicated subject-event rows.")
}

#5.2 Create and write demographic summaries for one analysis dataset
create_demographic_outputs <- function(
    analysis_data,
    dataset_key,
    dataset_label,
    results_directory) {

  demographic_data <- analysis_data %>%
    dplyr::mutate(
      subjectkey = as.character(subjectkey),
      eventname = as.character(eventname),
      group = as.character(group),
      sex = as.character(sex),
      site_name = as.character(site_name),
      scanner_model = dplyr::case_when(
        startsWith(site_name, "S") ~ "Siemens",
        startsWith(site_name, "P") ~ "Philips Healthcare",
        startsWith(site_name, "G") ~ "GE Healthcare",
        TRUE ~ NA_character_)) %>%
    dplyr::left_join(
      ethnicity_data,
      by = "subjectkey") %>%
    dplyr::left_join(
      cbcl_data,
      by = c(
        "subjectkey" = "src_subject_id",
        "eventname" = "eventname"))

  subject_level_data <- demographic_data %>%
    dplyr::arrange(subjectkey, eventname) %>%
    dplyr::group_by(subjectkey, group) %>%
    dplyr::summarise(
      sex = dplyr::first(sex),
      race_ethnicity = dplyr::first(race_ethnicity),
      .groups = "drop")

  grouped_subject_summary <- subject_level_data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n_subjects = dplyr::n(),
      n_female = sum(sex == "F", na.rm = TRUE),
      percent_female = 100 * n_female / n_subjects,
      n_white = sum(race_ethnicity == "White", na.rm = TRUE),
      percent_white = 100 * n_white / n_subjects,
      n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
      percent_black = 100 * n_black / n_subjects,
      n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
      percent_hispanic = 100 * n_hispanic / n_subjects,
      n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
      percent_asian = 100 * n_asian / n_subjects,
      n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
      percent_other = 100 * n_other / n_subjects,
      n_unsure = sum(race_ethnicity == "Unsure", na.rm = TRUE),
      percent_unsure = 100 * n_unsure / n_subjects,
      .groups = "drop")

  whole_subject_summary <- subject_level_data %>%
    dplyr::summarise(
      group = "Whole Sample",
      n_subjects = dplyr::n(),
      n_female = sum(sex == "F", na.rm = TRUE),
      percent_female = 100 * n_female / n_subjects,
      n_white = sum(race_ethnicity == "White", na.rm = TRUE),
      percent_white = 100 * n_white / n_subjects,
      n_black = sum(race_ethnicity == "Black", na.rm = TRUE),
      percent_black = 100 * n_black / n_subjects,
      n_hispanic = sum(race_ethnicity == "Hispanic", na.rm = TRUE),
      percent_hispanic = 100 * n_hispanic / n_subjects,
      n_asian = sum(race_ethnicity == "Asian", na.rm = TRUE),
      percent_asian = 100 * n_asian / n_subjects,
      n_other = sum(race_ethnicity == "Other", na.rm = TRUE),
      percent_other = 100 * n_other / n_subjects,
      n_unsure = sum(race_ethnicity == "Unsure", na.rm = TRUE),
      percent_unsure = 100 * n_unsure / n_subjects)

  subject_summary <- dplyr::bind_rows(
    grouped_subject_summary,
    whole_subject_summary) %>%
    dplyr::mutate(
      dataset = dataset_key,
      dataset_label = dataset_label,
      .before = 1)

  grouped_visit_summary <- demographic_data %>%
    dplyr::group_by(group, eventname) %>%
    dplyr::summarise(
      n_subjects = dplyr::n_distinct(subjectkey),
      mean_age = mean(age_in_years, na.rm = TRUE),
      sd_age = stats::sd(age_in_years, na.rm = TRUE),
      mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
      sd_fd = stats::sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
      mean_internalizing_cbcl = mean(
        cbcl_scr_syn_internal_t,
        na.rm = TRUE),
      sd_internalizing_cbcl = stats::sd(
        cbcl_scr_syn_internal_t,
        na.rm = TRUE),
      mean_externalizing_cbcl = mean(
        cbcl_scr_syn_external_t,
        na.rm = TRUE),
      sd_externalizing_cbcl = stats::sd(
        cbcl_scr_syn_external_t,
        na.rm = TRUE),
      mean_anxdep_cbcl = mean(
        cbcl_scr_syn_anxdep_t,
        na.rm = TRUE),
      sd_anxdep_cbcl = stats::sd(
        cbcl_scr_syn_anxdep_t,
        na.rm = TRUE),
      mean_dsm5_anx_cbcl = mean(
        cbcl_scr_dsm5_anxdisord_t,
        na.rm = TRUE),
      sd_dsm5_anx_cbcl = stats::sd(
        cbcl_scr_dsm5_anxdisord_t,
        na.rm = TRUE),
      n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
      percent_ge_scanner = 100 * n_ge_scanner / n_subjects,
      n_philips_scanner = sum(
        scanner_model == "Philips Healthcare",
        na.rm = TRUE),
      percent_philips_scanner = 100 * n_philips_scanner / n_subjects,
      n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
      percent_siemens_scanner = 100 * n_siemens_scanner / n_subjects,
      .groups = "drop")

  whole_visit_summary <- demographic_data %>%
    dplyr::group_by(eventname) %>%
    dplyr::summarise(
      group = "Whole Sample",
      n_subjects = dplyr::n_distinct(subjectkey),
      mean_age = mean(age_in_years, na.rm = TRUE),
      sd_age = stats::sd(age_in_years, na.rm = TRUE),
      mean_fd = mean(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
      sd_fd = stats::sd(rsfmri_c_ngd_meanmotion, na.rm = TRUE),
      mean_internalizing_cbcl = mean(
        cbcl_scr_syn_internal_t,
        na.rm = TRUE),
      sd_internalizing_cbcl = stats::sd(
        cbcl_scr_syn_internal_t,
        na.rm = TRUE),
      mean_externalizing_cbcl = mean(
        cbcl_scr_syn_external_t,
        na.rm = TRUE),
      sd_externalizing_cbcl = stats::sd(
        cbcl_scr_syn_external_t,
        na.rm = TRUE),
      mean_anxdep_cbcl = mean(
        cbcl_scr_syn_anxdep_t,
        na.rm = TRUE),
      sd_anxdep_cbcl = stats::sd(
        cbcl_scr_syn_anxdep_t,
        na.rm = TRUE),
      mean_dsm5_anx_cbcl = mean(
        cbcl_scr_dsm5_anxdisord_t,
        na.rm = TRUE),
      sd_dsm5_anx_cbcl = stats::sd(
        cbcl_scr_dsm5_anxdisord_t,
        na.rm = TRUE),
      n_ge_scanner = sum(scanner_model == "GE Healthcare", na.rm = TRUE),
      percent_ge_scanner = 100 * n_ge_scanner / n_subjects,
      n_philips_scanner = sum(
        scanner_model == "Philips Healthcare",
        na.rm = TRUE),
      percent_philips_scanner = 100 * n_philips_scanner / n_subjects,
      n_siemens_scanner = sum(scanner_model == "Siemens", na.rm = TRUE),
      percent_siemens_scanner = 100 * n_siemens_scanner / n_subjects,
      .groups = "drop")

  visit_summary <- dplyr::bind_rows(
    grouped_visit_summary,
    whole_visit_summary) %>%
    dplyr::mutate(
      dataset = dataset_key,
      dataset_label = dataset_label,
      .before = 1)

  subject_summary_wide <- subject_summary %>%
    dplyr::select(-dataset, -dataset_label) %>%
    tidyr::pivot_longer(
      cols = -group,
      names_to = "Variable",
      values_to = "Value") %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = Value) %>%
    dplyr::select(Variable, dplyr::everything())

  visit_summary_wide <- visit_summary %>%
    dplyr::select(-dataset, -dataset_label) %>%
    tidyr::pivot_longer(
      cols = -c(group, eventname),
      names_to = "Variable",
      values_to = "Value") %>%
    dplyr::mutate(
      group_event = paste(group, eventname, sep = "__")) %>%
    dplyr::select(-group, -eventname) %>%
    tidyr::pivot_wider(
      names_from = group_event,
      values_from = Value) %>%
    dplyr::select(Variable, dplyr::everything())

  readr::write_csv(
    subject_summary,
    file.path(results_directory, "demographics_subject_level_long.csv"))

  readr::write_csv(
    visit_summary,
    file.path(results_directory, "demographics_visit_level_long.csv"))

  readr::write_csv(
    subject_summary_wide,
    file.path(results_directory, "demographics_subject_level_wide.csv"))

  readr::write_csv(
    visit_summary_wide,
    file.path(results_directory, "demographics_visit_level_wide.csv"))

  writexl::write_xlsx(
    list(
      subject_long = subject_summary,
      visit_long = visit_summary,
      subject_wide = subject_summary_wide,
      visit_wide = visit_summary_wide),
    path = file.path(
      results_directory,
      "repeated_measures_demographics.xlsx"))

  list(
    subject_summary = subject_summary,
    visit_summary = visit_summary,
    subject_summary_wide = subject_summary_wide,
    visit_summary_wide = visit_summary_wide)
}


## Read and Validate Analysis Datasets ##

#6. Prepare the primary QC-complete analysis dataset
primary_data <- prepare_analysis_data(
  data_path = primary_data_path,
  dataset_name = "Primary QC-complete repeated measures dataset",
  expected_n_subjects = 1767L,
  expected_n_rows = 3534L,
  expected_group_counts = expected_primary_group_counts)

#6.1 Prepare the exact published-cohort sensitivity dataset
sensitivity_data <- prepare_analysis_data(
  data_path = sensitivity_data_path,
  dataset_name = "Published-cohort corrected sensitivity dataset",
  expected_n_subjects = 2074L,
  expected_n_rows = 4148L,
  expected_group_counts = expected_sensitivity_group_counts)


## Run Repeated Measures Analyses ##

#7.1 Run all 20 repeated measures models in the primary QC-complete cohort
primary_analysis_results <- run_dataset_analysis(
  analysis_data = primary_data,
  dataset_key = "primary",
  dataset_label = "Primary QC-complete cohort",
  results_directory = primary_results_dir)

#7.2 Run all 20 repeated measures models in the published-cohort sensitivity sample
sensitivity_analysis_results <- run_dataset_analysis(
  analysis_data = sensitivity_data,
  dataset_key = "published_cohort_sensitivity",
  dataset_label = "Published-cohort corrected sensitivity",
  results_directory = sensitivity_results_dir)


## Compare Primary and Sensitivity Results ##

#8. Join the primary and sensitivity omnibus interaction results
interaction_comparison <- primary_analysis_results$interaction %>%
  dplyr::select(
    analysis_order,
    DV,
    DV_Summary,
    primary_f_value = f_value,
    primary_numerator_df = numerator_df,
    primary_denominator_df = denominator_df,
    primary_p_value = p_value,
    primary_p_fdr = p_fdr_20_outcomes,
    primary_nominal = nominal_interaction,
    primary_fdr_significant = fdr_significant_interaction) %>%
  dplyr::full_join(
    sensitivity_analysis_results$interaction %>%
      dplyr::select(
        analysis_order,
        DV,
        DV_Summary,
        sensitivity_f_value = f_value,
        sensitivity_numerator_df = numerator_df,
        sensitivity_denominator_df = denominator_df,
        sensitivity_p_value = p_value,
        sensitivity_p_fdr = p_fdr_20_outcomes,
        sensitivity_nominal = nominal_interaction,
        sensitivity_fdr_significant = fdr_significant_interaction),
    by = c(
      "analysis_order",
      "DV",
      "DV_Summary")) %>%
  dplyr::arrange(analysis_order) %>%
  dplyr::mutate(
    same_nominal_conclusion =
      primary_nominal == sensitivity_nominal,
    same_fdr_conclusion =
      primary_fdr_significant == sensitivity_fdr_significant)

readr::write_csv(
  interaction_comparison,
  file.path(
    comparison_results_dir,
    "primary_vs_published_cohort_interaction_comparison.csv"))

writexl::write_xlsx(
  list(
    interaction_comparison = interaction_comparison,
    primary_interactions = primary_analysis_results$interaction,
    sensitivity_interactions = sensitivity_analysis_results$interaction,
    primary_posthocs = primary_analysis_results$posthocs,
    sensitivity_posthocs = sensitivity_analysis_results$posthocs),
  path = file.path(
    comparison_results_dir,
    "primary_vs_published_cohort_results.xlsx"))


## Generate Demographic and Clinical Summaries ##

#9. Generate demographic summaries for the primary QC-complete cohort
primary_demographic_results <- create_demographic_outputs(
  analysis_data = primary_data,
  dataset_key = "primary",
  dataset_label = "Primary QC-complete cohort",
  results_directory = primary_results_dir)

#9.1 Generate demographic summaries for the published-cohort sensitivity sample
sensitivity_demographic_results <- create_demographic_outputs(
  analysis_data = sensitivity_data,
  dataset_key = "published_cohort_sensitivity",
  dataset_label = "Published-cohort corrected sensitivity",
  results_directory = sensitivity_results_dir)


## Final Output Summary ##

#10. Print the primary and sensitivity interaction results
print(primary_analysis_results$interaction, n = Inf)
print(sensitivity_analysis_results$interaction, n = Inf)

#10.1 Print completion messages
message(
  "Repeated measures analysis completed successfully. Results are in: ",
  results_root)
