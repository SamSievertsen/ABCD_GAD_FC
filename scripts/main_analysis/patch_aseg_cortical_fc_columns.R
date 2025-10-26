## Set-Up ## 

# Load packages for loading, wrangling, and mutating data
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Set path to (corrected) column mapping
column_mapping_path <- "./data_raw/rsfmri_gpnet_aseg_correction.csv"

# Set path to raw rsfMRI data
ABCD_rsfMRI_Data_path <- "./data_raw/ABCD_rsfMRI_Data.csv"

# Set path to backup of previous version of subset data (unmodified and used as a reference for what corrected imaging data should be retained)
base_subset_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_BACKUP_before_patch.csv"

# Final corrected subset to write
patched_out_path <- "./data_processed/main_analysis/subset_qcd_imaging_data.csv"

# Validation report path to write
validation_report_path <- "./data_processed/main_analysis/gpnet_aseg_rename_validation_report.csv"

# Backup the previous subset output if it exists
if (file.exists(patched_out_path)) {
  file.copy(
    patched_out_path,
    sub("\\.csv$", paste0("_backup_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"), patched_out_path),
    overwrite = TRUE)
  message("Backed up previous output: ", patched_out_path)
}


## Patch imaging data column names ##

# Read in the mapping of old to correct column names
map_raw <- read_csv(column_mapping_path, show_col_types = FALSE)

# Define mapping column headers and assert presence of variables
col_prev <- "name_nda (previously used)"
col_corr <- "name_nda (correct)"

# Normalize and sanitize the mapping table
mapping <- map_raw %>%
  dplyr::select(`name_nda (previously used)`, `name_nda (correct)`) %>%
  filter(!is.na(`name_nda (previously used)`), !is.na(`name_nda (correct)`)) %>%
  distinct() %>%
  rename(old = `name_nda (previously used)`,
    correct = `name_nda (correct)`) %>%
  mutate(old = trimws(old), correct = trimws(correct))

# Enforce one to one transformation on old column names
dup_old <- mapping %>% count(old) %>% filter(n > 1)
if (nrow(dup_old) > 0) stop("Mapping has duplicate OLD names. Fix the CSV before proceeding.")

# Read raw imaging table and remove the variable description row
ABCD_rsfMRI_Data <- read_csv(ABCD_rsfMRI_Data_path, show_col_types = FALSE) %>% slice(-1)

# Define key columns and assert their presence in ABCD_rsfMRI_Data
key_cols <- c("subjectkey", "eventname")
if (!all(key_cols %in% names(ABCD_rsfMRI_Data))) stop("ABCD_rsfMRI_Data is missing subjectkey/eventname.")

# Assert uniqueness of subject ID's in ABCD_rsfMRI_Data
dups <- ABCD_rsfMRI_Data %>% count(across(all_of(key_cols))) %>% filter(n > 1)
if (nrow(dups) > 0) stop("ABCD_rsfMRI_Data has duplicate subjectkey+eventname rows; cannot safely join.")

# Put ABCD_rsfMRI_Data into correct name space using the mapping dataframe
old_to_correct <- setNames(mapping$correct, mapping$old)
nm <- old_to_correct[names(ABCD_rsfMRI_Data)]
names(ABCD_rsfMRI_Data)[!is.na(nm)] <- nm[!is.na(nm)]

# Read in the old subset of imaging data to derive the column names we want to keep
base <- read_csv(base_subset_path, show_col_types = FALSE)
if (!all(key_cols %in% names(base))) stop("Previous imaging data subset is missing subjectkey/eventname.")

# Restrict mapping to columns present in the old imaging data subset
map_used <- mapping %>% filter(old %in% names(base))

# Split mapping into truly misnamed and identity rows
misnamed_rows <- map_used %>% filter(old != correct)
identity_rows <- map_used %>% filter(old == correct)

# Collect misnamed old headers and all target correct headers
misnamed_olds_present <- misnamed_rows$old
targets_needed <- map_used$correct %>% unique()

# Assert ABCD_rsfMRI_Data contains all target correct headers
missing_targets <- setdiff(targets_needed, names(ABCD_rsfMRI_Data))
if (length(missing_targets)) {
  stop("These CORRECT targets are missing in ABCD_rsfMRI_Data (after header standardization): ",
    paste(missing_targets, collapse = ", "))
}

# Dynamic, informative message about the planned change operations
message(
  "Will drop ", length(misnamed_olds_present), " misnamed columns and refresh ",
    nrow(identity_rows), " identity columns from ABCD_rsfMRI_Data.")

# Build the patched subset by removing misnamed olds and preexisting targets then joining correct targets from ABCD_rsfMRI_Data
patched <- base %>%
  dplyr::select(-all_of(misnamed_olds_present)) %>%
  dplyr::select(-any_of(targets_needed)) %>%
  left_join(ABCD_rsfMRI_Data %>% select(all_of(c(key_cols, targets_needed))), by = key_cols)

# Prepare alignment index for validation of changes made to column names
row_idx <- match(paste(base$subjectkey, base$eventname),
            paste(ABCD_rsfMRI_Data$subjectkey, ABCD_rsfMRI_Data$eventname))

# HCreate a helper for numeric coercion and tolerance for parity
to_num <- function(x) suppressWarnings(as.numeric(x))
tol <- 1e-10

# Define a list to accumulate validation rows for reporting
val_rows <- list()

# Validate each target column name by comparing patched to ABCD_rsfMRI_Data and to one old reference where applicable
for (cor_nm in targets_needed) {
  olds_for_target <- misnamed_rows %>% filter(correct == cor_nm) %>% pull(old)
  old_ref <- if (length(olds_for_target)) olds_for_target[1] else NA_character_
  
  new_num <- to_num(patched[[cor_nm]])
  ABCD_rsfMRI_Data_num <- to_num(ABCD_rsfMRI_Data[[cor_nm]][row_idx])
  
  if (!is.na(old_ref) && old_ref %in% names(base)) {
    old_num <- to_num(base[[old_ref]])
    n_equal_old_vs_ABCD_rsfMRI_Data <- sum(abs(old_num - ABCD_rsfMRI_Data_num) <= tol | (is.na(old_num) & is.na(ABCD_rsfMRI_Data_num)), na.rm = TRUE)
    n_changed <- sum((abs(old_num - new_num) > tol) | xor(is.na(old_num), is.na(new_num)), na.rm = TRUE)
  } else {
    n_equal_old_vs_ABCD_rsfMRI_Data <- NA_integer_
    n_changed <- NA_integer_
  }
  
  n_equal_new_vs_ABCD_rsfMRI_Data <- sum(abs(new_num - ABCD_rsfMRI_Data_num) <= tol | (is.na(new_num) & is.na(ABCD_rsfMRI_Data_num)), na.rm = TRUE)
  
  val_rows[[length(val_rows) + 1]] <- tibble::tibble(
    target_correct_name = cor_nm,
    one_old_example = old_ref,
    n_rows = nrow(patched),
    n_equal_new_vs_ABCD_rsfMRI_Data = n_equal_new_vs_ABCD_rsfMRI_Data,
    n_equal_old_vs_ABCD_rsfMRI_Data = n_equal_old_vs_ABCD_rsfMRI_Data,
    n_values_changed_old_to_new = n_changed
  )
}

# Write the validation report outlining which columns were changed and to what
write_csv(bind_rows(val_rows), validation_report_path)


## Integrity checks that only flag truly misnamed leftovers ##

# Some old names also appear as correct targets and are allowed after the join; define what those are
allowable_old_names <- intersect(misnamed_olds_present, targets_needed)

# Define true misnamed olds that should not remain
strict_olds <- setdiff(misnamed_olds_present, targets_needed)

# Detect any strict leftovers still present
leftover_strict <- intersect(names(patched), strict_olds)

# Fail if any strict misnamed columns remain
if (length(leftover_strict)) {
  stop("These misnamed OLD columns remain but should have been removed: ",
    paste(leftover_strict, collapse = ", "))
}

# Note benign cases that are both old in some rows and correct in others
benign <- intersect(names(patched), allowable_old_names)
if (length(benign)) {
  message("Note: these names appear as old in the mapping but are kept because they are also correct targets: ",
    paste(benign, collapse = ", "))
}

# Save the corrected subset of data once validation checks show passable output
write_csv(patched, patched_out_path)

# Perform an additional spot check A that shows that the corrected column now matches ABCD_rsfMRI_Data and corresponds to the old mislabel from the old subset backup file
probeA_old <- "rsfmri_cor_ngd_cerc_scs_thprh"
probeA_new <- "rsfmri_cor_ngd_dsa_scs_cdelh"
if (probeA_new %in% names(patched)) {
  message("Spot-check A (", probeA_old, " -> ", probeA_new, "): first 10 rows")
  print(
    tibble::tibble(
      subjectkey = base$subjectkey,
      eventname = base$eventname,
      old_value = if (probeA_old %in% names(base)) base[[probeA_old]] else NA_real_,
      new_value = patched[[probeA_new]],
      ABCD_rsfMRI_Data_value  = ABCD_rsfMRI_Data[[probeA_new]][row_idx]) %>% 
      dplyr::slice_head(n = 10)
  )
}

# Perform an additional spot check B for another explicit pair
probeB_old <- "rsfmri_cor_ngd_fopa_scs_crcxrh"
probeB_new <- "rsfmri_cor_ngd_cerc_scs_aalh"
if (probeB_new %in% names(patched)) {
  message("Spot-check B (", probeB_old, " -> ", probeB_new, "): first 10 rows")
  print(
    tibble::tibble(
      subjectkey = base$subjectkey,
      eventname = base$eventname,
      old_value = if (probeB_old %in% names(base)) base[[probeB_old]] else NA_real_,
      new_value = patched[[probeB_new]],
      ABCD_rsfMRI_Data_value = ABCD_rsfMRI_Data[[probeB_new]][row_idx]
    ) %>% dplyr::slice_head(n = 10)
  )
}
