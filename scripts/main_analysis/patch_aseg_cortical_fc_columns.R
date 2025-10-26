## Set-Up ## 

# Load packages for loading, wrangling, and mutating data
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Set path to mapping file for old -> correct column names
column_mapping_path <- "./data_raw/rsfmri_gpnet_aseg_correction.csv"

# Set path to raw rsfMRI data
raw_path <- "./data_raw/ABCD_rsfMRI_Data.csv" 

# Set path to backup of imaging data to use as reference for selected columns/order
base_subset_path <- "./data_processed/main_analysis/subset_qcd_imaging_data_BACKUP_before_patch.csv"

# Set path to patched output (to write)
patched_out_path <- "./data_processed/main_analysis/subset_qcd_imaging_data.csv"

# Set path to audit report (to write)
validation_path <- "./data_processed/main_analysis/gpnet_aseg_rename_validation_report.csv"

# Back up the previous output if present with a timestamp
if (file.exists(patched_out_path)) {
  file.copy(
    patched_out_path,
    sub("\\.csv$", paste0("_backup_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv"), patched_out_path),
    overwrite = TRUE)
  message("Backed up previous output: ", patched_out_path)
}

# Read in mapping dataframe
map_raw <- read_csv(column_mapping_path, show_col_types = FALSE)

# Set column labels for old names
col_prev <- "name_nda (previously used)"

# Set column labels for correct names
col_corr <- "name_nda (correct)"                                     

# Clean mapping dataframe
mapping <- map_raw %>%
  dplyr::select(`name_nda (previously used)`, `name_nda (correct)`) %>%    
  filter(!is.na(`name_nda (previously used)`), !is.na(`name_nda (correct)`)) %>% 
  distinct() %>%
  rename(old = `name_nda (previously used)`,
         correct = `name_nda (correct)`) %>%
  mutate(old = trimws(old), correct = trimws(correct))

# Enforce unique OLD names
dup_old <- mapping %>% count(old) %>% filter(n > 1)       

# Guard against ambiguous map structure if duplicates exist after cleaning of mapping
if (nrow(dup_old) > 0) stop("Mapping has duplicate OLD names. Fix the CSV before proceeding.")  

# Build vector lookups for mapping old columns to new (correct) names
old_to_correct <- setNames(mapping$correct, mapping$old)
correct_set <- unique(mapping$correct)
old_set <- unique(mapping$old)

# Read in raw ABCD rsfMRI data and remove description row
raw <- read_csv(raw_path, show_col_types = FALSE) %>% slice(-1)

# Standardize subject IDs and timepoint indicators in raw ABCD rsfMRI data
raw <- raw %>%
  mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)))

# Assert & enforce unique merging keys in Raw ABCD rsfMRI data
raw_dups <- raw %>% count(subjectkey, eventname) %>% filter(n > 1)
if (nrow(raw_dups) > 0) stop("Raw ABCD rsfMRI data has duplicate subjectkey+eventname rows; cannot safely join.")  

# Rename Raw ABCD rsfMRI data headers from OLD to CORRECT once
nm <- old_to_correct[names(raw)] 
names(raw)[!is.na(nm)] <- nm[!is.na(nm)]

# Read in the base subset of imaging data that defines the column schema to use
base <- read_csv(base_subset_path, show_col_types = FALSE)

# Standardize & assert merging keys in base
if (!all(c("subjectkey","eventname") %in% names(base))) stop("Base is missing subjectkey/eventname.")
base <- base %>%
  mutate(
    subjectkey = trimws(as.character(subjectkey)),
    eventname = trimws(as.character(eventname)))

# Build merging key alignment vector
row_idx <- match(paste(base$subjectkey, base$eventname), paste(raw$subjectkey,  raw$eventname))

# Report merging key match rate to check for total rows in base data vs matched rows
n_total <- nrow(base) 
n_match <- sum(!is.na(row_idx))
message(sprintf("Key alignment: matched %d / %d rows (%.1f%%).",
                n_match, n_total, 100 * n_match / max(n_total, 1)))

# Initialize patched data with base subset
patched <- base

# Define columns to consider for overwrite with new correct names
base_cols <- setdiff(names(base), c("subjectkey","eventname"))

# Create a helper that converts to numeric quietly
to_num <- function(x) suppressWarnings(as.numeric(x))     

# Create a helper for numeric tolerance for parity across merged columns
tol <- 1e-10      

# Create a holder for per column validation records
val_rows <- list()      

# Create a counter for updated columns
updated_count <- 0L     

# Create a tracker for mapped columns missing in raw
skipped_missing_in_raw <- character()                                 


## Change Names of Incorrectly Labeled Data Where Relevant ##

# Overwrite each base column from Raw ABCD rsfMRI data using sequential merging logic
for (b in base_cols) {
  
  # If b is a correct name and exists in Raw ABCD rsfMRI data-after-rename -> use Raw ABCD rsfMRI data[b]
  if (b %in% correct_set && b %in% names(raw)) {
    src <- b
    src_reason <- "correct_name_in_raw"
    
    # Else if b is an old label and its correct target exists in Raw ABCD rsfMRI data -> use Raw ABCD rsfMRI data[mapped correct]
  } else if (b %in% names(old_to_correct) && old_to_correct[[b]] %in% names(raw)) {
    src <- old_to_correct[[b]]
    src_reason <- "mapped_correct"
    
    # Else if b happens to exist in Raw ABCD rsfMRI data (rare, but keep as fallback) -> use Raw ABCD rsfMRI data[b]
  } else if (b %in% names(raw)) {
    src <- b
    src_reason <- "identity_in_raw_fallback"
    
    # Else leave as-is (e.g., site_name, ksads_10_869_p)
  } else {
    if (b %in% names(old_to_correct) && !(old_to_correct[[b]] %in% names(raw))) {
      
      # Record mapped target missing in raw
      skipped_missing_in_raw <- c(skipped_missing_in_raw, b)         
    }
    
    # Skip overwrite for this column
    next                                                             
  }
  
  # Capture previous values for comparison
  old_vals <- patched[[b]]      
  
  # Pull aligned values from chosen raw source
  new_vals <- raw[[src]][row_idx] 
  
  # Write new values into patched
  patched[[b]] <- new_vals   
  
  # Increment update counter
  updated_count <- updated_count + 1L                                 
  
  # Numeric version of old values
  old_num <- to_num(old_vals)
  
  # Numeric version of new values
  new_num <- to_num(new_vals)     
  
  # Numeric version of raw source values
  raw_num <- to_num(raw[[src]][row_idx])      
  
  # Append per column validation summary
  val_rows[[length(val_rows) + 1L]] <- tibble::tibble(               
    output_column_name = b,
    raw_source_column_used = src,
    source_reason = src_reason,
    n_rows = nrow(patched),
    n_equal_new_vs_raw = sum(abs(new_num - raw_num) <= tol | (is.na(new_num) & is.na(raw_num)), na.rm = TRUE),
    n_equal_old_vs_raw = sum(abs(old_num - raw_num) <= tol | (is.na(old_num) & is.na(raw_num)), na.rm = TRUE),
    n_values_changed_old_to_new = sum((abs(old_num - new_num) > tol) | xor(is.na(old_num), is.na(new_num)), na.rm = TRUE))
}

# Log summary of updated columns
message(sprintf(
  "Updated %d/%d base columns from Raw ABCD rsfMRI data.",
  updated_count,
  length(base_cols)))


## Initial Spot Checks of Re-Naming Process ##

# Warn about any unmapped target columns during merging
if (length(skipped_missing_in_raw)) {
  warning(
    "These base columns are in the mapping but their correct targets are not in Raw ABCD rsfMRI data; left unchanged: ",
    paste(unique(skipped_missing_in_raw), collapse = ", "))                         
}

# Write validation report
if (length(val_rows)) {
  readr::write_csv(dplyr::bind_rows(val_rows), validation_path)
  message("Validation report written: ", validation_path)
} else {
  readr::write_csv(tibble::tibble(
    output_column_name = character(),
    raw_source_column_used = character(),
    n_rows = integer(),
    n_equal_new_vs_raw = integer(),
    n_equal_old_vs_raw = integer(),
    n_values_changed_old_to_new = integer()), 
    validation_path)
  message("Validation report written (no columns updated).")
}


## Incorporate Data Quality QC ##

# Read recommended imaging inclusion flags
recommended_imaging <- readr::read_csv("./data_raw/abcd_imgincl01.csv", show_col_types = FALSE)

# Drop the description row
recommended_imaging <- recommended_imaging[-1, ]

# Keep/normalize only required QC fields and standardize keys
recommended_imaging <- recommended_imaging %>%
  dplyr::select(src_subject_id, eventname, imgincl_rsfmri_include) %>%
  dplyr::mutate(
    src_subject_id = trimws(as.character(src_subject_id)),
    eventname = trimws(as.character(eventname)),
    imgincl_rsfmri_include = suppressWarnings(as.numeric(imgincl_rsfmri_include)))

# Attach QC flags to patched data
patched <- patched %>%
  dplyr::left_join(recommended_imaging, by = c("subjectkey" = "src_subject_id", "eventname" = "eventname"))

# Write final patched dataset
readr::write_csv(patched, patched_out_path)


## Integrity Checks ##

# Choose raw source for a base column using same rules in env
decide_src <- function(b, old_to_correct, raw_names, correct_set) {
  if (b %in% correct_set && b %in% raw_names) {
    list(src = b, reason = "correct_name_in_raw")
  } else if (b %in% names(old_to_correct) && old_to_correct[[b]] %in% raw_names) {
    list(src = old_to_correct[[b]], reason = "mapped_correct")
  } else if (b %in% raw_names) {
    list(src = b, reason = "identity_in_raw_fallback")
  } else {
    list(src = NA_character_, reason = "left_as_is")
  }
}

# Audit that the above logic matches the actual processing loop; first recompute list of data columns
base_cols <- setdiff(names(base), c("subjectkey","eventname"))

# Define tolerance and numeric equality helpers for audit
tol  <- 1e-10; to_num <- function(x) suppressWarnings(as.numeric(x))
all_equal_num <- function(a, b) all(abs(a - b) <= tol | (is.na(a) & is.na(b)), na.rm = TRUE)]

# Decide sources and audit the validity of them + the resulting data merged into the new df
src_info <- lapply(base_cols, function(b) decide_src(b, old_to_correct, names(raw), correct_set))
audit <- tibble::tibble(
  column = base_cols,
  src = vapply(src_info, function(x) x$src,    character(1)),
  reason = vapply(src_info, function(x) x$reason, character(1))) %>%
  
  # Flag columns that were updated & patched, with comparisons between src and putput (i.e., does patched == raw? How about base?)
  dplyr::mutate(
    updated = !is.na(src) & reason != "left_as_is",
    base_vals_num = purrr::map(column, ~ to_num(base[[.x]])),
    patched_vals_num = purrr::map(column, ~ to_num(patched[[.x]])),
    raw_vals_num = purrr::map(src,
                              ~ if (is.na(.x)) rep(NA_real_, nrow(base)) else to_num(raw[[.x]][row_idx])),
    eq_patched_raw = purrr::map2_lgl(patched_vals_num, raw_vals_num,  all_equal_num),
    eq_base_raw = purrr::map2_lgl(base_vals_num,    raw_vals_num,  all_equal_num),
    eq_patched_base = purrr::map2_lgl(patched_vals_num, base_vals_num, all_equal_num),
    n_changed = purrr::pmap_int(                                  
      list(base_vals_num, patched_vals_num),
      ~ sum(!(abs(..1 - ..2) <= tol | (is.na(..1) & is.na(..2))), na.rm = TRUE)),
    n_rows = nrow(base))

# Print summary of audit
print(audit %>% summarise(                                       
  n_base_cols = n(),
  n_updated = sum(updated),
  n_all_equal_new_vs_raw = sum(eq_patched_raw[updated]),
  n_any_values_changed = sum(n_changed > 0 & updated)))

# Check pairwise assertion that respects the same logic rules above
present_pairs <- mapping %>%
  
  # Pairs present in both places?
  dplyr::filter(old %in% names(base), old != correct, correct %in% names(raw)) %>%
  dplyr::mutate(
    src_used = vapply(old, function(b) decide_src(b, old_to_correct, names(raw), correct_set)$src, character(1)),
    
    # Numeric identity check to chosen source
    ok = mapply(function(o, s) {                                             
      a <- to_num(patched[[o]])
      b <- to_num(raw[[s]][row_idx])
      all(abs(a - b) <= tol | (is.na(a) & is.na(b)), na.rm = TRUE)
    }, old, src_used))

# Show any true mismatches (should be none)
print(present_pairs %>% dplyr::filter(!ok))

# For each correct metric present in base, list patched[CORRECT] == Raw ABCD rsfMRI data[CORRECT]
correct_in_base <- intersect(correct_set, names(base))

# Assert patched equals raw for correct columns
ok_correct <- vapply(correct_in_base, function(nm) {                
  a <- suppressWarnings(as.numeric(patched[[nm]]))
  b <- suppressWarnings(as.numeric(raw[[nm]][row_idx]))
  all(abs(a - b) <= 1e-10 | (is.na(a) & is.na(b)), na.rm = TRUE)
}, logical(1))

# Map output columns to sources
source_map <- audit %>% filter(updated) %>% select(column, src, reason)   

# Find any shared raw sources
collisions <- source_map %>% count(src) %>% filter(n > 1)    

# Print message with dynamic detail if collisions occurred during merging
if (nrow(collisions)) {
  message("Multiple output columns draw from the same Raw ABCD rsfMRI data source:")
  print(source_map %>% semi_join(collisions, by = "src") %>% arrange(src))
}
