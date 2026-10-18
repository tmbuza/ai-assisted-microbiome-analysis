# 1. Set up the session ----------------------------------------------------

# Run from the repository root. Both directory arguments are optional.
# Usage: Rscript scripts/R/01-check-inputs.R INPUT_DIRECTORY OUTPUT_DIRECTORY
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[1] else "data/example"
output_dir <- if (length(args) >= 2) args[2] else "results"

# 2. Define input-checking helpers -----------------------------------------

# Stop with a readable message when a required check fails.
check <- function(ok, message) {
  if (!isTRUE(ok)) {
    stop(message, call. = FALSE)
  }
}

# Read tab-separated inputs without silently repairing duplicate names.
read_input <- function(name) {
  path <- file.path(input_dir, name)
  check(file.exists(path), paste("Missing input:", path))

  read_tsv(path, name_repair = "minimal", show_col_types = FALSE)
}

# 3. Read the three input tables -------------------------------------------

# x: feature counts; m: sample metadata; t: feature taxonomy.
x <- read_input("feature-counts.tsv")
m <- read_input("sample-metadata.tsv")
t <- read_input("taxonomy.tsv")

# 4. Check table structure -------------------------------------------------

# Verify required columns and reject duplicate column names.
check(
  ncol(x) >= 2 && names(x)[1] == "feature_id",
  "Counts need feature_id followed by sample columns."
)
check(all(c("sample_id", "group") %in% names(m)), "Metadata need sample_id and group.")
check(
  all(c("feature_id", "taxon") %in% names(t)),
  "Taxonomy needs feature_id and taxon."
)
check(nrow(x) > 0 && nrow(m) > 0, "Input tables must contain rows.")
check(
  !anyDuplicated(names(x)) && !anyDuplicated(names(m)) && !anyDuplicated(names(t)),
  "Duplicate column names."
)

# 5. Check identifiers and labels ------------------------------------------

# IDs must be present, nonblank, and unique within each table.
valid_ids <- function(v) {
  !anyNA(v) &&
    all(nzchar(trimws(as.character(v)))) &&
    !anyDuplicated(v)
}

check(valid_ids(x$feature_id), "Feature IDs must be unique and nonblank.")
check(valid_ids(m$sample_id), "Sample IDs must be unique and nonblank.")
check(valid_ids(t$feature_id), "Taxonomy feature IDs must be unique and nonblank.")
check(
  !anyNA(m$group) && all(nzchar(trimws(as.character(m$group)))),
  "Groups must be nonblank."
)
check(
  !anyNA(t$taxon) && all(nzchar(trimws(as.character(t$taxon)))),
  "Taxon labels must be nonblank."
)

# 6. Match the tables ------------------------------------------------------

# Require the same sample and feature IDs; row order may differ.
check(
  setequal(names(x)[-1], m$sample_id),
  "Count sample IDs must match metadata exactly."
)
check(
  setequal(x$feature_id, t$feature_id),
  "Count feature IDs must match taxonomy exactly."
)

# 7. Check count values and group sizes ------------------------------------

# a is the numeric feature-by-sample matrix used by the analysis.
check(all(vapply(x[-1], is.numeric, logical(1))), "Counts must be numeric.")
a <- as.matrix(x[-1])
check(
  !anyNA(a) && all(is.finite(a)) && all(a >= 0) && all(a == floor(a)),
  "Counts must be finite nonnegative integers without missing values."
)
check(all(colSums(a) > 0), "Every sample must have a positive count total.")
check(n_distinct(m$group) >= 2, "Community comparison needs at least two groups.")

# These conditions warrant review but do not prevent the analysis from running.
if (any(rowSums(a) == 0)) {
  warning("All-zero features are present; document any filtering.")
}
if (any(table(m$group) < 3)) {
  warning("Some groups have fewer than three samples; summaries are limited.")
}

# 8. Save a compact validation report --------------------------------------

# Create the tables directory if it does not already exist.
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

write_tsv(
  tibble(
    check = c("samples", "features", "groups", "validation"),
    value = c(nrow(m), nrow(x), n_distinct(m$group), "PASS")
  ),
  file.path(output_dir, "tables", "input-checks.tsv")
)

message("Input checks passed.")
