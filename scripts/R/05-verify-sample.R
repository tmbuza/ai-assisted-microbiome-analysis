# Select the input, output, and sample to check.
input_dir <- "data/example"
output_dir <- "results"
sample_id <- "A01"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3L) {
  input_dir <- args[[1]]
  output_dir <- args[[2]]
  sample_id <- args[[3]]
} else if (length(args) != 0L) {
  stop("Usage: Rscript scripts/R/05-verify-sample.R [input_dir output_dir sample_id]",
       call. = FALSE)
}

# Read the original counts and the saved analysis results.
counts <- read.delim(
  file.path(input_dir, "feature-counts.tsv"),
  check.names = FALSE
)

reported <- read.delim(
  file.path(output_dir, "tables", "alpha-diversity.tsv"),
  check.names = FALSE
)

# Stop if the requested sample is missing or has multiple result rows.
stopifnot(sample_id %in% names(counts))
reported_sample <- reported[which(reported$sample_id == sample_id), , drop = FALSE]
stopifnot(nrow(reported_sample) == 1)

# Extract one sample's counts after the input-validation step.
sample_counts <- counts[[sample_id]]
stopifnot(
  is.numeric(sample_counts),
  !anyNA(sample_counts),
  all(is.finite(sample_counts)),
  all(sample_counts >= 0),
  all(sample_counts == floor(sample_counts)),
  sum(sample_counts) > 0
)

# 1. Add all feature counts.
feature_total <- sum(sample_counts)

# 2. Count features present in this sample.
observed_features <- sum(sample_counts > 0)

# 3. Express each count as a proportion of the sample total.
proportions <- sample_counts / feature_total

# Check that the proportions sum to one.
stopifnot(abs(sum(proportions) - 1) < 1e-10)

# 4. Exclude zeros before taking logarithms.
# Zero-count features contribute zero to Shannon diversity.
positive_proportions <- proportions[proportions > 0]

# 5. Calculate Shannon diversity using natural logarithms.
shannon <- -sum(
  positive_proportions * log(positive_proportions)
)

# Compare the manual values with the saved table.
manual <- c(feature_total, observed_features, shannon)
saved <- c(
  reported_sample$feature_total,
  reported_sample$observed_features,
  reported_sample$shannon
)

comparison <- data.frame(
  measure = c("feature_total", "observed_features", "shannon"),
  manual = manual,
  saved = saved,
  agrees = abs(manual - saved) < 1e-6
)

print(comparison, digits = 7, row.names = FALSE)
stopifnot(all(comparison$agrees))
message("Verification passed for sample ", sample_id, ".")
