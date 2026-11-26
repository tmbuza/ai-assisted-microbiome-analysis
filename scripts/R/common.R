# Shared helpers. Run entry-point scripts from the folder containing _quarto.yml.
fail <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
read_text <- function(path) {
  fail(file.exists(path), paste("Missing file:", path))
  x <- read.delim(path, check.names = FALSE, colClasses = "character",
                  na.strings = c("NA", ""), quote = "", comment.char = "")
  fail(!anyDuplicated(names(x)) && all(nzchar(trimws(names(x)))),
       paste("Duplicate or empty column names:", path))
  x
}
valid_ids <- function(x) !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
read_inputs <- function(input) {
  x <- read_text(file.path(input, "feature-counts.tsv"))
  m <- read_text(file.path(input, "sample-metadata.tsv"))
  tx <- read_text(file.path(input, "taxonomy.tsv"))
  fail(ncol(x) >= 2 && names(x)[1] == "feature_id" && nrow(x) > 1,
       "Counts need feature_id and sample columns, with at least two features.")
  fail(all(c("sample_id", "group") %in% names(m)), "Metadata need sample_id and group.")
  fail(all(c("feature_id", "taxon") %in% names(tx)), "Taxonomy needs feature_id and taxon.")
  fail(valid_ids(x$feature_id) && valid_ids(m$sample_id) && valid_ids(tx$feature_id),
       "Identifiers must be unique, nonmissing and nonblank.")
  fail(!anyNA(m$group) && all(nzchar(trimws(m$group))), "Group labels are missing.")
  fail(!anyNA(tx$taxon), "Taxonomy labels are missing; use feature IDs if unassigned.")
  fail(setequal(names(x)[-1], m$sample_id), "Counts and metadata sample IDs differ.")
  fail(setequal(x$feature_id, tx$feature_id), "Counts and taxonomy feature IDs differ.")
  a <- as.matrix(x[-1])
  suppressWarnings(storage.mode(a) <- "double")
  fail(all(is.finite(a)) && all(a >= 0) && all(a == floor(a)),
       "Counts must be finite nonnegative integers. Do not round signals into counts.")
  fail(all(colSums(a) > 0), "Every sample needs a positive count total.")
  fail(ncol(a) >= 4 && length(unique(m$group)) >= 2,
       "This comparison needs at least four samples and two groups.")
  dimnames(a) <- list(x$feature_id, names(x)[-1])
  m <- m[match(colnames(a), m$sample_id), , drop = FALSE]
  rownames(m) <- m$sample_id
  tx <- tx[match(rownames(a), tx$feature_id), , drop = FALSE]
  list(a = a, m = m, tx = tx)
}
read_settings <- function(path) {
  z <- read_text(path)
  fail(identical(names(z), c("parameter", "value")) && valid_ids(z$parameter),
       "Settings need unique parameter/value rows.")
  cfg <- setNames(z$value, z$parameter)
  required <- c("reference", "comparison", "covariates", "design_confirmed",
                "run_da", "top_n", "seed", "depth", "depth_repeats", "permutations", "prevalence", "alpha")
  fail(all(required %in% names(cfg)), "Settings are missing required parameters.")
  fail(!anyNA(cfg[required]), "Required settings cannot be empty or missing.")
  fail(cfg[["design_confirmed"]] %in% c("TRUE", "FALSE") &&
       cfg[["run_da"]] %in% c("TRUE", "FALSE"), "Boolean settings must be TRUE or FALSE.")
  for (n in c( "top_n", "seed", "depth", "depth_repeats", "permutations", "prevalence", "alpha")) {
    v <- suppressWarnings(as.numeric(cfg[[n]]))
    fail(length(v) == 1 && is.finite(v), paste("Invalid setting:", n))
  }
  for (n in c( "top_n", "seed", "depth", "depth_repeats", "permutations"))
    fail(as.numeric(cfg[[n]]) == floor(as.numeric(cfg[[n]])), paste(n, "must be an integer."))
  fail(as.numeric(cfg[["top_n"]]) >= 1 && as.numeric(cfg[["depth"]]) >= 0 &&
       as.numeric(cfg[["depth_repeats"]]) >= 2 && as.numeric(cfg[["permutations"]]) >= 99,
       "Use positive top_n, nonnegative depth, >=2 repeats and >=99 permutations.")
  fail(as.numeric(cfg[["prevalence"]]) > 0 && as.numeric(cfg[["prevalence"]]) < 1 &&
       as.numeric(cfg[["alpha"]]) > 0 && as.numeric(cfg[["alpha"]]) < 1,
       "Prevalence and alpha must lie strictly between zero and one.")
  cfg
}
sha_inputs <- function(input, config) {
  paths <- c(file.path(input, c("feature-counts.tsv", "sample-metadata.tsv", "taxonomy.tsv")), config)
  # Base-R MD5 records exact bytes for freshness checks, not security.
  data.frame(file = normalizePath(paths), md5 = unname(tools::md5sum(paths)))
}
write_table <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE,
                                           quote = FALSE, na = "NA")
shannon <- function(a) {
  p <- sweep(a, 2, colSums(a), "/")
  q <- p; q[q == 0] <- 1
  -colSums(p * log(q))
}
feature_labels <- function(tx) {
  label <- vapply(strsplit(tx$taxon, ";", fixed = TRUE), function(z) {
    z <- trimws(sub("^[dkpcofgs]__|^D_[0-9]+__", "", trimws(z)))
    z <- z[nzchar(z)]
    if (length(z)) tail(z, 1) else "Unclassified"
  }, character(1))
  if ("taxon_label" %in% names(tx)) {
    supplied <- !is.na(tx$taxon_label) & nzchar(trimws(tx$taxon_label))
    label[supplied] <- tx$taxon_label[supplied]
  }
  # Full IDs prevent collisions when several OTUs share a genus assignment.
  setNames(paste0(label, " [", tx$feature_id, "]"), tx$feature_id)
}
entry_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  fail(length(args) == 3, "Usage: Rscript SCRIPT input_dir output_dir config.tsv")
  fail(file.exists("_quarto.yml"), "Run from the main project folder.")
  args
}
verify_fresh <- function(input, config, stored) {
  now <- sha_inputs(input, config)
  fail(identical(now$md5, stored$md5), "Inputs or settings changed. Rerun exploration first.")
}
model_support <- function(a, design, feature_ids, available) {
  positive <- a[feature_ids, , drop = FALSE] > 0
  ranks <- apply(positive, 1, function(z) qr(design[z, , drop = FALSE])$rank)
  data.frame(feature_id = feature_ids, positive_samples = rowSums(positive),
    positive_model_rank = ranks, design_columns = ncol(design),
    positive_support_full_rank = ranks == ncol(design), estimate_available = available)
}
