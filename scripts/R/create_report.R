# Create an editable draft only once. Reports are independent of the book.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% c("synthetic", "zackular", "my-study"))
  stop("Usage: Rscript scripts/R/create_report.R synthetic|zackular|my-study", call. = FALSE)
if (!file.exists("_quarto.yml")) stop("Run from the main project folder.", call. = FALSE)
dataset <- args[1]
folder <- dataset
completion <- file.path("results", folder, "exploration-status.txt")
if (!file.exists(completion) || !identical(readLines(completion), "COMPLETE"))
  stop("Complete exploration before creating a report with its figures.", call. = FALSE)
source("scripts/R/common.R")
e <- readRDS(file.path("results", folder, "exploration.rds"))
verify_fresh(file.path("data", folder), file.path("config", paste0(dataset, ".tsv")), e$checksums)
label <- c(synthetic = "Synthetic microbiome practice", zackular = "Zackular human stool study",
            `my-study` = "My microbiome study")[[dataset]]
dir.create("reports", showWarnings = FALSE)
if (!file.exists("reports/_quarto.yml"))
  writeLines(c("project:", "  type: default", "  output-dir: ."), "reports/_quarto.yml")
path <- file.path("reports", paste0(dataset, "-report.qmd"))
if (file.exists(path)) {
  message("Existing draft kept: ", path)
  quit(save = "no", status = 0)
}
s <- readLines("templates/report.qmd", warn = FALSE)
s <- gsub("__LABEL__", label, s, fixed = TRUE)
s <- gsub("__FOLDER__", folder, s, fixed = TRUE)
s <- gsub("__DATASET__", dataset, s, fixed = TRUE)
inference <- c("## Statistical evidence", "",
               "Inference was not complete when this draft was created. Do not invent test results.", "")
status <- file.path("results", folder, "inference-status.txt")
if (file.exists(status) && identical(readLines(status), "COMPLETE")) {
  inference <- c("## Statistical evidence", "",
    "Report PERMANOVA's group R-squared and p-value, the adjusted covariates and the separate dispersion check.",
    "Compare proportion-based and equal-depth sensitivity results; explain any disagreement.", "",
    paste0("![Dispersion check.](../results/", folder, "/figures/dispersion.png)"), "")
  if (file.exists(file.path("results", folder, "tables/differential-abundance.tsv")))
    inference <- c(inference, "## Differential abundance", "",
      "Name the OTU, contrast, natural-log effect, uncertainty, BH-adjusted p-value and sensitivity result.",
      "Explain the selection of displayed features and distinguish OTU associations from genus-level claims.", "",
      paste0("![ANCOM-BC2 volcano.](../results/", folder, "/figures/da-volcano.png)"), "",
      paste0("![ANCOM-BC2 effects and intervals.](../results/", folder, "/figures/da-effects.png)"), "")
}
i <- match("__INFERENCE__", s)
if (!is.na(i)) s <- append(s[-i], inference, after = i - 1)
writeLines(s, path)
message("Created ", path, ". Replace prompts with checked results; existing writing will be kept.")
