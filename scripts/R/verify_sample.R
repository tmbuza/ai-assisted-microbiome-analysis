source("scripts/R/common.R")
args <- commandArgs(trailingOnly = TRUE)
fail(length(args) == 4, "Usage: Rscript scripts/R/verify_sample.R input output config sample_id")
e <- readRDS(file.path(args[2], "exploration.rds"))
fail(identical(readLines(file.path(args[2], "exploration-status.txt")), "COMPLETE"), "Exploration is incomplete.")
verify_fresh(args[1], args[3], e$checksums)
d <- read_inputs(args[1]); id <- args[4]
fail(id %in% colnames(d$a), "Requested sample is absent.")
v <- d$a[, id]
manual <- c(sum(v), sum(v > 0), shannon(matrix(v, ncol = 1)))
reported <- read.delim(file.path(args[2], "tables/alpha-diversity.tsv"), check.names = FALSE)
reported <- reported[reported$sample_id == id, ]
fail(nrow(reported) == 1, "Expected exactly one saved sample row.")
saved <- as.numeric(reported[1, c("feature_total", "observed_features", "shannon")])
comparison <- data.frame(measure = c("feature_total", "observed_features", "shannon"),
                          manual = manual, saved = saved, agrees = abs(manual - saved) < 1e-6)
print(comparison, row.names = FALSE, digits = 9)
fail(all(comparison$agrees), "Saved values disagree with the current input.")
message("Verified ", id, ". This checks arithmetic and freshness, not biological validity.")
