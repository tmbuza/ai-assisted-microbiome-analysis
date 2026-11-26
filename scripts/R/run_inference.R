# Independent two-group inference, after successful exploration.
source("scripts/R/common.R")
args <- entry_args(); input <- args[1]; output <- args[2]; config <- args[3]
if (dir.exists(output)) writeLines("INCOMPLETE: inference has started but not completed.", file.path(output, "inference-status.txt"))

warning_file <- file.path(output, "inference-warnings.txt")
unlink(warning_file)
withCallingHandlers({
fail(requireNamespace("vegan", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE),
     "Install vegan and ggplot2; see setup.")
library(ggplot2)
cfg <- read_settings(config)
fail(identical(cfg[["design_confirmed"]], "TRUE"),
     "Confirm independent sample units and an appropriate two-group design before inference.")
status <- file.path(output, "inference-status.txt")
fail(file.exists(file.path(output, "exploration-status.txt")) &&
     identical(readLines(file.path(output, "exploration-status.txt")), "COMPLETE"),
     "Complete exploration first.")
e <- readRDS(file.path(output, "exploration.rds"))
verify_fresh(input, config, e$checksums)
writeLines("INCOMPLETE: do not interpret inferential outputs until this run finishes.", status)
# Remove only files owned by this stage, so disabled DA cannot expose old results.
unlink(file.path(output, "ancombc2-fit.rds"))
unlink(file.path(tab <- file.path(output, "tables"), c("permanova.tsv", "dispersion-test.tsv",
  "dispersion-distances.tsv", "da-filter.tsv", "ancombc2-all-coefficients.tsv",
  "differential-abundance.tsv", "da-discoveries.tsv", "da-effect-display.tsv", "feature-model-support.tsv")))
unlink(file.path(output, "figures", c("dispersion.png", "da-volcano.png", "da-effects.png")))

a <- e$a; m <- e$m
fail(all(table(m$group) >= 3), "Inference requires at least three independent samples per group.")
covs <- if (cfg[["covariates"]] == "none") character() else trimws(strsplit(cfg[["covariates"]], ",")[[1]])
fail(all(covs %in% names(m)) && !anyDuplicated(covs) && !"group" %in% covs &&
     all(make.names(covs) == covs), "Covariates must be distinct, valid metadata column names.")
for (v in covs) {
  fail(!anyNA(m[[v]]) && all(nzchar(trimws(m[[v]]))), paste("Missing covariate:", v))
  m[[v]] <- type.convert(m[[v]], as.is = TRUE)
  if (!is.numeric(m[[v]])) m[[v]] <- factor(m[[v]])
  fail(length(unique(m[[v]])) > 1, paste("Covariate has only one value:", v))
}
# Always encode reference/comparison as safe internal names. The output records
# their original clinical labels, so spaces or punctuation cannot reverse direction.
m$group <- factor(ifelse(as.character(m$group) == cfg[["reference"]], "Reference", "Comparison"),
                   levels = c("Reference", "Comparison"))
form <- reformulate(c(covs, "group"))
design <- model.matrix(form, m)
fail(qr(design)$rank == ncol(design), "Model is rank deficient; inspect confounding and covariates.")
fail(nrow(design) - ncol(design) >= 5, "Fewer than five residual degrees of freedom; simplify the design.")
tab <- file.path(output, "tables"); fig <- file.path(output, "figures")
write_table(data.frame(sample_id = rownames(design), design, check.names = FALSE),
            file.path(tab, "model-matrix.tsv"))
write_table(data.frame(sample_id = rownames(m), m, check.names = FALSE)[, !duplicated(c("sample_id", names(m)))],
            file.path(tab, "inference-metadata.tsv"))
seed <- as.numeric(cfg[["seed"]]); permutations <- as.numeric(cfg[["permutations"]])
permanova <- list(); dispersions <- list(); dispersion_rows <- list()
for (method in c("proportions", "equal-depth-sensitivity")) {
  dis <- if (method == "proportions") e$distance else e$depth_distance
  # Matching Lingoes geometry is used for PCoA, PERMANOVA and dispersion.
  # by='margin' tests each term conditional on the others.
  environment(form) <- environment()
  full_form <- as.formula(paste("dis", paste(deparse(form), collapse = " ")))
  set.seed(seed)
  test <- vegan::adonis2(full_form, data = m, permutations = permutations,
                         by = "margin", add = "lingoes")
  permanova[[method]] <- data.frame(method = method, term = rownames(test),
                                     as.data.frame(test), check.names = FALSE)
  bd <- vegan::betadisper(dis, m$group, type = "median", bias.adjust = TRUE, add = "lingoes")
  set.seed(seed)
  bt <- vegan::permutest(bd, permutations = permutations)
  dispersions[[method]] <- data.frame(method = method, term = rownames(bt$tab),
                                      as.data.frame(bt$tab), check.names = FALSE)
  dispersion_rows[[method]] <- data.frame(method = method, sample_id = names(bd$distances),
      group = e$m$group[match(names(bd$distances), e$m$sample_id)], distance_to_median = bd$distances)
}
write_table(do.call(rbind, permanova), file.path(tab, "permanova.tsv"))
write_table(do.call(rbind, dispersions), file.path(tab, "dispersion-test.tsv"))
dp <- do.call(rbind, dispersion_rows)
write_table(dp, file.path(tab, "dispersion-distances.tsv"))
p <- ggplot(dp, aes(group, distance_to_median, colour = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = .12, height = 0, seed = 42)) +
  facet_wrap(~method) + theme_minimal() + guides(colour = "none") +
  labs(x = "Group", y = "Distance to group spatial median",
       subtitle = "Dispersion is a separate check; this plot is not a disease-effect estimate")
ggsave(file.path(fig, "dispersion.png"), p, width = 10, height = 5, dpi = 180, bg = "white")
if (cfg[["run_da"]] != "TRUE") {
  writeLines(c("Community tests complete. Differential abundance disabled in settings.",
    paste("Model:", paste(deparse(form), collapse = " "))), file.path(output, "inference-log.txt"))
  capture.output(sessionInfo(), file = file.path(output, "inference-session-info.txt"))
  writeLines("COMPLETE", status)
  message("Community inference complete; DA disabled for this dataset.")
  quit(save = "no", status = 0)
}
fail(requireNamespace("ANCOMBC", quietly = TRUE), "Install ANCOMBC from Bioconductor; see setup.")
fail(utils::packageVersion("ANCOMBC") >= numeric_version("2.16.0"),
     "This script requires ANCOMBC >= 2.16.0. Use a compatible R/Bioconductor release; see setup.")
# Prevalence filtering applies only to the feature model, not to exploration.
prevalence <- rowMeans(a > 0)
keep <- prevalence >= as.numeric(cfg[["prevalence"]]) & rowSums(a) > 0
write_table(data.frame(feature_id = rownames(a), prevalence = prevalence, retained_for_da = keep),
            file.path(tab, "da-filter.tsv"))
fail(sum(keep) >= 2, "Fewer than two features pass the DA prevalence filter.")
set.seed(seed)
fit <- ANCOMBC::ancombc2(data = a[keep, , drop = FALSE], meta_data = m,
    taxa_are_rows = TRUE, fix_formula = paste(c(covs, "group"), collapse = " + "),
    rand_formula = NULL, p_adj_method = "BH", pseudo_sens = TRUE, conservative = TRUE,
    prv_cut = 0, lib_cut = 0, group = "group", struc_zero = FALSE,
    neg_lb = FALSE, alpha = as.numeric(cfg[["alpha"]]), n_cl = 1,
    global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE, verbose = TRUE)
saveRDS(fit, file.path(output, "ancombc2-fit.rds"))
r <- as.data.frame(fit$res)
write_table(r, file.path(tab, "ancombc2-all-coefficients.tsv"))
suffix <- "groupComparison"
required <- paste0(c("lfc_", "se_", "p_", "q_", "diff_", "passed_ss_"), suffix)
fail(all(c("taxon", required) %in% names(r)),
     "Unexpected ANCOMBC result columns. Inspect package version and the full coefficient table.")
res <- data.frame(feature_id = r$taxon, lfc = r[[paste0("lfc_", suffix)]],
    se = r[[paste0("se_", suffix)]], pvalue = r[[paste0("p_", suffix)]],
    padj = r[[paste0("q_", suffix)]], declared_da = r[[paste0("diff_", suffix)]],
    passed_sensitivity = r[[paste0("passed_ss_", suffix)]])
res$estimate_available <- is.finite(res$lfc) & is.finite(res$se) & res$se > 0
write_table(model_support(a, design, res$feature_id, res$estimate_available),
            file.path(tab, "feature-model-support.tsv"))
# ANCOMBC assigns p=q=1 to non-estimable tests. Preserve its raw table,
# but make unavailable effects explicit in the learner-facing summary.
res$pvalue[!res$estimate_available] <- NA_real_
res$padj[!res$estimate_available] <- NA_real_
res$robust_discovery <- !is.na(res$padj) & res$padj < as.numeric(cfg[["alpha"]]) &
                        !is.na(res$passed_sensitivity) & res$passed_sensitivity == TRUE &
                        !is.na(res$declared_da) & res$declared_da == TRUE
res$ci_low <- res$lfc - 1.96 * res$se
res$ci_high <- res$lfc + 1.96 * res$se
res$fold_change <- exp(res$lfc)
res$comparison <- cfg[["comparison"]]; res$reference <- cfg[["reference"]]
labels <- feature_labels(e$tx); res$label <- unname(labels[res$feature_id])
res <- res[order(is.na(res$padj), res$padj, res$feature_id), ]
write_table(res, file.path(tab, "differential-abundance.tsv"))
write_table(res[res$robust_discovery, ], file.path(tab, "da-discoveries.tsv"))
vol <- res[is.finite(res$lfc) & is.finite(res$padj), ]
vol$minus_log10_padj <- -log10(pmax(vol$padj, .Machine$double.xmin))
p <- ggplot(vol, aes(lfc, minus_log10_padj, colour = robust_discovery)) +
  geom_point(alpha = .7) + geom_hline(yintercept = -log10(as.numeric(cfg[["alpha"]])), linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey60") + theme_minimal() +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#0072B2")) +
  labs(x = paste("Natural-log fold change:", cfg[["comparison"]], "versus", cfg[["reference"]]),
       y = "−log10(BH-adjusted p-value)", colour = "Sensitivity-supported",
       subtitle = paste(sum(res$robust_discovery), "features meet both criteria; no effect-size cutoff"))
ggsave(file.path(fig, "da-volcano.png"), p, width = 9, height = 5.5, dpi = 180, bg = "white")
ef <- head(res[is.finite(res$lfc) & is.finite(res$se) & !is.na(res$padj), ], 20)
ef$label <- factor(ef$label, levels = rev(ef$label))
write_table(ef, file.path(tab, "da-effect-display.tsv"))
p <- ggplot(ef, aes(lfc, label, colour = robust_discovery)) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label)) + geom_point() +
  geom_vline(xintercept = 0, linetype = 2) + theme_minimal(base_size = 11) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#0072B2")) +
  labs(x = "Natural-log fold change with approximate pointwise 95% interval", y = "Feature",
       colour = "Sensitivity-supported", subtitle = "Up to 20 features ordered by adjusted p-value; intervals are not multiplicity-adjusted")
ggsave(file.path(fig, "da-effects.png"), p, width = 12, height = 8, dpi = 180, bg = "white")
writeLines(c(paste("Comparison:", cfg[["comparison"]], "versus", cfg[["reference"]]),
    paste("Model:", paste(deparse(form), collapse = " ")),
    "ANCOM-BC2: raw supplied integer counts; BH correction; pseudo-count sensitivity enabled.",
    "Structural-zero calling disabled; no manual pseudo-count added.",
    paste("Features retained:", sum(keep)),
    paste("Sensitivity-supported discoveries:", sum(res$robust_discovery)),
    paste("Unavailable group effects or standard errors:", sum(!res$estimate_available)),
    "PERMANOVA: marginal terms, Lingoes correction, independent-sample permutations.",
    "Dispersion is unadjusted for covariates; a nonsignificant test does not prove equal dispersion."),
    file.path(output, "inference-log.txt"))
capture.output(sessionInfo(), file = file.path(output, "inference-session-info.txt"))
writeLines("COMPLETE", status)
message("Inference complete. Inspect model, dispersion, sensitivity and DA tables before interpreting.")

}, warning = function(w) {
  cat(conditionMessage(w), "\nCall: ", paste(deparse(conditionCall(w)), collapse = " "), "\n", file = warning_file, append = TRUE)
  if (grepl("converg|singular", conditionMessage(w), ignore.case = TRUE))
    stop("Model warning needs review: ", conditionMessage(w), call. = FALSE)
  invokeRestart("muffleWarning")
})
if (file.exists(warning_file)) message("Review recorded warnings: ", warning_file)
