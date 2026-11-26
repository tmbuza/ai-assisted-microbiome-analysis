# Descriptive exploration: no p-values are computed here.
source("scripts/R/common.R")
args <- entry_args(); input <- args[1]; output <- args[2]; config <- args[3]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
writeLines("INCOMPLETE: exploration has started but not completed.", file.path(output, "exploration-status.txt"))
unlink(file.path(output, "inference-status.txt"))

warning_file <- file.path(output, "exploration-warnings.txt")
unlink(warning_file)
withCallingHandlers({
fail(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("vegan", quietly = TRUE),
     "Install ggplot2 and vegan using the setup page.")
library(ggplot2)
cfg <- read_settings(config); d <- read_inputs(input)
a <- d$a; m <- d$m; tx <- d$tx
fail(setequal(unique(m$group), c(cfg[["reference"]], cfg[["comparison"]])),
     "This two-group workflow requires the reference and comparison labels in settings.")
m$group <- factor(m$group, levels = c(cfg[["reference"]], cfg[["comparison"]]))
tab <- file.path(output, "tables"); fig <- file.path(output, "figures")
dir.create(tab, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
status <- file.path(output, "exploration-status.txt")
writeLines("INCOMPLETE: do not interpret this folder until exploration finishes.", status)
unlink(file.path(output, "inference-status.txt"))
checksums <- sha_inputs(input, config)
write_table(checksums, file.path(tab, "input-checksums.tsv"))
write_table(data.frame(parameter = names(cfg), value = unname(cfg)), file.path(tab, "parameters.tsv"))
write_table(data.frame(check = c("samples", "features", "groups", "validation"),
                       value = c(ncol(a), nrow(a), nlevels(m$group), "PASS")),
            file.path(tab, "input-checks.tsv"))
total <- colSums(a); prop <- sweep(a, 2, total, "/")
alpha <- data.frame(sample_id = colnames(a), group = m$group, feature_total = total,
                    observed_features = colSums(a > 0), shannon = shannon(a))
fail(isTRUE(all.equal(unname(alpha$shannon),
    unname(vegan::diversity(t(a), index = "shannon")), tolerance = 1e-10)),
    "Shannon formula and vegan disagree.")
depth <- as.numeric(cfg[["depth"]]); if (depth == 0) depth <- min(total)
fail(depth <= min(total) && depth > 0,
     "Configured depth exceeds a sample total. Review exclusions or choose a lower depth.")
alpha$expected_richness_at_depth <- as.numeric(vegan::rarefy(t(a), sample = depth))
alpha$standard_depth <- depth
write_table(alpha, file.path(tab, "alpha-diversity.tsv"))
summary <- do.call(rbind, lapply(levels(m$group), function(g) {
  z <- alpha[alpha$group == g, ]
  data.frame(group = g, samples = nrow(z), median_shannon = median(z$shannon),
             min_shannon = min(z$shannon), max_shannon = max(z$shannon),
             min_total = min(z$feature_total), max_total = max(z$feature_total),
             median_observed = median(z$observed_features),
             median_expected_richness = median(z$expected_richness_at_depth))
}))
write_table(summary, file.path(tab, "group-summary.tsv"))
labels <- feature_labels(tx)
top <- head(rownames(a)[order(-rowMeans(prop), rownames(a))], as.numeric(cfg[["top_n"]]))
relative <- data.frame(feature_id = rep(rownames(a), ncol(a)),
                       sample_id = rep(colnames(a), each = nrow(a)),
                       count = as.vector(a), proportion = as.vector(prop))
relative$group <- m$group[match(relative$sample_id, m$sample_id)]
relative$taxon <- tx$taxon[match(relative$feature_id, tx$feature_id)]
write_table(relative, file.path(tab, "relative-abundance.tsv"))
display <- relative
display$display_feature <- ifelse(display$feature_id %in% top,
                                  labels[display$feature_id], "Other")
display <- aggregate(proportion ~ sample_id + group + display_feature, display, sum)
display$display_feature <- factor(display$display_feature,
    levels = c(unname(labels[top]), if (any(display$display_feature == "Other")) "Other"))
sample_order <- m$sample_id[order(m$group, m$sample_id)]
display$sample_id <- factor(display$sample_id, levels = sample_order)
write_table(display, file.path(tab, "composition-display.tsv"))
pal <- setNames(grDevices::hcl.colors(length(top), "Dark 3"), labels[top])
pal <- c(pal, Other = "#D9D9D9")
group_pal <- setNames(c("#009E73", "#D55E00"), levels(m$group))
theme_set(theme_minimal(base_size = 12))
save_plot <- function(p, name, width = 8, height = 5.5) {
  ggsave(file.path(fig, paste0(name, ".png")), p, width = width, height = height, dpi = 180, bg = "white")
}
save_plot(ggplot(display, aes(sample_id, proportion, fill = display_feature)) +
    geom_col(width = .85) + facet_wrap(~group, scales = "free_x") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = function(z) paste0(round(z * 100), "%")) +
    labs(x = "Sample", y = "Feature proportion", fill = "Feature",
         subtitle = paste("Top", length(top), "features by mean sample proportion; remainder is Other")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)),
    "composition", max(12, min(24, ncol(a) * .35)), 7)
save_plot(ggplot(alpha, aes(group, shannon, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = .5) +
    geom_point(position = position_jitter(width = .12, height = 0, seed = 42)) +
    scale_colour_manual(values = group_pal) + guides(colour = "none") +
    labs(x = "Group", y = "Shannon diversity (natural log)",
         subtitle = "Each point is a sample; all supplied features included"), "shannon-diversity")
rich <- rbind(data.frame(sample_id = alpha$sample_id, group = alpha$group,
                        measure = "Observed in supplied table", richness = alpha$observed_features),
              data.frame(sample_id = alpha$sample_id, group = alpha$group,
                         measure = paste("Expected at", depth, "reads"),
                         richness = alpha$expected_richness_at_depth))
write_table(rich, file.path(tab, "richness-display.tsv"))
save_plot(ggplot(rich, aes(group, richness, colour = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = .12, height = 0, seed = 42)) +
    facet_wrap(~measure) + scale_colour_manual(values = group_pal) + guides(colour = "none") +
    labs(x = "Group", y = "Feature richness"), "richness", 10)
save_plot(ggplot(alpha, aes(feature_total, shannon, colour = group)) + geom_point(size = 2.5) +
    scale_colour_manual(values = group_pal) +
    labs(x = "Feature-table read total", y = "Shannon diversity (natural log)", colour = "Group"),
    "depth-diversity")
# Exact expected richness, conditional on each supplied count vector.
curve <- do.call(rbind, lapply(seq_len(ncol(a)), function(i) {
  sizes <- unique(round(seq(1, total[i], length.out = 30)))
  data.frame(sample_id = colnames(a)[i], group = m$group[i], reads = sizes,
      expected_richness = vapply(sizes, function(n) as.numeric(vegan::rarefy(a[, i], n)), numeric(1)))
}))
write_table(curve, file.path(tab, "rarefaction.tsv"))
save_plot(ggplot(curve, aes(reads, expected_richness, group = sample_id, colour = group)) +
    geom_line(alpha = .45) + geom_vline(xintercept = depth, linetype = 2) +
    scale_colour_manual(values = group_pal) +
    labs(x = "Subsample read total", y = "Expected detected features", colour = "Group",
         subtitle = "Conditional on the supplied processed table; dashed line is comparison depth"),
    "rarefaction")
distance <- vegan::vegdist(t(prop), method = "bray")
fail(any(as.vector(distance) > 0), "Samples have identical compositions; ordination is undefined.")
ordination <- function(dis, name, subtitle) {
  ord <- vegan::wcmdscale(dis, k = 2, eig = TRUE, add = "lingoes")
  fail(ncol(ord$points) >= 2, "Fewer than two positive ordination axes.")
  percent <- 100 * ord$eig[1:2] / sum(ord$eig[ord$eig > 0])
  scores <- data.frame(sample_id = rownames(ord$points), PCoA1 = ord$points[, 1],
                        PCoA2 = ord$points[, 2])
  scores$group <- m$group[match(scores$sample_id, m$sample_id)]
  write_table(scores, file.path(tab, paste0(name, "-scores.tsv")))
  write_table(data.frame(axis = seq_along(ord$eig), eigenvalue = ord$eig,
                         percent = 100 * ord$eig / sum(ord$eig[ord$eig > 0])),
              file.path(tab, paste0(name, "-variance.tsv")))
  save_plot(ggplot(scores, aes(PCoA1, PCoA2, colour = group)) + geom_point(size = 2.5) +
      scale_colour_manual(values = group_pal) + coord_equal() +
      labs(x = sprintf("PCoA1 (%.1f%%)", percent[1]), y = sprintf("PCoA2 (%.1f%%)", percent[2]),
           colour = "Group", subtitle = subtitle), name)
  ord
}
ord <- ordination(distance, "pcoa", "Bray–Curtis on proportions; Lingoes-corrected PCoA")
# Sensitivity analysis only: never overwrite raw counts used for ANCOM-BC2.
set.seed(as.numeric(cfg[["seed"]]))
repeats <- as.numeric(cfg[["depth_repeats"]]); accum <- matrix(0, ncol(a), ncol(a))
rare_alpha <- matrix(NA_real_, ncol(a), repeats)
for (i in seq_len(repeats)) {
  rare <- vegan::rrarefy(t(a), sample = depth)
  accum <- accum + as.matrix(vegan::vegdist(rare, "bray"))
  rare_alpha[, i] <- vegan::diversity(rare, "shannon")
}
dimnames(accum) <- list(colnames(a), colnames(a))
depth_distance <- as.dist(accum / repeats)
invisible(ordination(depth_distance, "pcoa-depth-sensitivity",
           paste("Mean Bray–Curtis across", repeats, "subsamples at", depth, "reads; Lingoes PCoA")))
write_table(data.frame(sample_id = colnames(a), group = m$group,
    depth = depth, repeats = repeats, median_shannon = apply(rare_alpha, 1, median),
    min_shannon = apply(rare_alpha, 1, min), max_shannon = apply(rare_alpha, 1, max)),
    file.path(tab, "depth-sensitivity-alpha.tsv"))
for (nm in c("distance", "depth_distance")) {
  dm <- as.matrix(get(nm))
  write_table(data.frame(sample_id = rownames(dm), dm, check.names = FALSE),
              file.path(tab, paste0(nm, ".tsv")))
}
heat_ids <- head(rownames(a)[order(-rowMeans(prop), rownames(a))], min(25, nrow(a)))
heat <- relative[relative$feature_id %in% heat_ids, ]
heat$sample_id <- factor(heat$sample_id, levels = sample_order)
heat$label <- factor(labels[heat$feature_id], levels = rev(labels[heat_ids]))
heat$log10_proportion_plus_offset <- log10(heat$proportion + 1e-6)
write_table(heat, file.path(tab, "heatmap-display.tsv"))
save_plot(ggplot(heat, aes(sample_id, label, fill = log10_proportion_plus_offset)) +
    geom_tile() + facet_wrap(~group, scales = "free_x") +
    scale_fill_viridis_c(name = "log10(p + 1e-6)") +
    labs(x = "Sample", y = "Feature", subtitle = "Top 25 features by mean proportion; no row scaling") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 8)),
    "abundance-heatmap", 13, 8)
prev <- do.call(rbind, lapply(levels(m$group), function(g) {
  j <- which(m$group == g)
  data.frame(feature_id = rownames(a), group = g, prevalence = rowMeans(a[, j, drop = FALSE] > 0),
             mean_proportion = rowMeans(prop[, j, drop = FALSE]))
}))
write_table(prev, file.path(tab, "prevalence.tsv"))
pv <- prev[prev$feature_id %in% heat_ids, ]; pv$label <- factor(labels[pv$feature_id], levels = rev(labels[heat_ids]))
save_plot(ggplot(pv, aes(prevalence, label, colour = group, size = mean_proportion)) +
    geom_point(position = position_dodge(width = .4)) + scale_colour_manual(values = group_pal) +
    scale_size_area(max_size = 7) +
    labs(x = "Fraction of samples with count > 0", y = "Feature", colour = "Group", size = "Mean proportion",
         subtitle = "Same top 25 features as the heatmap; detection depends on the supplied table"),
    "prevalence", 11, 8)
saveRDS(list(a = a, m = m, tx = tx, config = cfg, checksums = checksums,
              distance = distance, depth_distance = depth_distance, depth = depth),
        file.path(output, "exploration.rds"))
capture.output(sessionInfo(), file = file.path(output, "session-info.txt"))
writeLines("COMPLETE", status)
message("Exploration complete. Tables and nine figures saved in ", output)

}, warning = function(w) {
  cat(conditionMessage(w), "\n", file = warning_file, append = TRUE)
  if (grepl("converg|singular", conditionMessage(w), ignore.case = TRUE))
    stop("Model warning needs review: ", conditionMessage(w), call. = FALSE)
  invokeRestart("muffleWarning")
})
if (file.exists(warning_file)) message("Review recorded warnings: ", warning_file)
