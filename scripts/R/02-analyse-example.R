# 1. Validate inputs and load packages -------------------------------------

# Run from the repository root. Both directory arguments are optional.
# Usage: Rscript scripts/R/02-analyse-example.R INPUT_DIRECTORY OUTPUT_DIRECTORY
# The validation script reads inputs, defines x, m, t, and a, and creates tables/.
source("scripts/R/01-check-inputs.R")
suppressPackageStartupMessages(library(vegan))

# 2. Prepare the figure directory ------------------------------------------

figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Calculate within-sample diversity -------------------------------------

# vegan expects samples in rows, so transpose the count matrix.
# Shannon diversity uses all features and the natural logarithm.
alpha <- tibble(
  sample_id = names(x)[-1],
  feature_total = colSums(a),
  observed_features = colSums(a > 0),
  shannon = as.numeric(diversity(t(a), index = "shannon"))
) |>
  left_join(m, by = "sample_id")

# 4. Calculate relative abundance ------------------------------------------

# Divide each count by its own sample total, then attach taxonomy and metadata.
relative <- x |>
  pivot_longer(-feature_id, names_to = "sample_id", values_to = "count") |>
  group_by(sample_id) |>
  mutate(proportion = count / sum(count)) |>
  ungroup() |>
  left_join(t, by = "feature_id") |>
  left_join(m, by = "sample_id")

# 5. Prepare the composition display ---------------------------------------

# Rank features by mean sample proportion, with feature ID breaking ties.
# Retain eight features and combine the remainder as Other for display only.
top_ids <- relative |>
  group_by(feature_id) |>
  summarise(mean_proportion = mean(proportion), .groups = "drop") |>
  arrange(desc(mean_proportion), feature_id) |>
  slice_head(n = 8) |>
  pull(feature_id)

# Include feature IDs in labels to distinguish features sharing a taxon name.
labels <- paste(t$taxon, paste0("[", t$feature_id, "]"))
names(labels) <- t$feature_id

plot_data <- relative |>
  mutate(
    display_feature = if_else(
      feature_id %in% top_ids, labels[feature_id], "Other"
    )
  ) |>
  group_by(sample_id, group, display_feature) |>
  summarise(proportion = sum(proportion), .groups = "drop")
plot_data$display_feature <- factor(
  plot_data$display_feature,
  levels = c(
    unname(labels[top_ids]),
    if (any(plot_data$display_feature == "Other")) "Other"
  )
)

# Keep samples ordered consistently by group, then sample ID.
sample_order <- m |>
  arrange(group, sample_id) |>
  pull(sample_id)
plot_data$sample_id <- factor(plot_data$sample_id, levels = sample_order)

# 6. Summarise groups and save tables --------------------------------------

# Report sample counts, Shannon ranges, and count-total ranges for each group.
summary <- alpha |>
  group_by(group) |>
  summarise(
    samples = n(),
    median_shannon = median(shannon),
    min_shannon = min(shannon),
    max_shannon = max(shannon),
    min_total = min(feature_total),
    max_total = max(feature_total),
    .groups = "drop"
  )

write_tsv(alpha, file.path(output_dir, "tables", "alpha-diversity.tsv"))
write_tsv(relative, file.path(output_dir, "tables", "relative-abundance.tsv"))
write_tsv(plot_data, file.path(output_dir, "tables", "composition-display.tsv"))
write_tsv(summary, file.path(output_dir, "tables", "group-summary.tsv"))

# 7. Plot community composition --------------------------------------------

# Each bar shows one sample; displayed proportions sum to one.
p1 <- ggplot(plot_data, aes(sample_id, proportion, fill = display_feature)) +
  geom_col(width = .85) +
  facet_wrap(~group, scales = "free_x") +
  scale_y_continuous(labels = function(z) paste0(round(z * 100), "%")) +
  labs(
    x = "Sample",
    y = "Feature proportion",
    fill = "Feature",
    subtitle = "Top eight features by mean sample proportion; remainder grouped as Other"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# 8. Plot Shannon diversity by group ---------------------------------------

# Show every sample. A fixed jitter seed keeps point positions reproducible.
# Hide boxplot outlier symbols because all samples are plotted separately.
p2 <- ggplot(alpha, aes(group, shannon, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = .5) +
  geom_point(
    position = position_jitter(width = .12, height = 0, seed = 42),
    size = 2.5
  ) +
  labs(
    x = "Group",
    y = "Shannon diversity (natural log)",
    subtitle = "Each point is a sample; all features included"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# 9. Inspect count totals alongside diversity ------------------------------

# Use this plot to inspect the relationship between count total and diversity.
p3 <- ggplot(alpha, aes(feature_total, shannon, colour = group)) +
  geom_point(size = 2.5) +
  labs(
    x = "Feature-table count total",
    y = "Shannon diversity (natural log)",
    colour = "Group"
  ) +
  theme_minimal(base_size = 12)

# 10. Save figures and session information ---------------------------------

# Allow wider figures for more samples or groups, within the existing limits.
# Review readability for very large studies.
ggsave(
  file.path(figure_dir, "composition.png"),
  p1,
  width = max(10, min(24, nrow(m) * .35)),
  height = 5.5,
  dpi = 180
)
ggsave(
  file.path(figure_dir, "shannon-diversity.png"),
  p2,
  width = max(7, n_distinct(m$group) * 1.5),
  height = 5,
  dpi = 180
)
ggsave(
  file.path(figure_dir, "depth-diversity.png"),
  p3,
  width = 7,
  height = 5,
  dpi = 180
)

# Record package and R versions alongside the outputs.
capture.output(sessionInfo(), file = file.path(output_dir, "session-info.txt"))
message("Outputs saved to ", output_dir)
