# Expanded community comparison

Updated dataset: 24 samples, 15 fictional features, varying sample totals, within-group variation, and overlapping diversity values.

Scripts accept custom input/output directories and group labels. Composition shows the top eight features plus Other; diversity uses all features. A third figure displays count totals versus diversity.

The blank template is restored in templates/learner-report.qmd. Root learner-report.qmd and the book example report are filled examples. CDI message and note blocks emphasize findings and limits. The supplied renv.lock is unchanged.

Extract this folder and copy its contents into your existing repository, allowing replacement of matching files. Preserve your .Rprofile and renv directory. Run the analysis before rendering. Existing docs may contain obsolete pages from earlier layouts; check the newly rendered navigation.

Verification: table matching, integer counts, positive totals, proportion sums, Shannon values, group overlap, YAML parsing, and document links checked independently. R and Quarto execution were unavailable. Included preview figures were generated with Python; R regenerates them with its plotting style.
