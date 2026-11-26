# AI-Assisted Microbiome Analysis with R

Three examples follow the same seven steps: synthetic practice, the Zackular
human stool study, and the learner's own study. Exploration and inference are
separate commands with recorded settings and input freshness checks.

Start at `00-setup.qmd`. Run exploration and inference as documented on each
example's Run page. New outputs replace named files; retain copies when
comparing alternative settings. Never interpret an incomplete run.

| Example | Inputs | Outputs | Settings |
|---|---|---|---|
| Synthetic | data/synthetic | results/synthetic | config/synthetic.tsv |
| Human stool | data/zackular | results/zackular | config/zackular.tsv |
| Your study | data/my-study | results/my-study | config/my-study.tsv |

Exploration creates nine plots. Inference adds dispersion and, for the human
example, ANCOM-BC2 differential abundance plots. Sample verification is a
separate command. See `01-run-the-workflow.qmd` for the output map and
`99-appendix.qmd` for assumptions, depth sensitivity, natural-log effects,
multiple testing, provenance and references.

Reports are created by `scripts/R/create_report.R` in ignored `reports/`.
A draft is never overwritten. Render a report from inside that folder.

This update overlays the existing project. Back it up first, then copy files
at their matching relative paths. Navigation now lists only the three retained
examples. Existing DietSwap/Mouse59 source or data folders are not deleted;
remove them separately later if desired. Older scripts remain for legacy
examples. Do not retain old generated `docs/` pages as evidence of the new book.
Render the updated book with `quarto render` after completing the documented
runs. Keep existing project-specific deployment settings and styling if needed.
The update includes CDI styling adapted from the RNA-seq guide, with Quarto's cosmo theme. The stylesheet has its own name so it does not overwrite an existing cdi.css.

Do not overwrite your existing .gitignore: add the patterns in
`gitignore-additions.txt` to it. This keeps reports and learner inputs/results
out of version control. Keep supplied example data and provenance available.

See `VALIDATION.md` for checks performed in the delivery environment and any
remaining limitations. No private learner dataset is included.

Warnings are saved in `exploration-warnings.txt` or `inference-warnings.txt` when present. Read them before accepting results; model convergence or singularity warnings stop inference. Synthetic rarefaction warnings about the lack of count-one features reflect the supplied teaching table and limit biological interpretation.
