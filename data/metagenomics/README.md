# CDI Shotgun Metagenomics Prototype

**Case question**  
How do gut microbial composition, diversity, and selected functional features differ across blood pressure groups?

## Run the synthetic example

Use Python 3.10 or later. Open a terminal in this folder.

```bash
python3 verify_inputs.py example
python3 run_workflow.py example results/example
```

The first command is an optional preflight. The second always runs the same checks before analysis. A failed check exits with status 2 and leaves the requested result directory unwritten. Choose a new output directory for each run.

## Prepare your study inputs

Use UTF-8 tab-separated files with a header and one row per sample. Use exactly the same sample IDs in all three files. Row order can differ.

| File | Columns |
|---|---|
| sample-metadata.tsv | sample_id, group, and optional metadata |
| taxon-abundance.tsv | sample_id followed by taxa at a single rank, with counts or linear abundances from one profiling method |
| functional-features.tsv | sample_id followed by selected functional features with documented comparable units within each feature |

The verifier checks filenames, headers, row widths, unique sample IDs, exact sample matching, group sizes, finite nonnegative values, and positive taxon totals. The design assumes independent samples and at least two samples in each of at least two groups. Review sample independence, profiling method, taxonomic rank, feature units, and normalization before interpreting the output. Missing values require an explicit upstream decision. Zero represents a supplied numeric zero.

## Apply it to the hypertension cohort

Start from BioProject PRJEB13870 and the corresponding SPIRE taxonomic profile. Use a verified run-to-sample mapping before converting profiles to the required format. Use SampleID and Disease_State from the supplied CSV to construct metadata, then move the 14 functional columns into the functional table. Keep Enterotype as optional metadata and retain the supplied Shannon index separately for provenance. Document the original and retained counts, all exclusions, and the functional feature definitions and units before running the real case. The public-profile join and functional units are pending verification.

## Check the outputs

- verified-groups.tsv records the matched sample counts.
- composition-by-group.tsv gives the median within-sample relative abundance for each taxon by group. Taxon medians need not sum to one.
- shannon-by-sample.tsv and shannon-by-group.tsv give natural-log Shannon diversity and its group median. S01 in the synthetic example should be approximately 0.693147.
- functional-features-by-group.tsv reports group medians in the supplied feature units.

The bundled input is synthetic. The workflow starts from processed shotgun profiles and produces descriptive summaries. Report observed group differences as associations and describe the study design, potential confounders, and feature definitions alongside them.
