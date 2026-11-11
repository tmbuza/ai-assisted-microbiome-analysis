# Mouse59 input preparation

Source: user-supplied Mouse59GutOTUs.csv. Publication, processing history, treatment definitions and sample independence have not been verified.

Extract data/my-study/ into the repository root. Back up existing my-study inputs before replacing them.

- feature-counts.tsv: transposed original counts, 2,027 OTUs by 59 samples. No filtering, normalization or rarefaction.
- sample-metadata.tsv: SampleID renamed sample_id; Class renamed group. Labels retained verbatim: HighFatCa (30), HighFat20 (29). Do not infer treatment meanings from these labels alone.
- taxonomy.tsv: structural placeholder. taxon repeats each OTU ID. No biological taxonomy was supplied or inferred. Composition plots describe OTUs, not named taxa. Replace with verified matching taxonomy when available.

All samples have positive totals, all counts are nonnegative integers, sample and feature IDs are unique and aligned, and all features have positive totals. R workflow execution was not tested here.

Lowest sample totals (count, sample, group):
(99, 'SRR4125132', 'HighFat20')
(182, 'SRR4125107', 'HighFatCa')
(4420, 'SRR4125126', 'HighFat20')
(4550, 'SRR4125145', 'HighFat20')
(8991, 'SRR4125149', 'HighFat20')

Investigate the lowest count sample, sequencing depth differences, cage/animal IDs, repeated sampling and treatment definitions before making biological conclusions.

Run from the repository root:

    Rscript scripts/R/01-check-inputs.R data/my-study results/my-study
    Rscript scripts/R/02-analyse-example.R data/my-study results/my-study
