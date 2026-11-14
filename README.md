# AI-Assisted Microbiome Analysis

Follow a consistent descriptive workflow from three tidy TSV inputs to composition and diversity figures, checked tables, and a short report.

## Learning path

1. Set up the R environment: `00-setup.qmd`.
2. Practise with synthetic data: `examples/synthetic/`.
3. Follow the selected DietSwap samples: `examples/dietswap/`.
4. Apply the workflow to Mouse59 OTUs: `examples/mouse59/`.
5. Prepare your own analysis: `examples/my-study/`.

Every dataset follows the same seven pages, from question to report. Your own study uses the same steps after you prepare its inputs.

Each dataset contains its exact analysis and optional verification commands. Its final page includes a report outline with matching figure paths and a render command.

## Data and outputs

| Dataset | Inputs | Outputs |
|---|---|---|
| Synthetic | `data/example/` | `results/example/` |
| DietSwap | `data/dietswap/` | `results/dietswap/` |
| Mouse59 | `data/mouse59/` | `results/mouse59/` |
| Learner's study | `data/my-study/` | `results/my-study/` |

Each input folder needs `feature-counts.tsv`, `sample-metadata.tsv`, and `taxonomy.tsv`. Keep learner folders empty until supplying your own inputs; create them when needed because Git does not retain empty directories.

Mouse59 contains 59 samples and 2,027 OTUs. Its taxonomy file uses OTU IDs as display labels, and its original Class labels are preserved as groups. Counts were transposed from the supplied CSV without filtering or normalization. It is a workflow demonstration with limited study context.

## Scripts

The main entry point is `scripts/R/02-analyse-example.R`. It calls `01-check-inputs.R` and generates figures and tables. `05-verify-sample.R` is an optional learning exercise. Keep `renv.lock` for package restoration.

This guide starts from prepared feature tables. Raw sequencing processing and statistical inference require additional workflows.

## Render the guide

Generate results for the three supplied datasets using their Run the analysis pages, then run `quarto render` from the repository root. The book uses existing outputs; rendering does not run the analysis. Generated pages are written to `docs/`.
