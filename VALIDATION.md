# Delivery validation — 26 November 2026

## Executed here

- R 4.5.0; vegan 2.6.10; ggplot2 3.5.1; ANCOMBC 2.16.0.
- Both nine-figure exploration workflows completed in R.
- Synthetic community inference and human community plus ANCOM-BC2 inference
  completed. Outputs contain 10 synthetic and 12 human figures including
  inference diagnostics; both run-status files say COMPLETE for each example.
- Shannon formula agrees with vegan for every sample. Human Shannon values
  agree with the earlier independent Python calculation to 1e-12.
- A01 and Healthy1-2027 passed independent sample verification.
- Per-sample full/display proportions sum to one; saved distance matrices
  are symmetric, have zero diagonals and match their sample IDs.
- Freshness guard rejected changed configuration bytes. Fractional-count
  validation rejected altered noninteger counts. Own-study inference is
  blocked by default until the learner confirms an appropriate design.
- Report creation preserved an existing draft. The human report rendered as
  standalone HTML from the isolated reports project.
- The 26-page Quarto book rendered with Quarto 1.11.9. Local guide links and
  figure targets were checked. PCoA, heatmap and DA figures were visually
  inspected. PNGs use an explicit white background.
- The provenance export script reproduced all three cleaned human TSV inputs
  byte for byte in a separate directory.

## Human reference results and limitations

The selected table contains 60 samples and 2,078 OTUs. The adjusted
proportion-based PERMANOVA group term has R-squared 0.0199196 and p=0.165;
the equal-depth sensitivity gives R-squared 0.0198492 and p=0.172.
Unadjusted dispersion p-values are 0.420 and 0.428. These are not proof of
identical community profiles or equal dispersion.

The 10% prevalence filter retained 1,915 OTUs. No OTU met both adjusted-p and
conservative pseudo-count sensitivity criteria. The package's raw coefficient
table is retained alongside the learner-facing summary.

A numerical warning occurs in t-distribution p-value calculation. One hundred
OTUs lack finite group effects or standard errors. All 100 have rank-deficient
positive-count subsets for the full adjustment model. They are flagged as
unavailable, with missing p-values in the learner-facing summary. The raw
ANCOMBC table uses replacement p=q=1 for unavailable tests; that must not be
interpreted as evidence of equality. Diagnostics and the guide explicitly
explain this limitation. The two-group/multi-group-deactivation warning is
expected for this two-group example. No inference-stage warnings are hidden.

Synthetic rarefaction warns that the supplied teaching counts lack count-one
features. Its curves are conditional illustrations with only 15 fictional
features; they are not evidence of adequate biological sampling.

Upstream human count processing was not reconstructed from raw reads. Sample
IDs, disease labels and ages match the original spreadsheet; its participant-ID
field is empty. The omitted two adenoma records have an unresolved omission
reason. Original and cleaned sex/site labels are retained. Biological validity,
full confounder control and clinical screening performance are not established.

## What to check on your computer

Install the documented packages and rerun using your own R environment before
publishing. Versions can affect numerical or plotting details. Read warning
files, the run log and model diagnostics, not only completion markers.
Own-study data are not included and cannot be validated in advance.

The package is an overlay, not a deletion operation. Earlier DietSwap and
Mouse59 source/data folders are left for later cleanup; the new book navigation
uses only synthetic, Zackular and own-study sections. Existing report drafts,
existing .gitignore patterns, deployment files and local styling need to be
preserved when applying the update. Merge the supplied ignore additions.

The generated reports used for validation are intentionally excluded from the
archive. The report template and nonoverwriting creation helper are included.

## Naming update

Synthetic inputs and outputs now consistently use data/synthetic and results/synthetic. Commands, figure resources, report creation and legacy defaults were updated. Input count bytes and analysis settings were unchanged; input path records were migrated without recomputing scientific results.

After renaming, the sample verifier passed for A01 and report creation produced valid synthetic figure paths. Script syntax, input freshness, source links and rendered HTML image paths passed. Input TSV and figure bytes matched the previous archive exactly.
