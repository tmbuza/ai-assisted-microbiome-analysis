# Metagenomics workflow test

Copy the three TSVs into data/metagenomics/my-study/ from the repository root.

- sample-metadata.tsv preserves the CSV sample IDs and maps Disease_State to group. The supplied Shannon index is included for provenance.
- functional-features.tsv preserves all 14 functional columns and their values. Their units and derivation remain unverified.
- taxon-abundance.tsv contains 12 simulated taxa generated with a fixed random seed (42). The CSV contains no taxon abundances. The simulated values were generated without using group, enterotype, Shannon, or functional values.

Use the bundle to test validation, analysis, and plotting. Composition and Shannon results from the workflow describe simulated taxa. Replace the synthetic taxon table with a verified sample-matched profile before biological interpretation. The Shannon values computed by the workflow are distinct from the source Shannon values in metadata.

Run from the repo root with the Python scripts already installed in scripts/python/.

```bash
python3 scripts/python/verify_inputs.py data/metagenomics/my-study
python3 scripts/python/run_workflow.py data/metagenomics/my-study results/metagenomics/my-study
```
