# Zackular source records

See the guide appendix for the human-study citation, dataset selection,
metadata harmonisation and remaining provenance limits.

The included zackular2014.rda was downloaded from microsud/microbiomeutilities,
repository HEAD 046a9f928f35fc2a02dd6cfb9175a0c7ddf65c4b. Its SHA-256 is
48d3b9e2ff118f6f2353484b599fe48eb9f7a9aae09709d7d13d73018ad4228c.

EDRN_MIMARKS.xlsx is the original 90-row study spreadsheet from
https://mothur.s3.us-east-2.amazonaws.com/data/MicrobiomeBiomarkerCRC/EDRN_MIMARKS.xlsx
SHA-256: 25312980e24c128f6748542cfc4dc84a5d5963d1c77cc9a22a084465f68b6d10.

prepare_inputs.py reconstructs the three data/zackular TSVs from these files,
including the documented sex and location label cleaning. It requires Python
with rdata, pandas, numpy and openpyxl; it is a maintainer preparation tool,
not part of the learner's R workflow. Do not rerun it to repair an unexplained
count discrepancy. A reproduction should agree before replacing existing data.

All 88 packaged sample IDs match original records, and ages and disease
labels agree. The 60 Healthy/Carcinoma records are all retained. Two absent
adenoma records are Adenoma1-2065 and Adenoma12-2799; the omission reason
remains unresolved. Original host_subject_id values are empty. No independent
participant crosswalk or full upstream processing audit is claimed.
