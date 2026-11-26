"""Rebuild the teaching inputs from the included source files.

Requires: python -m pip install rdata pandas numpy openpyxl
Run: python prepare_inputs.py
This reads the packaged counts; it does not process raw sequencing reads.
"""
from pathlib import Path
from dataclasses import replace
import hashlib
import json
import warnings
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent
warnings.filterwarnings("ignore", message="Missing constructor for R class")
p = rdata.parser.parse_file(ROOT / "zackular2014.rda")
x = rdata.conversion.convert(p)["zackular2014"]
def convert(obj):
    return rdata.conversion.convert(replace(p, object=obj))
slots = {}
node = p.object.value[0].attributes
while node is not None and node.value is not None:
    slots[str(convert(node.tag))] = node.value[0]
    node = node.value[1]
tax_attrs = convert(slots["tax_table"].attributes)
meta_attrs = convert(slots["sam_data"].attributes)
features = list(map(str, x.otu_table.dim_0.values))
samples = list(map(str, x.otu_table.dim_1.values))
counts = pd.DataFrame(x.otu_table.values, index=features, columns=samples)
tax = pd.DataFrame(np.asarray(x.tax_table).reshape(tuple(tax_attrs["dim"]), order="F"),
                   index=tax_attrs["dimnames"][0], columns=tax_attrs["dimnames"][1])
meta = pd.DataFrame(x.sam_data, index=list(map(str, meta_attrs["row.names"])))
original = pd.read_excel(ROOT / "EDRN_MIMARKS.xlsx").set_index("sample_id")
assert counts.index.is_unique and counts.columns.is_unique
assert tax.index.is_unique and meta.index.is_unique
assert set(counts.index) == set(tax.index)
assert set(samples) == set(meta.index)
assert set(samples) <= set(original.index)
a = counts.to_numpy()
assert np.isfinite(a).all() and (a >= 0).all() and (a == np.floor(a)).all()
assert (a.sum(axis=0) > 0).all()
expected = meta.DiseaseState.map({"H": "normal", "CRC": "carcinoma", "nonCRC": "adenoma"})
assert (expected == original.loc[meta.index, "disease_stat"]).all()
original_age = original.loc[meta.index, "age"].astype(str).str.strip().astype(float)
assert np.allclose(meta.age.astype(float), original_age)
assert original.loc[meta.index, "target_subfragment"].eq("V4").all()
selected = sorted(meta.index[meta.DiseaseState.isin(["H", "CRC"])])
selected_counts = counts.loc[:, selected]
out_meta = pd.DataFrame({"sample_id": selected,
    "group": meta.loc[selected, "DiseaseState"].map({"H": "Healthy", "CRC": "Carcinoma"}).values,
    "original_disease_state": meta.loc[selected, "DiseaseState"].values,
    "age": meta.loc[selected, "age"].values,
    "sex": original.loc[selected, "sex"].values,
    "collection_location": original.loc[selected, "geo_loc_name"].values,
    "fobt_result": meta.loc[selected, "FOBT.result"].values})
assert out_meta.group.value_counts().to_dict() == {"Carcinoma": 30, "Healthy": 30}
assert len(selected) == 60
out_tax = pd.DataFrame({"feature_id": features})
out_tax["taxon"] = ["; ".join(map(str, tax.loc[f].values)) for f in features]
out = ROOT.parents[2] / "data/zackular"
out.mkdir(parents=True, exist_ok=True)
selected_counts.rename_axis("feature_id").to_csv(out / "feature-counts.tsv", sep="\t")
out_meta["sex_original"] = out_meta.sex
out_meta["collection_location_original"] = out_meta.collection_location
out_meta.sex = out_meta.sex.str.strip().map({"f": "Female", "m": "Male"})
out_meta.collection_location = out_meta.collection_location.str.strip().replace({"Ann Arbor,. MI": "Ann Arbor, MI"})
assert not out_meta.sex.isna().any()
out_meta.to_csv(out / "sample-metadata.tsv", sep="\t", index=False)
out_tax.to_csv(out / "taxonomy.tsv", sep="\t", index=False)
assert np.array_equal(pd.read_csv(out / "feature-counts.tsv", sep="\t", index_col=0).values,
                      selected_counts.values)
b = selected_counts.values.astype(float)
proportions = b / b.sum(axis=0)
terms = np.zeros_like(proportions)
positive = proportions > 0
terms[positive] = proportions[positive] * np.log(proportions[positive])
alpha = pd.DataFrame({"sample_id": selected, "feature_total": b.sum(axis=0).astype(int),
    "observed_features": (b > 0).sum(axis=0), "shannon": -terms.sum(axis=0)})
alpha = alpha.merge(out_meta[["sample_id", "group"]], on="sample_id", validate="one_to_one")
review = ROOT / "export-audit"
review.mkdir(exist_ok=True)
alpha.to_csv(review / "alpha-diversity-python-check.tsv", sep="\t", index=False)
summary = alpha.groupby("group").agg(samples=("sample_id", "size"),
    min_total=("feature_total", "min"), max_total=("feature_total", "max"),
    median_shannon=("shannon", "median"), min_shannon=("shannon", "min"),
    max_shannon=("shannon", "max"))
summary.to_csv(review / "group-summary-python-check.tsv", sep="\t")
audit = {"source_samples": len(samples), "source_features": len(features),
    "source_groups": meta.DiseaseState.value_counts().to_dict(),
    "absent_from_package": sorted(set(original.index) - set(samples)),
    "selected_samples": len(selected), "retained_features": len(features),
    "all_zero_features_in_selection": int((b.sum(axis=1) == 0).sum()),
    "original_subject_id_nonmissing": int(original.host_subject_id.notna().sum()),
    "unique_numeric_tokens_in_selected_sample_names": len({s.split('-')[1] for s in selected}),
    "source_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                      for f in [ROOT / "zackular2014.rda", ROOT / "EDRN_MIMARKS.xlsx"]}}
(review / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
print(summary.to_string())
print(json.dumps(audit, indent=2))
