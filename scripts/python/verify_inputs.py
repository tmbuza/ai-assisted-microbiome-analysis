"""Validate processed shotgun metagenomics inputs."""
import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def read_table(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, [])
        if (not set(required).issubset(header) or len(set(header)) != len(header)
                or any(not name or name != name.strip() for name in header)):
            raise ValueError(f"{path.name}: invalid header. Required: {', '.join(required)}")
        rows = {}
        for line, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise ValueError(f"{path.name}:{line}: wrong number of fields")
            row = dict(zip(header, values))
            sample = row.pop("sample_id")
            if not sample or sample != sample.strip() or sample in rows:
                raise ValueError(f"{path.name}:{line}: blank, padded or duplicate sample ID {sample!r}")
            rows[sample] = row
    if not rows:
        raise ValueError(f"{path.name}: no sample rows")
    return rows


def numeric_table(path, samples, require_positive_total=False):
    rows = read_table(path, ["sample_id"])
    missing = sorted(samples - rows.keys())
    extra = sorted(rows.keys() - samples)
    if missing or extra:
        raise ValueError(f"{path.name}: missing {len(missing)} {missing[:5]}; extra {len(extra)} {extra[:5]}")
    for sample, row in rows.items():
        if not row:
            raise ValueError(f"{path.name}: at least one feature column is required")
        for feature, value in row.items():
            try:
                number = float(value)
            except ValueError:
                number = math.nan
            if not math.isfinite(number) or number < 0:
                raise ValueError(f"{path.name}: {sample}, {feature}: expected a finite nonnegative number")
            row[feature] = number
        if require_positive_total and (not math.isfinite(sum(row.values())) or not sum(row.values())):
            raise ValueError(f"{path.name}: {sample}: taxon total must be positive and finite")
    return rows


def verify_inputs(folder):
    metadata = read_table(folder / "sample-metadata.tsv", ["sample_id", "group"])
    groups = {sample: row["group"] for sample, row in metadata.items()}
    if any(not group or group != group.strip() for group in groups.values()):
        raise ValueError("sample-metadata.tsv: blank or padded group label")
    counts = Counter(groups.values())
    if len(counts) < 2 or min(counts.values()) < 2:
        raise ValueError("sample-metadata.tsv: require at least two groups with two samples each")
    taxa = numeric_table(folder / "taxon-abundance.tsv", set(groups), True)
    functions = numeric_table(folder / "functional-features.tsv", set(groups))
    return groups, taxa, functions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    args = parser.parse_args()
    try:
        groups, _, _ = verify_inputs(args.inputs)
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(2, f"Input check failed: {error}\n")
    print(f"Verified {len(groups)} matched samples. Groups: {dict(Counter(groups.values()))}")
