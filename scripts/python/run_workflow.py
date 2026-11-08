"""Run descriptive comparisons after validating all three input files."""
import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median
from tempfile import TemporaryDirectory


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


def write_figures(folder, groups, proportions, diversity, functions):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ValueError("Install Matplotlib with python3 -m pip install 'matplotlib>=3.7'") from error

    folder.mkdir()
    by_group = defaultdict(list)
    for sample, group in groups.items():
        by_group[group].append(sample)
    samples = [sample for members in by_group.values() for sample in members]
    positions = list(range(len(samples)))
    colors = plt.get_cmap("tab20").colors

    taxa = sorted(proportions[samples[0]],
                  key=lambda taxon: (-mean(proportions[s][taxon] for s in samples), taxon))
    shown = taxa[:10]
    series = [(taxon, [proportions[s][taxon] for s in samples]) for taxon in shown]
    if taxa[10:]:
        series.append(("Remaining taxa (pooled)",
                       [sum(proportions[s][t] for t in taxa[10:]) for s in samples]))

    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    bottom = [0.0] * len(samples)
    for index, (taxon, values) in enumerate(series):
        ax.bar(positions, values, bottom=bottom, width=0.9, color=colors[index], label=taxon)
        bottom = [total + value for total, value in zip(bottom, values)]
    offset = 0
    centers = []
    for members in by_group.values():
        if offset:
            ax.axvline(offset - 0.5, color="0.35", linewidth=0.8)
        centers.append(offset + (len(members) - 1) / 2)
        offset += len(members)
    ax.set_xticks(centers, list(by_group))
    ax.set(ylabel="Relative abundance", xlabel="Group (one bar per sample)",
           title="Taxonomic composition", ylim=(0, 1))
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    fig.savefig(folder / "composition.png", dpi=180)
    plt.close(fig)

    def plot_points(ax, profiles, feature):
        for index, members in enumerate(by_group.values()):
            values = [profiles[s][feature] for s in members]
            offsets = [0.36 * (i / (len(values) - 1) - 0.5) for i in range(len(values))]
            ax.scatter([index + offset for offset in offsets], values, s=24, alpha=0.75)
            ax.hlines(median(values), index - 0.24, index + 0.24, color="black", linewidth=1.5)
        ax.set_xticks(range(len(by_group)), list(by_group), rotation=20, ha="right")
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)

    fig, ax = plt.subplots(figsize=(6, 4.5), layout="constrained")
    plot_points(ax, diversity, "shannon")
    ax.set(title="Shannon diversity", ylabel="Shannon index (natural log)")
    fig.savefig(folder / "shannon-diversity.png", dpi=180)
    plt.close(fig)

    features = list(functions[samples[0]])
    columns = min(3, len(features))
    rows = (len(features) + columns - 1) // columns
    fig, axes = plt.subplots(rows, columns, squeeze=False, layout="constrained",
                             figsize=(4.2 * columns, 3.3 * rows))
    for ax, feature in zip(axes.flat, features):
        plot_points(ax, functions, feature)
        ax.set(title=feature.replace("_", " "), ylabel="Value (input units)")
    for ax in list(axes.flat)[len(features):]:
        ax.set_visible(False)
    fig.suptitle("Selected functional features")
    fig.savefig(folder / "functional-features.png", dpi=180)
    plt.close(fig)


def write_table(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def summarize(groups, profiles):
    values = defaultdict(list)
    for sample, features in profiles.items():
        for feature, value in features.items():
            values[groups[sample], feature].append(value)
    for (group, feature), observations in sorted(values.items()):
        yield group, feature, len(observations), median(observations)


def run(inputs, output):
    groups, taxa, functions = verify_inputs(inputs)
    if output.exists() or output.is_symlink():
        raise ValueError(f"Choose a new output directory: {output}")
    proportions = {}
    diversity = {}
    for sample, row in taxa.items():
        total = sum(row.values())
        profile = {taxon: value / total for taxon, value in row.items()}
        proportions[sample] = profile
        shannon = -sum(p * math.log(p) for p in profile.values() if p > 0)
        diversity[sample] = {"shannon": shannon}

    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as staging:
        results = Path(staging) / "results"
        results.mkdir()
        write_table(results / "verified-groups.tsv", ["group", "n_samples"],
                    sorted(Counter(groups.values()).items()))
        write_table(results / "shannon-by-sample.tsv", ["sample_id", "group", "shannon"],
                    ((sample, groups[sample], row["shannon"]) for sample, row in diversity.items()))
        for filename, profiles in [
            ("composition-by-group.tsv", proportions),
            ("shannon-by-group.tsv", diversity),
            ("functional-features-by-group.tsv", functions),
        ]:
            write_table(results / filename, ["group", "feature", "n_samples", "median"],
                        summarize(groups, profiles))
        write_figures(results / "figures", groups, proportions, diversity, functions)
        results.rename(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, help="Directory containing the three input TSVs")
    parser.add_argument("output", type=Path, nargs="?", help="New directory for tables and figures")
    parser.add_argument("--verify-only", action="store_true", help="Check input files and exit")
    args = parser.parse_args()
    if args.verify_only == (args.output is not None):
        parser.error("Supply an output directory or use --verify-only.")
    try:
        if args.verify_only:
            groups, _, _ = verify_inputs(args.inputs)
            print(f"Verified {len(groups)} matched samples. Groups: {dict(Counter(groups.values()))}")
        else:
            run(args.inputs, args.output)
            print(f"Complete: {args.output}")
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(2, f"Workflow stopped: {error}\n")
