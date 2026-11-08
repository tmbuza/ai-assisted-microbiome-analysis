"""Draw sample-level figures from validated metagenomics profiles."""
from collections import defaultdict
from statistics import mean, median


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
