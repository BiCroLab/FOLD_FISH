#!/usr/bin/env python3

### This script computes chromosome territory volumes by convex hull for each nucleus and channel in the provided datasets. It also generates summary statistics and plots.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import ConvexHull, QhullError


# =========================
# SETTINGS
# =========================

DATASETS = [
    {
        "name": "100X",
        "dots_file": "/group/bienko/projects/FFISH/Giulia_20251217_15min/pk1_15min_Allchannels_NewAcq_100X/05.final/pk1_15min_Allchannels_NewAcq_100X_all_channels_cwdots_20260429.tsv",
    },
    {
        "name": "100X001",
        "dots_file": "/group/bienko/projects/FFISH/Giulia_20251217_15min/pk1_15min_Allchannels_NewAcq_100X001/05.final/pk1_15min_Allchannels_NewAcq_100X001_all_channels_cwdots_20260429.tsv",
    },
    {
        "name": "100X002",
        "dots_file": "/group/bienko/projects/FFISH/Giulia_20251217_15min/pk1_15min_Allchannels_NewAcq_100X002/05.final/pk1_15min_Allchannels_NewAcq_100X002_all_channels_cwdots_20260429.tsv",
    },
]

OUTPUT_DIR = "/group/bienko/projects/FFISH/Giulia_20251217_15min/CT_cluster_volumes"
OUT_PREFIX = "CT_cluster_volumes"

CHANNEL_ORDER = ["chr17", "chr18", "chr19"]

CHANNEL_COLORS = {
    "chr17": "#0000ff",
    "chr18": "#bb0a1e",
    "chr19": "#008000",
}

CHANNEL_MAP = {
    "dw_594": "chr17",
    "dw_561": "chr18",
    "dw_630": "chr19",
}

X_UM = 0.065
Y_UM = 0.065
Z_UM = 0.200

REMOVE_VOLUME_OUTLIERS = True
OUTLIER_IQR_FACTOR = 1.5

REMOVE_N_DOTS_OUTLIERS = True
N_DOTS_MAX = 100


# =========================
# HELPERS
# =========================

def normalize_channel_name(ch):
    ch = str(ch)
    parts = ch.split("_")

    if len(parts) >= 2 and parts[0] == "dw":
        return "_".join(parts[:2])

    return ch


def map_channel_to_chromosome(ch):
    base = normalize_channel_name(ch)
    return CHANNEL_MAP.get(base, np.nan)


def channel_sort_key(ch):
    ch = str(ch)
    if ch in CHANNEL_ORDER:
        return CHANNEL_ORDER.index(ch)
    return len(CHANNEL_ORDER)


def ordered_channels(values):
    return sorted([str(v) for v in values], key=channel_sort_key)


# =========================
# VOLUME COMPUTATION
# =========================

def convex_hull_volume_um3(points_um):
    points_um = np.asarray(points_um, dtype=float)
    points_um = points_um[np.all(np.isfinite(points_um), axis=1)]

    if points_um.shape[0] < 4:
        return np.nan

    try:
        hull = ConvexHull(points_um)
        return float(hull.volume)
    except QhullError:
        return np.nan


def compute_cluster_volumes_for_dataset(dataset):
    name = dataset["name"]
    dots_file = dataset["dots_file"]

    if not os.path.isfile(dots_file):
        print(f"Missing dots file, skipping: {dots_file}")
        return None

    print(f"\nProcessing dataset: {name}")
    print(f"Dots file: {dots_file}")

    D = pd.read_csv(dots_file, sep="\t")

    required_cols = ["x", "y", "z", "fov", "channel", "nuc_id"]
    missing = [c for c in required_cols if c not in D.columns]

    if missing:
        raise ValueError(f"Missing required columns in {dots_file}: {missing}")

    D = D.copy()

    D["dataset"] = name
    D["fov"] = D["fov"].astype(str)
    D["nuc_id"] = D["nuc_id"].astype(int)

    D["raw_channel"] = D["channel"].astype(str)

    print("Raw channel counts before mapping:")
    print(D["raw_channel"].value_counts().sort_index())

    D["channel"] = D["raw_channel"].map(map_channel_to_chromosome)

    before = len(D)
    unmapped = D.loc[D["channel"].isna(), "raw_channel"].value_counts().sort_index()

    if len(unmapped) > 0:
        print("Unmapped raw channels, dropped:")
        print(unmapped)

    D = D.dropna(subset=["channel"]).copy()
    after = len(D)

    print(f"Mapped dots to chromosomes: {after}/{before}")
    print(f"Chromosomes found: {ordered_channels(D['channel'].unique())}")

    D["x_um"] = D["x"].astype(float) * X_UM
    D["y_um"] = D["y"].astype(float) * Y_UM
    D["z_um"] = D["z"].astype(float) * Z_UM

    rows = []

    group_cols = [
        "dataset",
        "fov",
        "nuc_id",
        "channel",
    ]

    for key, group in D.groupby(group_cols):
        dataset_name, fov, nuc_id, channel = key

        points_um = group[["x_um", "y_um", "z_um"]].to_numpy()
        volume_um3 = convex_hull_volume_um3(points_um)

        rows.append({
            "dataset": dataset_name,
            "fov": fov,
            "nuc_id": nuc_id,
            "channel": channel,
            "n_dots": len(group),
            "cluster_volume_um3": volume_um3,
        })

    V = pd.DataFrame(rows)

    V["nucleus_global_id"] = (
        V["dataset"].astype(str)
        + "_fov"
        + V["fov"].astype(str).str.zfill(3)
        + "_nuc"
        + V["nuc_id"].astype(int).astype(str)
    )

    return V


def compute_all_cluster_volumes():
    tables = []

    for ds in DATASETS:
        V = compute_cluster_volumes_for_dataset(ds)
        if V is not None:
            tables.append(V)

    if not tables:
        raise RuntimeError("No cluster-volume tables were created.")

    return pd.concat(tables, ignore_index=True)


# =========================
# OUTLIER REMOVAL
# =========================

def flag_volume_outliers(V):
    V = V.copy()
    V["volume_outlier"] = False

    for ch, idx in V.groupby("channel").groups.items():
        values = V.loc[idx, "cluster_volume_um3"].dropna()

        if len(values) < 4:
            continue

        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        lower = q1 - OUTLIER_IQR_FACTOR * iqr
        upper = q3 + OUTLIER_IQR_FACTOR * iqr

        outlier_idx = idx[
            (V.loc[idx, "cluster_volume_um3"] < lower) |
            (V.loc[idx, "cluster_volume_um3"] > upper)
        ]

        V.loc[outlier_idx, "volume_outlier"] = True

        print(
            f"{ch}: volume outliers = {len(outlier_idx)} "
            f"outside [{lower:.3f}, {upper:.3f}] um3"
        )

    return V


def flag_n_dots_outliers(V):
    V = V.copy()
    V["n_dots_outlier"] = V["n_dots"] > N_DOTS_MAX

    print(
        f"n_dots outliers = {int(V['n_dots_outlier'].sum())} "
        f"with n_dots > {N_DOTS_MAX}"
    )

    return V


def filter_outliers(V):
    V = V.copy()

    if REMOVE_VOLUME_OUTLIERS:
        V = flag_volume_outliers(V)
    else:
        V["volume_outlier"] = False

    if REMOVE_N_DOTS_OUTLIERS:
        V = flag_n_dots_outliers(V)
    else:
        V["n_dots_outlier"] = False

    keep = ~V["volume_outlier"] & ~V["n_dots_outlier"]
    return V.loc[keep].copy(), V


# =========================
# SUMMARY
# =========================

def make_summary(V):
    summary = (
        V
        .groupby("channel")
        .agg(
            n_clusters=("channel", "count"),
            n_datasets=("dataset", "nunique"),
            n_nuclei=("nucleus_global_id", "nunique"),
            mean_n_dots=("n_dots", "mean"),
            median_n_dots=("n_dots", "median"),
            mean_cluster_volume_um3=("cluster_volume_um3", "mean"),
            median_cluster_volume_um3=("cluster_volume_um3", "median"),
        )
        .reset_index()
    )

    summary["channel"] = pd.Categorical(
        summary["channel"],
        categories=CHANNEL_ORDER,
        ordered=True,
    )

    return summary.sort_values("channel").reset_index(drop=True)


# =========================
# PLOTTING
# =========================

def plot_combined_boxplot(
    V,
    value_col,
    ylabel,
    title,
    outname,
):
    plot_df = V.dropna(subset=[value_col]).copy()

    if plot_df.empty:
        print(f"No valid data for: {value_col}")
        return

    plot_df["channel"] = plot_df["channel"].astype(str)
    channels = ordered_channels(plot_df["channel"].unique())

    data_by_channel = {
        ch: plot_df.loc[plot_df["channel"] == ch, value_col].to_numpy(dtype=float)
        for ch in channels
    }

    fig, ax = plt.subplots(figsize=(7, 5))

    data = [data_by_channel[ch] for ch in channels]

    bp = ax.boxplot(
        data,
        showfliers=True,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2),
        boxprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.5),
        capprops=dict(color="black", linewidth=1.5),
        flierprops=dict(
            marker="o",
            markersize=3.0,
            markerfacecolor="black",
            markeredgecolor="black",
            alpha=0.6,
        ),
    )

    for patch, ch in zip(bp["boxes"], channels):
        patch.set_facecolor(CHANNEL_COLORS.get(ch, "gray"))
        patch.set_alpha(0.35)

    unit = "μm³" if value_col == "cluster_volume_um3" else "dots"

    labels = []
    for ch in channels:
        vals = data_by_channel[ch]
        labels.append(
            f"{ch}\n"
            f"$n$ = {len(vals)}\n"
            f"Median = {np.nanmedian(vals):.2f} {unit}"
        )

    ax.set_xticks(range(1, len(channels) + 1))
    ax.set_xticklabels(labels, fontsize=10)

    ax.set_xlabel("Chromosome", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(outname, format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved plot: {outname}")


def make_plots(V):
    plot_combined_boxplot(
        V,
        "cluster_volume_um3",
        "Cluster volume (μm³)",
        "Chromosome territory volume by convex hull",
        os.path.join(
            OUTPUT_DIR,
            f"{OUT_PREFIX}_convex_hull_volume_um3_combined.pdf",
        ),
    )

    plot_combined_boxplot(
        V,
        "n_dots",
        "Number of dots per chromosome territory",
        "Number of dots per chromosome territory",
        os.path.join(
            OUTPUT_DIR,
            f"{OUT_PREFIX}_n_dots_per_chromosome_territory_combined.pdf",
        ),
    )


# =========================
# MAIN
# =========================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Computing chromosome territory volumes by convex hull...")

    V_raw = compute_all_cluster_volumes()

    raw_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_raw.tsv",
    )
    V_raw.to_csv(raw_out, sep="\t", index=False)
    print(f"Saved raw volume table: {raw_out}")

    V, V_with_flags = filter_outliers(V_raw)

    # Keep only territories with valid convex-hull volume.
    # This makes the volume and n_dots plots use the same n.
    V = V.dropna(subset=["cluster_volume_um3"]).copy()

    flags_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_with_outlier_flags.tsv",
    )
    V_with_flags.to_csv(flags_out, sep="\t", index=False)
    print(f"Saved table with outlier flags: {flags_out}")

    print(f"\nTotal chromosome territories after filtering: {len(V)}")
    print(f"Datasets: {sorted(V['dataset'].unique())}")
    print(f"Nuclei: {V['nucleus_global_id'].nunique()}")
    print(f"Channels: {ordered_channels(V['channel'].unique())}")

    filtered_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_filtered.tsv",
    )
    V.to_csv(filtered_out, sep="\t", index=False)
    print(f"Saved filtered volume table: {filtered_out}")

    summary = make_summary(V)

    summary_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_summary_by_channel.tsv",
    )
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"Saved summary: {summary_out}")

    make_plots(V)

    print("Done!")


if __name__ == "__main__":
    main()
