#!/usr/bin/env python3

### This script computes the Euclidean distance transform (EDT) of chromosome territory centroids from the nuclear periphery using provided nuclear mask files. 
### It processes multiple datasets, adds EDT measurements to territory tables, generates summary statistics, and creates boxplots with significance testing.


import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import distance_transform_edt
from scipy.stats import mannwhitneyu


# =========================
# SETTINGS
# =========================

DATASETS = [
    {
        "name": "pps_001",
        "territory_file": "/group/bienko/projects/FFISH/Giulia_20251103_pepsin/001/06.chrom_territories/CT_territories.tsv",
        "nuclear_mask_pattern": "/group/bienko/projects/FFISH/Giulia_20251103_pepsin/001/04.stardist/375_{fov:03d}_mask2.npy",
    },

    {
        "name": "pps_002",
        "territory_file": "/group/bienko/projects/FFISH/Giulia_20251103_pepsin/002/06.chrom_territories/CT_territories.tsv",
        "nuclear_mask_pattern": "/group/bienko/projects/FFISH/Giulia_20251103_pepsin/002/04.stardist/375_{fov:03d}_mask2.npy",
    }

]

OUTPUT_DIR = "/group/bienko/projects/FFISH/Giulia_20251103_pepsin/CT_all_EDT"
OUT_PREFIX = "CT_all_EDT"

CHANNEL_ORDER = ["chr17", "chr18", "chr19"]

CHANNEL_COLORS = {
    "chr17": "#0000ff",
    "chr18": "#bb0a1e",
    "chr19": "#008000",
}


# pixel size in microns
X_UM = 0.1083
Y_UM = 0.1083
Z_UM = 0.300


# =========================
# COLUMN DETECTION
# =========================

def find_first_existing_column(df, candidates, required=True):
    for col in candidates:
        if col in df.columns:
            return col

    if required:
        raise ValueError(
            "None of these columns were found:\n"
            + "\n".join([f"  - {c}" for c in candidates])
            + "\nAvailable columns are:\n"
            + "\n".join([f"  - {c}" for c in df.columns])
        )

    return None


def get_centroid_columns_2d(df):
    y_col = find_first_existing_column(
        df,
        [
            "territory_y_px",
            "territory_projected_y_px",
            "centroid_y",
            "territory_centroid_y",
            "cluster_centroid_y",
            "y_centroid",
            "centroid_y_px",
            "centroid_y_vox",
            "center_y",
            "y",
        ],
    )

    x_col = find_first_existing_column(
        df,
        [
            "territory_x_px",
            "territory_projected_x_px",
            "centroid_x",
            "territory_centroid_x",
            "cluster_centroid_x",
            "x_centroid",
            "centroid_x_px",
            "centroid_x_vox",
            "center_x",
            "x",
        ],
    )

    return y_col, x_col


# =========================
# HELPERS
# =========================

def channel_sort_key(ch):
    ch = str(ch)
    if ch in CHANNEL_ORDER:
        return CHANNEL_ORDER.index(ch)
    return len(CHANNEL_ORDER)


def ordered_channels(values):
    values = [str(v) for v in values]
    return sorted(values, key=channel_sort_key)


def p_to_stars(p):
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_significance_bars(ax, data_by_channel, channels):
    pairs = list(itertools.combinations(range(len(channels)), 2))

    all_values = np.concatenate([
        np.asarray(data_by_channel[ch], dtype=float)
        for ch in channels
        if len(data_by_channel[ch]) > 0
    ])

    all_values = all_values[~np.isnan(all_values)]

    if all_values.size == 0:
        return

    y_max = float(np.nanmax(all_values))
    y_min = float(np.nanmin(all_values))
    y_range = y_max - y_min

    if y_range == 0:
        y_range = 1.0

    bar_height = y_max + 0.08 * y_range
    step = 0.12 * y_range

    for k, (i, j) in enumerate(pairs):
        ch1 = channels[i]
        ch2 = channels[j]

        v1 = np.asarray(data_by_channel[ch1], dtype=float)
        v2 = np.asarray(data_by_channel[ch2], dtype=float)

        v1 = v1[~np.isnan(v1)]
        v2 = v2[~np.isnan(v2)]

        if len(v1) < 2 or len(v2) < 2:
            label = "n/a"
        else:
            _, p = mannwhitneyu(v1, v2, alternative="two-sided")
            label = p_to_stars(p)

        x1 = i + 1
        x2 = j + 1
        y = bar_height + k * step
        h = 0.03 * y_range

        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], linewidth=1,color="black")
        ax.text(
            (x1 + x2) / 2,
            y + h,
            label,
            ha="center",
            va="bottom",
            fontsize=9,
        )

    ax.set_ylim(top=bar_height + len(pairs) * step + 0.15 * y_range)


# =========================
# EDT COMPUTATION
# =========================

def load_nuclear_mask(mask_file):
    mask = np.load(mask_file)

    if mask.ndim != 2:
        raise ValueError(
            f"Expected 2D mask2.npy, got shape {mask.shape}: {mask_file}"
        )

    return mask


def sample_edt_for_fov(group, nuclear_mask, y_col, x_col):
    result = group.copy()

    result["centroid_edt_distance_um"] = np.nan
    result["normalized_centroid_edt_distance"] = np.nan
    result["centroid_inside_nucleus_mask"] = False

    for nuc_id, nuc_df in result.groupby("nuc_id"):
        nucleus_mask = nuclear_mask == int(nuc_id)

        if not np.any(nucleus_mask):
            continue

        edt_um = distance_transform_edt(
            nucleus_mask,
            sampling=(Y_UM, X_UM),
        )

        max_edt_um = float(np.nanmax(edt_um))

        for idx, row in nuc_df.iterrows():
            y = int(round(row[y_col]))
            x = int(round(row[x_col]))

            if (
                y < 0 or y >= edt_um.shape[0] or
                x < 0 or x >= edt_um.shape[1]
            ):
                continue

            d = float(edt_um[y, x])

            result.loc[idx, "centroid_edt_distance_um"] = d
            result.loc[idx, "centroid_inside_nucleus_mask"] = bool(
                nucleus_mask[y, x]
            )

            if max_edt_um > 0:
                result.loc[idx, "normalized_centroid_edt_distance"] = (
                    d / max_edt_um
                )

    return result


def add_edt_to_territory_table(
    T,
    dataset_name,
    nuclear_mask_pattern,
):
    y_col, x_col = get_centroid_columns_2d(T)

    print(f"\nUsing centroid columns for {dataset_name}:")
    print(f"  y: {y_col}")
    print(f"  x: {x_col}")

    updated_groups = []

    for fov, group in T.groupby("fov"):
        fov_int = int(fov)

        mask_file = nuclear_mask_pattern.format(fov=fov_int)

        print(f"\nFOV {fov_int}")
        print(f"Mask file: {mask_file}")

        if not os.path.isfile(mask_file):
            print("  Missing mask2.npy -> EDT values set to NaN")

            g = group.copy()
            g["centroid_edt_distance_um"] = np.nan
            g["normalized_centroid_edt_distance"] = np.nan
            g["centroid_inside_nucleus_mask"] = False

            updated_groups.append(g)
            continue

        nuclear_mask = load_nuclear_mask(mask_file)

        g = sample_edt_for_fov(
            group,
            nuclear_mask=nuclear_mask,
            y_col=y_col,
            x_col=x_col,
        )

        n_valid = g["centroid_edt_distance_um"].notna().sum()

        print(f"  Territories: {len(g)}")
        print(f"  Valid EDT measurements: {n_valid}")

        updated_groups.append(g)

    return pd.concat(updated_groups, ignore_index=True)


def process_all_datasets():
    all_tables = []

    updated_dir = os.path.join(
        OUTPUT_DIR,
        "updated_CT_territories",
    )

    os.makedirs(updated_dir, exist_ok=True)

    for ds in DATASETS:
        name = ds["name"]

        territory_file = ds["territory_file"]
        nuclear_mask_pattern = ds["nuclear_mask_pattern"]

        if not os.path.isfile(territory_file):
            print(f"Missing territory file, skipping: {territory_file}")
            continue

        print("\n================================================")
        print(f"Processing dataset: {name}")
        print("================================================")

        T = pd.read_csv(territory_file, sep="\t")

        T["dataset"] = name
        T["source_file"] = territory_file

        T = add_edt_to_territory_table(
            T,
            dataset_name=name,
            nuclear_mask_pattern=nuclear_mask_pattern,
        )

        updated_out = os.path.join(
            updated_dir,
            f"{name}_CT_territories_with_EDT.tsv",
        )

        T.to_csv(updated_out, sep="\t", index=False)

        print(f"\nSaved updated EDT territory file:")
        print(updated_out)

        all_tables.append(T)

    if not all_tables:
        raise RuntimeError("No datasets were processed.")

    return pd.concat(all_tables, ignore_index=True)


# =========================
# SUMMARY
# =========================

def make_summary(T):
    agg_dict = {
        "n_territories": ("territory_id", "count"),
        "n_datasets": ("dataset", "nunique"),
        "n_nuclei": ("nucleus_global_id", "nunique"),
    }

    optional_metrics = [
        "centroid_edt_distance_um",
        "normalized_centroid_edt_distance",
    ]

    for col in optional_metrics:
        if col in T.columns:
            agg_dict[f"mean_{col}"] = (col, "mean")
            agg_dict[f"median_{col}"] = (col, "median")

    summary = (
        T
        .groupby("channel")
        .agg(**agg_dict)
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

def plot_boxplot_with_points(
    T,
    value_col,
    ylabel,
    title,
    outname,
    ylim=None,
):
    if value_col not in T.columns:
        print(f"Column missing, skipping EDT plot: {value_col}")
        return

    plot_df = T.dropna(subset=[value_col]).copy()

    if plot_df.empty:
        print(f"No valid EDT data for: {value_col}")
        return

    plot_df["channel"] = plot_df["channel"].astype(str)

    channels = ordered_channels(
        plot_df["channel"].unique()
    )

    data_by_channel = {
        ch: plot_df.loc[
            plot_df["channel"] == ch,
            value_col
        ].to_numpy()
        for ch in channels
    }

    nuclei_by_channel = (
        plot_df[["channel", "nucleus_global_id"]]
        .drop_duplicates()
        .groupby("channel")
        .size()
        .to_dict()
    )

    total_nuclei = plot_df["nucleus_global_id"].nunique()

    fig, ax = plt.subplots(figsize=(7, 5))

    data = [data_by_channel[ch] for ch in channels]

    bp = ax.boxplot(
        data,
        showfliers=False,
        patch_artist=True,
    )

    for patch, ch in zip(bp["boxes"], channels):
        patch.set_facecolor(
            CHANNEL_COLORS.get(ch, "gray")
        )
        patch.set_alpha(0.35)

    rng = np.random.default_rng(1)

    for i, ch in enumerate(channels, start=1):
        y = data_by_channel[ch]

        x = rng.normal(
            i,
            0.04,
            size=len(y),
        )

        ax.scatter(
            x,
            y,
            s=18,
            alpha=0.6,
            color=CHANNEL_COLORS.get(ch, "gray"),
        )

    ax.set_xticks(range(1, len(channels) + 1))

    ax.set_xticklabels(
        [
            (
                f"{ch}\n"
                f"nuclei={nuclei_by_channel.get(ch, 0)}\n"
                f"territories={len(data_by_channel[ch])}"
            )
            for ch in channels
        ],
        fontsize=8,
    )

    ax.set_xlabel("Chromosome")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    ax.text(
        0.98,
        0.98,
        f"Total nuclei = {total_nuclei}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.7,
        ),
    )

    if ylim is not None:
        ax.set_ylim(*ylim)

    add_significance_bars(
        ax,
        data_by_channel,
        channels,
    )

    plt.tight_layout()
    plt.savefig(outname)
    plt.close(fig)

    print(f"Saved EDT plot: {outname}")


def make_all_combined_plots(T):
    plot_boxplot_with_points(
        T,
        "centroid_edt_distance_um",
        "EDT distance from nuclear periphery (um)",
        "Territory centroid EDT distance from nuclear periphery",
        os.path.join(
            OUTPUT_DIR,
            f"{OUT_PREFIX}_centroid_EDT_distance_from_nuclear_periphery_um.pdf",
        ),
    )

    plot_boxplot_with_points(
        T,
        "normalized_centroid_edt_distance",
        "Normalized EDT distance from nuclear periphery",
        "Normalized territory centroid EDT distance from nuclear periphery",
        os.path.join(
            OUTPUT_DIR,
            f"{OUT_PREFIX}_normalized_centroid_EDT_distance_from_nuclear_periphery.pdf",
        ),
        ylim=(-0.05, 1.05),
    )


# =========================
# MAIN
# =========================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(
        "\nComputing territory centroid EDT distances "
        "from nuclear periphery using mask2.npy files..."
    )

    T = process_all_datasets()

    T["nucleus_global_id"] = (
        T["dataset"].astype(str)
        + "_fov"
        + T["fov"].astype(int).astype(str).str.zfill(3)
        + "_nuc"
        + T["nuc_id"].astype(int).astype(str)
    )

    print("\n================================================")
    print("FINAL SUMMARY")
    print("================================================")

    print(f"Total pooled territories: {len(T)}")

    print(
        f"Nuclei in pooled EDT analysis: "
        f"{T['nucleus_global_id'].nunique()}"
    )

    print(
        f"Channels: "
        f"{ordered_channels(T['channel'].unique())}"
    )

    print(f"Saving EDT outputs to: {OUTPUT_DIR}")

    combined_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_with_EDT.tsv",
    )

    T.to_csv(combined_out, sep="\t", index=False)

    print(f"\nSaved pooled EDT table:")
    print(combined_out)

    summary = make_summary(T)

    summary_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_EDT_summary_by_channel.tsv",
    )

    summary.to_csv(summary_out, sep="\t", index=False)

    print(f"\nSaved EDT summary:")
    print(summary_out)

    make_all_combined_plots(T)

    print("\nDone!")


if __name__ == "__main__":
    main()