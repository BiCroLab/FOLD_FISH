#!/usr/bin/env python3

### This script analyzes chromosome territory data to compute and visualize the distribution of clusters per nucleus and dots per cluster across multiple datasets.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# SETTINGS
# =========================

TERRITORY_FILES = [
    "/group/bienko/projects/FFISH/Giulia_20260132_pk/005/06.chrom_territories/CT_territories.tsv",
    "/group/bienko/projects/FFISH/Giulia_20260132_pk/006/06.chrom_territories/CT_territories.tsv",
]

DATASET_NAMES = [
    "pk_005",
    "pk_006",
]

OUTPUT_DIR = "/group/bienko/projects/FFISH/Giulia_20260132_pk/CT_all_cluster_distributions"
OUT_PREFIX = "CT_all_cluster_distributions"

CHANNEL_ORDER = ["chr17", "chr18", "chr19"]

BAR_COLOR = "#6f746d"
EDGE_COLOR = "black"
MEDIAN_COLOR = "#ff9900"


# =========================
# LOAD DATA
# =========================

def load_all_territories(files, dataset_names):
    tables = []

    for fname, dataset in zip(files, dataset_names):
        if not os.path.isfile(fname):
            print(f"Missing file, skipping: {fname}")
            continue

        T = pd.read_csv(fname, sep="\t")
        T["dataset"] = dataset
        T["source_file"] = fname
        tables.append(T)

    if not tables:
        raise RuntimeError("No territory files were loaded.")

    return pd.concat(tables, ignore_index=True)


def add_global_ids(T):
    T = T.copy()

    T["nucleus_global_id"] = (
        T["dataset"].astype(str)
        + "_fov"
        + T["fov"].astype(int).astype(str).str.zfill(3)
        + "_nuc"
        + T["nuc_id"].astype(int).astype(str)
    )

    T["territory_global_id"] = (
        T["nucleus_global_id"]
        + "_"
        + T["channel"].astype(str)
        + "_territory"
        + T["territory_id"].astype(int).astype(str)
    )

    return T


# =========================
# PLOTS
# =========================

def plot_clusters_per_nucleus(T):
    counts = (
        T[["nucleus_global_id", "territory_global_id"]]
        .drop_duplicates()
        .groupby("nucleus_global_id")
        .size()
        .to_numpy()
    )

    n_nuclei = len(counts)

    if n_nuclei == 0:
        print("No nuclei found for clusters-per-nucleus plot.")
        return

    max_count = int(np.nanmax(counts))
    bins = np.arange(-0.5, max_count + 1.5, 1)

    weights = np.ones_like(counts, dtype=float) * 100.0 / n_nuclei

    fig, ax = plt.subplots(figsize=(4.2, 4.0))

    ax.hist(
        counts,
        bins=bins,
        weights=weights,
        color=BAR_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=1.2,
    )

    ax.axvline(
    np.nanmedian(counts),
    color=MEDIAN_COLOR,
    linestyle="--",
    linewidth=1.6,
    )  

    ax.set_xlabel("# of clusters per nucleus", fontsize=13)
    ax.set_ylabel("% of nuclei", fontsize=13)
    ax.set_title("", fontsize=14)

    ax.text(
        0.96,
        0.96,
        f"$n$ = {n_nuclei:,}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=13,
        style="italic",
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)
    ax.tick_params(axis="both", labelsize=12, width=1.5)

    ax.set_xticks(np.arange(0, max_count + 1, 1))

    out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_clusters_per_nucleus.pdf",
    )

    plt.tight_layout()
    plt.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out}")


def plot_dots_per_cluster(T):
    values = T["n_dots"].dropna().astype(float).to_numpy()

    n_clusters = len(values)

    if n_clusters == 0:
        print("No clusters found for dots-per-cluster plot.")
        return

    median_value = float(np.nanmedian(values))

    max_value = int(np.nanmax(values))
    bins = np.arange(0.5, max_value + 1.5, 1)

    weights = np.ones_like(values, dtype=float) * 100.0 / n_clusters

    fig, ax = plt.subplots(figsize=(4.2, 4.0))

    ax.hist(
        values,
        bins=bins,
        weights=weights,
        color=BAR_COLOR,
        edgecolor=EDGE_COLOR,
        linewidth=1.0,
    )

    ax.axvline(
        median_value,
        color=MEDIAN_COLOR,
        linestyle="--",
        linewidth=1.6,
    )

    ax.set_xlabel("# of dots per cluster", fontsize=13)
    ax.set_ylabel("% of clusters", fontsize=13)
    ax.set_title("", fontsize=14)

    ax.text(
        0.96,
        0.96,
        f"$n$ = {n_clusters:,}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=13,
        style="italic",
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)
    ax.tick_params(axis="both", labelsize=12, width=1.5)

    out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_dots_per_cluster.pdf",
    )

    plt.tight_layout()
    plt.savefig(out, format="pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {out}")


def save_tables(T):
    cluster_counts = (
        T[["nucleus_global_id", "territory_global_id"]]
        .drop_duplicates()
        .groupby("nucleus_global_id")
        .size()
        .reset_index(name="n_clusters_per_nucleus")
    )

    dots_per_cluster = (
        T[[
            "dataset",
            "fov",
            "nuc_id",
            "nucleus_global_id",
            "channel",
            "territory_id",
            "territory_global_id",
            "n_dots",
        ]]
        .drop_duplicates()
        .copy()
    )

    cluster_counts_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_clusters_per_nucleus.tsv",
    )

    dots_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_dots_per_cluster.tsv",
    )

    cluster_counts.to_csv(cluster_counts_out, sep="\t", index=False)
    dots_per_cluster.to_csv(dots_out, sep="\t", index=False)

    print(f"Saved: {cluster_counts_out}")
    print(f"Saved: {dots_out}")


# =========================
# MAIN
# =========================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading chromosome territory files...")
    T = load_all_territories(TERRITORY_FILES, DATASET_NAMES)
    T = add_global_ids(T)

    print(f"Total territories/clusters: {T['territory_global_id'].nunique()}")
    print(f"Total nuclei: {T['nucleus_global_id'].nunique()}")
    print(f"Channels: {sorted(T['channel'].astype(str).unique())}")
    print(f"Saving to: {OUTPUT_DIR}")

    pooled_out = os.path.join(
        OUTPUT_DIR,
        f"{OUT_PREFIX}_pooled_territories.tsv",
    )

    T.to_csv(pooled_out, sep="\t", index=False)
    print(f"Saved: {pooled_out}")

    save_tables(T)

    plot_clusters_per_nucleus(T)
    plot_dots_per_cluster(T)

    print("Done!")


if __name__ == "__main__":
    main()
