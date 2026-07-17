from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import save_figure


def add_qc_metrics(adata):
    import scanpy as sc

    percent_top = tuple(x for x in (10, 50, 100) if x < adata.n_vars) or None
    sc.pp.calculate_qc_metrics(adata, percent_top=percent_top, inplace=True)
    genes = adata.obs["n_genes_by_counts"].to_numpy(dtype=float)
    counts = adata.obs["total_counts"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        adata.obs["log10GenesPerUMI"] = np.where(
            (genes > 0) & (counts > 1),
            np.log10(genes) / np.log10(counts),
            np.nan,
        )
    return adata


def _plot_qc_panel(adata, out_dir, prefix):
    import matplotlib.pyplot as plt
    import scanpy as sc

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "log10GenesPerUMI"],
        groupby="sample",
        jitter=False,
        rotation=90,
        multi_panel=True,
        show=False,
    )
    for ax in plt.gcf().get_axes():
        ax.grid(False)
    save_figure(Path(out_dir) / f"{prefix}_nCount_nFeature_violin")

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="sample",
        title="",
        legend_loc="none",
        ax=axes[0],
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="log10GenesPerUMI",
        color="sample",
        title="",
        ax=axes[1],
        show=False,
    )
    for ax in axes:
        ax.grid(False)
    fig.tight_layout()
    save_figure(Path(out_dir) / f"{prefix}_nCount_nFeature_scatter", fig)


def data_qc(adata, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = add_qc_metrics(adata)
    _plot_qc_panel(adata, out_dir, "raw")
    adata.obs.to_csv(out_dir / "all_cell_features_raw.tsv", sep="\t")

    for sample in adata.obs["sample"].drop_duplicates():
        sample_dir = out_dir / str(sample)
        sample_dir.mkdir(parents=True, exist_ok=True)
        tmp = adata[adata.obs["sample"] == sample].copy()
        _plot_qc_panel(tmp, sample_dir, "raw")
        tmp.obs.to_csv(sample_dir / "all_cell_features_raw.tsv", sep="\t")
    return adata


def filter_parameter_record(
    adata,
    output_file,
    min_gene_quantile=0.05,
    max_gene_quantile=0.995,
    min_complexity=0.8,
    max_complexity=0.99,
):
    rows = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        rows.append(
            {
                "sample": sample,
                "min_genes": int(np.quantile(tmp.obs["n_genes_by_counts"], min_gene_quantile)),
                "max_genes": int(np.quantile(tmp.obs["n_genes_by_counts"], max_gene_quantile)),
                "min_GenesPerUMI": min_complexity,
                "max_GenesPerUMI": max_complexity,
            }
        )
    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)


def cell_filter(adata_qc, para_file, out_dir):
    import anndata

    thresholds = pd.read_csv(para_file, sep="\t")
    objects = []
    for row in thresholds.itertuples(index=False):
        tmp = adata_qc[adata_qc.obs["sample"] == row.sample]
        keep = (
            (tmp.obs["n_genes_by_counts"] > row.min_genes)
            & (tmp.obs["n_genes_by_counts"] < row.max_genes)
            & (tmp.obs["log10GenesPerUMI"] > row.min_GenesPerUMI)
            & (tmp.obs["log10GenesPerUMI"] < row.max_GenesPerUMI)
        )
        objects.append(tmp[keep].copy())

    adata = anndata.concat(objects, join="inner", merge="same", index_unique=None)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _plot_qc_panel(adata, out_dir, "filtered")
    adata.obs.to_csv(out_dir / "all_cell_features_filtered.tsv", sep="\t")

    for sample in adata.obs["sample"].drop_duplicates():
        sample_dir = out_dir / str(sample)
        sample_dir.mkdir(parents=True, exist_ok=True)
        _plot_qc_panel(adata[adata.obs["sample"] == sample].copy(), sample_dir, "filtered")
    return adata
