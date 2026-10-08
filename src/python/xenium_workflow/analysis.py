from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import rc_context

from .runtime import (
    NPG_COLORS,
    NPG_SEQUENTIAL,
    categorical_palette,
    group_palette,
    save_figure,
    set_adata_palette,
    style_dotplot,
)


# 提取全体细胞的 UMAP 范围，供分组面板共用。
def _umap_limits(adata):
    coords = adata.obsm["X_umap"]
    return (coords[:, 0].min(), coords[:, 0].max()), (coords[:, 1].min(), coords[:, 1].max())


# 给一个或多个 UMAP 面板设置统一范围及等比例坐标。
def _fix_umap_axes(axes, xlim, ylim):
    for ax in np.atleast_1d(axes).ravel():
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")


# 功能：绘制并整理 降维聚类与 marker 分析结果（定义：cluster_plot）。
def cluster_plot(adata, res, out_dir, group_colors=None):
    import scanpy as sc
    import squidpy as sq

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xlim, ylim = _umap_limits(adata)
    set_adata_palette(adata, res)

    with rc_context({"figure.figsize": (6, 6)}):
        sc.pl.umap(
            adata,
            color=res,
            title="Cell Cluster",
            legend_loc="on data",
            size=8,
            legend_fontsize=9,
            legend_fontweight="normal",
            legend_fontoutline=2,
            show=False,
        )
        _fix_umap_axes(plt.gca(), xlim, ylim)
        save_figure(out_dir / "umap_dimplot_cluster")

    for sample in adata.obs["sample"].drop_duplicates():
        with rc_context({"figure.figsize": (8, 8)}):
            sq.pl.spatial_scatter(
                adata[adata.obs["sample"] == sample].copy(),
                color=res,
                title=str(sample),
                shape=None,
                size=1.5,
                img=False,
            )
            save_figure(out_dir / f"spatial_{sample}_cluster")

    groups = adata.obs["group"].drop_duplicates().tolist()
    if len(groups) > 1:
        fig, axes = plt.subplots(1, len(groups), figsize=(len(groups) * 5, 5), squeeze=False)
        for ax, group in zip(axes.ravel(), groups):
            sc.pl.umap(
                adata[adata.obs["group"] == group],
                color=res,
                title=str(group),
                ax=ax,
                size=8,
                legend_loc="on data",
                legend_fontsize=11,
                legend_fontweight="normal",
                legend_fontoutline=2,
                show=False,
            )
        _fix_umap_axes(axes, xlim, ylim)
        fig.tight_layout()
        save_figure(out_dir / "umap_split_group_cluster", fig)

    counts = pd.crosstab(adata.obs[res], adata.obs["sample"])
    counts.to_csv(out_dir / "cluster_statistics.tsv", sep="\t")
    adata.obs.to_csv(out_dir / "all_cluster_cell_features.tsv", sep="\t")
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(out_dir / "data_umap.tsv", sep="\t")

    proportions = pd.crosstab(adata.obs["sample"], adata.obs[res], normalize="index")
    palette = categorical_palette(proportions.columns)
    ax = proportions.plot(
        kind="bar",
        stacked=True,
        color=[palette[x] for x in proportions.columns],
        figsize=(max(6, len(proportions) * 0.75), 4.5),
        width=0.9,
    )
    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("Cell proportion")
    ax.set_title("Cluster cell proportions")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", ncol=2)
    save_figure(out_dir / "cluster_cell_proportion", ax.figure)
    return adata


# 基于 raw 表达检验 cluster marker，输出显著正向基因并可附加 GO 富集。
def marker_gene(adata, res, organism, out_dir, run_go=True, method="t-test"):
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.tl.rank_genes_groups(
        adata,
        groupby=res,
        use_raw=True,
        method=method,
        pts=True,
    )
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=res,
        standard_scale="var",
        n_genes=5,
        cmap=LinearSegmentedColormap.from_list("npg_dotplot", NPG_SEQUENTIAL),
        smallest_dot=18,
        show=False,
    )
    style_dotplot()
    save_figure(out_dir / "dotplot_top5_marker_gene")

    all_markers = []
    for cluster in adata.obs[res].cat.categories:
        cluster_dir = out_dir / f"cluster_{cluster}"
        cluster_dir.mkdir(parents=True, exist_ok=True)
        markers = sc.get.rank_genes_groups_df(adata, group=cluster)
        markers = markers[(markers["logfoldchanges"] > 0) & (markers["pvals_adj"] < 0.05)].copy()
        markers.to_csv(cluster_dir / "differentially_expressed_genes.tsv", sep="\t", index=False)

        if run_go and not markers.empty:
            try:
                _go_enrichment(markers["names"].tolist(), organism, cluster_dir)
            except Exception as exc:
                print(f"GO enrichment skipped for cluster {cluster}: {exc}")

        markers["cluster"] = str(cluster)
        all_markers.append(markers)

    all_markers = pd.concat(all_markers, ignore_index=True)
    all_markers = all_markers.rename(columns={"names": "Symbol"})
    all_markers.to_csv(out_dir / "all_cluster_markers.tsv", sep="\t", index=False)
    all_markers["cluster"].value_counts().sort_index().to_csv(
        out_dir / "cluster_marker_counts.tsv",
        sep="\t",
        header=["marker_count"],
    )
    return adata


# 调用 Enrichr 的 GO 2023 基因集，输出 BP/CC/MF 富集表及配对柱图。
def _go_enrichment(gene_list, organism, out_dir):
    import gseapy as gp

    background = "hsapiens_gene_ensembl" if organism == "Human" else "mmusculus_gene_ensembl"
    result = gp.enrichr(
        gene_list=gene_list,
        gene_sets=[
            "GO_Biological_Process_2023",
            "GO_Cellular_Component_2023",
            "GO_Molecular_Function_2023",
        ],
        background=background,
        organism=organism,
        cutoff=0.05,
        no_plot=True,
        outdir=str(out_dir),
        verbose=False,
    )
    table = result.results.copy()
    table["Term"] = table["Term"].str.replace(r" \(.*?\)", "", regex=True)
    table["Gene_set"] = table["Gene_set"].replace(
        {
            "GO_Biological_Process_2023": "BP",
            "GO_Cellular_Component_2023": "CC",
            "GO_Molecular_Function_2023": "MF",
        }
    )
    table = table.rename(columns={"Gene_set": "GO"})
    table.to_csv(Path(out_dir) / "GO_enrichment.tsv", sep="\t", index=False)
    gp.barplot(
        table,
        column="Adjusted P-value",
        group="GO",
        size=10,
        top_term=5,
        figsize=(3, 5),
        title="GO Enrichment",
        color=NPG_COLORS[:3],
        ofname=str(Path(out_dir) / "GO_enrichment.png"),
    )
    gp.barplot(
        table,
        column="Adjusted P-value",
        group="GO",
        size=10,
        top_term=5,
        figsize=(3, 5),
        title="GO Enrichment",
        color=NPG_COLORS[:3],
        ofname=str(Path(out_dir) / "GO_enrichment.pdf"),
    )


# 检查 cluster 映射完整性后写入细胞类型，输出 UMAP、空间图和样本组成。
def celltype_anno(adata, res, anno_dict, out_dir, group_colors=None):
    import scanpy as sc
    import squidpy as sq

    clusters = adata.obs[res].astype(str)
    missing = sorted(set(clusters) - set(anno_dict))
    if missing:
        raise ValueError(f"Missing cell-type annotations for clusters: {', '.join(missing)}")
    adata.obs["celltype"] = pd.Categorical(clusters.map(anno_dict))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xlim, ylim = _umap_limits(adata)
    set_adata_palette(adata, "celltype")
    with rc_context({"figure.figsize": (6, 6)}):
        sc.pl.umap(
            adata,
            color="celltype",
            title="Cell Types",
            legend_loc="right margin",
            legend_fontsize=9,
            size=8,
            show=False,
        )
        _fix_umap_axes(plt.gca(), xlim, ylim)
        save_figure(out_dir / "umap_dimplot_celltype")

    groups = adata.obs["group"].drop_duplicates().tolist()
    if len(groups) > 1:
        fig, axes = plt.subplots(1, len(groups), figsize=(len(groups) * 6, 6), squeeze=False)
        for ax, group in zip(axes.ravel(), groups):
            sc.pl.umap(
                adata[adata.obs["group"] == group],
                color="celltype",
                title=str(group),
                ax=ax,
                size=10,
                legend_loc="on data",
                legend_fontsize=11,
                legend_fontoutline=2,
                show=False,
            )
        _fix_umap_axes(axes, xlim, ylim)
        fig.tight_layout()
        save_figure(out_dir / "umap_dimplot_celltype_groups", fig)

    for sample in adata.obs["sample"].drop_duplicates():
        with rc_context({"figure.figsize": (8, 8)}):
            sq.pl.spatial_scatter(
                adata[adata.obs["sample"] == sample].copy(),
                color="celltype",
                title=str(sample),
                shape=None,
                size=1.5,
                img=False,
            )
            save_figure(out_dir / f"spatial_{sample}_celltype")

    adata.obs.to_csv(out_dir / "all_cell_features_celltype.tsv", sep="\t")
    proportions = pd.crosstab(adata.obs["sample"], adata.obs["celltype"], normalize="index")
    palette = categorical_palette(proportions.columns)
    ax = proportions.plot(
        kind="bar",
        stacked=True,
        color=[palette[x] for x in proportions.columns],
        figsize=(max(6, len(proportions) * 1.2), 4.5),
        width=0.9,
    )
    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("Cell proportion")
    ax.set_title("Cell type proportions")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(out_dir / "celltype_proportion", ax.figure)
    return adata
