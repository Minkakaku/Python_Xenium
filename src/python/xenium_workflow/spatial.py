from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

from matplotlib.colors import LinearSegmentedColormap

from .runtime import NPG_DIVERGING, group_palette, save_figure, set_adata_palette


# 返回细胞类型的稳定顺序，分类列使用既有类别顺序。
def _categories(values):
    if isinstance(values.dtype, pd.CategoricalDtype):
        return [str(x) for x in values.cat.categories]
    return sorted(pd.unique(values.astype(str)).tolist())


# 功能：计算并汇总 空间统计与邻域 指标（定义：_sample_group_table）。
def _sample_group_table(adata, sample_col="sample", group_col="group"):
    return adata.obs[[sample_col, group_col]].drop_duplicates(sample_col).set_index(sample_col)


# 将原始 P 值转换为图中的显著性符号；缺失值显示 NA。
def _p_label(p):
    if pd.isna(p):
        return "NA"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# 对传入的样本汇总值做两组 Mann–Whitney U 检验，返回未校正 P 值。
def _two_group_tests(table, value_col, key_cols, group_col="group"):
    groups = table[group_col].drop_duplicates().tolist()
    rows = []
    if len(groups) != 2:
        return pd.DataFrame(columns=[*key_cols, "group_1", "group_2", "p_value", "significance"])
    for keys, data in table.groupby(key_cols, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        x = data.loc[data[group_col] == groups[0], value_col].dropna()
        y = data.loc[data[group_col] == groups[1], value_col].dropna()
        p_value = mannwhitneyu(x, y).pvalue if len(x) and len(y) else np.nan
        rows.append(
            dict(zip(key_cols, keys))
            | {
                "group_1": groups[0],
                "group_2": groups[1],
                "p_value": p_value,
                "significance": _p_label(p_value),
            }
        )
    return pd.DataFrame(rows)


# 逐样本建立空间图并计算细胞类型中心性，保存样本层面的表和组别展示。
def spatial_statistics(adata, out_dir, cluster_key="celltype", group_colors=None):
    import seaborn as sns
    import squidpy as sq

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    group_map = _sample_group_table(adata)
    all_scores = []

    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample].copy()
        sq.gr.spatial_neighbors(tmp, coord_type="generic", n_neighs=10, delaunay=True)
        sq.gr.centrality_scores(tmp, cluster_key=cluster_key)
        sq.pl.centrality_scores(tmp, cluster_key=cluster_key, figsize=(16, 5))
        save_figure(out_dir / f"{sample}_centrality_scores")

        scores = tmp.uns[f"{cluster_key}_centrality_scores"].copy()
        scores.index.name = cluster_key
        scores = scores.reset_index()
        scores["sample"] = sample
        scores["group"] = group_map.loc[sample, "group"]
        all_scores.append(scores)

    scores = pd.concat(all_scores, ignore_index=True)
    scores.to_csv(out_dir / "centrality_scores_by_sample.tsv", sep="\t", index=False)
    metrics = [x for x in scores.columns if x not in {cluster_key, "sample", "group"}]
    palette = group_palette(scores["group"].drop_duplicates(), group_colors)
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(max(8, scores[cluster_key].nunique() * 0.7), 5))
        sns.boxplot(data=scores, x=cluster_key, y=metric, hue="group", palette=palette, ax=ax, showfliers=False)
        sns.stripplot(data=scores, x=cluster_key, y=metric, hue="group", palette=palette, dodge=True, ax=ax, legend=False)
        ax.grid(False)
        ax.tick_params(axis="x", rotation=75)
        ax.set_xlabel("")
        ax.legend(title="", loc="upper left", bbox_to_anchor=(1.02, 1))
        save_figure(out_dir / f"centrality_{metric}_by_group", fig)
    return scores


# 在各样本内部查询 k 个邻居，计算每个细胞邻域内各类型的比例。
def _knn_composition(adata, column, k=10, sample_col="sample"):
    categories = _categories(adata.obs[column])
    frames = []
    for sample in adata.obs[sample_col].drop_duplicates():
        tmp = adata[adata.obs[sample_col] == sample]
        coords = np.asarray(tmp.obsm["spatial"])
        if len(coords) < 2:
            continue
        # 多查询一个位置以去掉自身；小样本以实际细胞数限制邻居数。
        n_neighbors = min(k + 1, len(coords))
        indices = cKDTree(coords).query(coords, k=n_neighbors)[1]
        if indices.ndim == 1:
            indices = indices[:, None]
        indices = indices[:, 1:]
        onehot = pd.get_dummies(tmp.obs[column].astype(str)).reindex(columns=categories, fill_value=0).to_numpy()
        composition = onehot[indices].mean(axis=1)
        frame = pd.DataFrame(composition, index=tmp.obs_names, columns=categories)
        frame["sample"] = sample
        frames.append(frame)
    return pd.concat(frames), categories


# 聚类局部邻域组成得到 CN，输出类型富集、样本内 CN 比例及空间图。
def cell_neighbors(
    adata,
    column,
    k=10,
    n_clusters=10,
    out_dir=None,
    group_colors=None,
):
    import seaborn as sns
    import squidpy as sq
    from sklearn.cluster import MiniBatchKMeans

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    composition, categories = _knn_composition(adata, column, k=k)
    model = MiniBatchKMeans(n_clusters=n_clusters, random_state=0, n_init=10)
    labels = model.fit_predict(composition[categories])
    label_col = f"CNs_{n_clusters}"
    label_values = pd.Series([f"CN{x}" for x in labels], index=composition.index)
    adata.obs[label_col] = pd.Categorical(label_values.reindex(adata.obs_names))
    set_adata_palette(adata, label_col)

    composition.assign(**{label_col: label_values}).to_csv(
        out_dir / f"CNs_{n_clusters}_neighbor_composition.tsv",
        sep="\t",
    )
    centers = pd.DataFrame(model.cluster_centers_, columns=categories, index=[f"CN{x}" for x in range(n_clusters)])
    tissue = pd.get_dummies(adata.obs[column].astype(str)).reindex(columns=categories, fill_value=0).mean()
    # 用全对象的细胞类型比例作背景；伪计数避免零比例导致无穷值。
    fold_change = np.log2((centers + 1e-6).div(tissue + 1e-6, axis=1))
    fold_change.to_csv(out_dir / f"CNs_{n_clusters}_log2_enrichment.tsv", sep="\t")
    grid = sns.clustermap(
        fold_change,
        vmin=-2,
        vmax=2,
        cmap=LinearSegmentedColormap.from_list("npg_diverging", NPG_DIVERGING),
        row_cluster=False,
        linewidths=0,
    )
    save_figure(out_dir / f"CNs_{n_clusters}_log2_enrichment", grid.fig)

    sample_groups = _sample_group_table(adata).reset_index()
    proportions = pd.crosstab(adata.obs["sample"], adata.obs[label_col], normalize="index")
    long = proportions.reset_index().melt(id_vars="sample", var_name="CN", value_name="proportion")
    long = long.merge(sample_groups, on="sample", how="left")
    long.to_csv(out_dir / f"CNs_{n_clusters}_proportion_by_sample.tsv", sep="\t", index=False)

    palette = group_palette(long["group"].drop_duplicates(), group_colors)
    fig, ax = plt.subplots(figsize=(max(8, n_clusters * 0.7), 6))
    sns.boxplot(data=long, x="CN", y="proportion", hue="group", palette=palette, showfliers=False, ax=ax)
    sns.stripplot(data=long, x="CN", y="proportion", hue="group", palette=palette, dodge=True, ax=ax, legend=False)
    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("Cell proportion")
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(out_dir / f"CNs_{n_clusters}_cell_proportion", fig)

    for sample in adata.obs["sample"].drop_duplicates():
        sq.pl.spatial_scatter(
            adata[adata.obs["sample"] == sample].copy(),
            color=label_col,
            title=f"{sample} CNs",
            shape=None,
            img=False,
        )
        save_figure(out_dir / f"{sample}_spatial_CNs_{n_clusters}")
    return adata


# 功能：把细胞级邻近结果汇总到样本层面（定义：proximity_sample_summary）。
def proximity_sample_summary(
    adata,
    radius=50,
    distance_bins=(0, 50, 100, 150, 200),
    celltype_col="celltype",
):
    categories = _categories(adata.obs[celltype_col])
    radius_rows = []
    annulus_rows = []

    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        coords = np.asarray(tmp.obsm["spatial"])
        types = tmp.obs[celltype_col].astype(str).to_numpy()
        group = tmp.obs["group"].iloc[0]
        tree = cKDTree(coords)
        radius_neighbors = tree.query_ball_point(coords, radius)
        max_neighbors = tree.query_ball_point(coords, max(distance_bins))

        for source in categories:
            source_indices = np.flatnonzero(types == source)
            if not len(source_indices):
                continue
            radius_matrix = np.zeros((len(source_indices), len(categories)), dtype=float)
            annulus_matrices = {
                (low, high): np.zeros((len(source_indices), len(categories)), dtype=float)
                for low, high in zip(distance_bins[:-1], distance_bins[1:])
            }

            for row, source_idx in enumerate(source_indices):
                neighbors = [x for x in radius_neighbors[source_idx] if x != source_idx]
                if neighbors:
                    counts = pd.Series(types[neighbors]).value_counts()
                    radius_matrix[row] = [counts.get(target, 0) / len(neighbors) for target in categories]

                candidates = np.array([x for x in max_neighbors[source_idx] if x != source_idx], dtype=int)
                if not len(candidates):
                    continue
                distances = np.linalg.norm(coords[candidates] - coords[source_idx], axis=1)
                for low, high in zip(distance_bins[:-1], distance_bins[1:]):
                    selected = candidates[(distances > low) & (distances <= high)]
                    counts = pd.Series(types[selected]).value_counts()
                    annulus_matrices[(low, high)][row] = [counts.get(target, 0) for target in categories]

            for target, value in zip(categories, radius_matrix.mean(axis=0)):
                radius_rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "source_celltype": source,
                        "target_celltype": target,
                        "proportion": value,
                        "source_cells": len(source_indices),
                    }
                )
            for (low, high), matrix in annulus_matrices.items():
                for target, value in zip(categories, np.median(matrix, axis=0)):
                    annulus_rows.append(
                        {
                            "sample": sample,
                            "group": group,
                            "source_celltype": source,
                            "target_celltype": target,
                            "distance_range": f"{low}-{high}",
                            "distance_mid": (low + high) / 2,
                            "median_count": value,
                            "source_cells": len(source_indices),
                        }
                    )
    return pd.DataFrame(radius_rows), pd.DataFrame(annulus_rows)


# 汇总指定半径和环带内的邻近关系，在样本层面比较两组并绘图。
def proximity_analysis(
    adata,
    out_dir,
    radius=50,
    distance_bins=(0, 50, 100, 150, 200),
    group_colors=None,
):
    import seaborn as sns

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    radius_table, annulus_table = proximity_sample_summary(adata, radius, distance_bins)
    radius_table.to_csv(out_dir / f"proximity_{radius}um_by_sample.tsv", sep="\t", index=False)
    annulus_table.to_csv(out_dir / "proximity_distance_bins_by_sample.tsv", sep="\t", index=False)
    tests = _two_group_tests(radius_table, "proportion", ["source_celltype", "target_celltype"])
    tests.to_csv(out_dir / f"proximity_{radius}um_group_tests.tsv", sep="\t", index=False)
    palette = group_palette(radius_table["group"].drop_duplicates(), group_colors)

    for source in radius_table["source_celltype"].drop_duplicates():
        source_dir = out_dir / str(source)
        source_dir.mkdir(parents=True, exist_ok=True)
        data = radius_table[radius_table["source_celltype"] == source]
        order = data["target_celltype"].drop_duplicates().tolist()
        fig, ax = plt.subplots(figsize=(max(10, len(order) * 0.65), 6))
        sns.boxplot(data=data, x="target_celltype", y="proportion", hue="group", order=order, palette=palette, showfliers=False, ax=ax)
        sns.stripplot(data=data, x="target_celltype", y="proportion", hue="group", order=order, palette=palette, dodge=True, ax=ax, legend=False)
        source_tests = tests[tests["source_celltype"] == source].set_index("target_celltype")
        ymax = data["proportion"].max()
        for x, target in enumerate(order):
            if target in source_tests.index:
                ax.text(x, ymax * 1.08 if ymax else 0.02, source_tests.loc[target, "significance"], ha="center")
        ax.set_ylim(top=(ymax * 1.18 if ymax else 0.05))
        ax.set_xlabel("")
        ax.set_ylabel(f"Neighbor proportion within {radius} µm")
        ax.tick_params(axis="x", rotation=75)
        ax.grid(False)
        ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left")
        save_figure(source_dir / f"proximity_{radius}um", fig)

        annulus = annulus_table[annulus_table["source_celltype"] == source]
        grid = sns.relplot(
            data=annulus,
            x="distance_mid",
            y="median_count",
            hue="group",
            col="target_celltype",
            col_wrap=4,
            kind="line",
            marker="o",
            errorbar="se",
            palette=palette,
            facet_kws={"sharey": False},
        )
        grid.set_axis_labels("Distance from source cell (µm)", "Median neighboring-cell count")
        save_figure(source_dir / "proximity_distance_ranges", grid.fig)
    return radius_table, annulus_table, tests


# 在样本内寻找每个细胞的最近一个邻居，记录其类型和距离。
def nearest_neighbor_table(adata):
    rows = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        if tmp.n_obs < 2:
            continue
        coords = np.asarray(tmp.obsm["spatial"])
        distances, indices = cKDTree(coords).query(coords, k=2)
        source_types = tmp.obs["celltype"].astype(str).to_numpy()
        groups = tmp.obs["group"].astype(str).to_numpy()
        for i in range(tmp.n_obs):
            rows.append(
                {
                    "cell_id": tmp.obs_names[i],
                    "source_type": source_types[i],
                    "neighbor_type": source_types[indices[i, 1]],
                    "distance": distances[i, 1],
                    "sample": sample,
                    "group": groups[i],
                }
            )
    return pd.DataFrame(rows)


# 功能：绘制并整理 空间统计与邻域 分析结果（定义：nearest_neighbor_boxplot）。
def nearest_neighbor_boxplot(summary, tests, out_dir, group_colors=None):
    import seaborn as sns

    palette = group_palette(summary["group"].drop_duplicates(), group_colors)
    for source in summary["source_type"].drop_duplicates():
        data = summary[summary["source_type"] == source]
        targets = data["neighbor_type"].drop_duplicates().tolist()
        ncols = 4
        nrows = int(np.ceil(len(targets) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 4), squeeze=False)
        source_tests = tests[tests["source_type"] == source].set_index("neighbor_type")
        for ax, target in zip(axes.ravel(), targets):
            panel = data[data["neighbor_type"] == target]
            sns.boxplot(data=panel, x="group", y="median_distance", hue="group", palette=palette, showfliers=False, ax=ax, legend=False)
            sns.stripplot(data=panel, x="group", y="median_distance", hue="group", palette=palette, ax=ax, legend=False)
            ymax = panel["median_distance"].max()
            label = source_tests.loc[target, "significance"] if target in source_tests.index else "NA"
            ax.text(0.5, ymax * 1.08 if ymax else 0.2, label, ha="center")
            ax.set_ylim(top=(ymax * 1.18 if ymax else 1))
            ax.set_title(str(target))
            ax.set_xlabel("")
            ax.set_ylabel("Median nearest-neighbor distance (µm)")
            ax.grid(False)
        for ax in axes.ravel()[len(targets):]:
            ax.remove()
        fig.suptitle(f"Source: {source}", y=1.01)
        fig.tight_layout()
        source_dir = Path(out_dir) / str(source)
        source_dir.mkdir(parents=True, exist_ok=True)
        save_figure(source_dir / "nearest_neighbor_distance_by_group", fig)


# 按样本及来源/邻居类型汇总最近邻距离中位数，再做组间检验和展示。
def nearest_neighbor(adata, out_dir, group_colors=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    distances = nearest_neighbor_table(adata)
    summary = (
        distances.groupby(["sample", "group", "source_type", "neighbor_type"], observed=True)["distance"]
        .median()
        .rename("median_distance")
        .reset_index()
    )
    tests = _two_group_tests(summary, "median_distance", ["source_type", "neighbor_type"])
    distances.to_csv(out_dir / "nearest_neighbor_cell_distances.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "nearest_neighbor_sample_medians.tsv", sep="\t", index=False)
    tests.to_csv(out_dir / "nearest_neighbor_group_tests.tsv", sep="\t", index=False)
    nearest_neighbor_boxplot(summary, tests, out_dir, group_colors)
    return distances, summary, tests
