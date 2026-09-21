import pandas as pd

from .analysis import celltype_anno, cluster_plot, marker_gene
from .io import load_xenium_data
from .qc import cell_filter, data_qc, filter_parameter_record
from .runtime import prepare_project_dirs
from .spatial import cell_neighbors, nearest_neighbor, proximity_analysis, spatial_statistics


# 读取 Xenium 样本并按配置应用 ROI，保存合并后的原始 H5AD。
def step_load(project):
    prepare_project_dirs(project)
    adata = load_xenium_data(
        data_dir=project.RAW_DIR,
        sample_info=project.SAMPLE_INFO,
        roi_dir=project.ROI_DIR,
        apply_roi=project.APPLY_ROI,
        workers=project.LOAD_WORKERS,
        xenium_jobs=project.XENIUM_JOBS,
    )
    adata.write(project.ADATA_RAW)
    return adata


# 保存过滤前 QC 对象与图表，并生成供人工检查的逐样本阈值表。
def step_qc(project):
    import scanpy as sc

    prepare_project_dirs(project)
    adata = sc.read_h5ad(project.ADATA_RAW)
    adata = data_qc(adata, project.RESULTS_DIR / "01_Data_QC")
    adata.write(project.ADATA_QC)
    filter_parameter_record(
        adata,
        project.QC_THRESHOLDS,
        project.QC_MIN_GENE_QUANTILE,
        project.QC_MAX_GENE_QUANTILE,
        project.QC_MIN_COMPLEXITY,
        project.QC_MAX_COMPLEXITY,
    )
    return adata


# 读取 QC 对象和人工阈值表，保存过滤后的细胞 H5AD。
def step_filter(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_QC)
    adata = cell_filter(adata, project.QC_THRESHOLDS, project.RESULTS_DIR / "02_Cell_filter")
    adata.write(project.ADATA_FILTERED)
    return adata


# 读取独立聚类入口的结果，输出聚类展示、marker 和可选 GO 富集。
def step_cluster_plots_and_markers(project, run_go=True):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    adata = cluster_plot(
        adata,
        project.RESOLUTION_KEY,
        project.RESULTS_DIR / "03_Dimension_reduction_cluster",
        project.GROUP_COLORS,
    )
    adata = marker_gene(
        adata,
        project.RESOLUTION_KEY,
        project.ORGANISM,
        project.RESULTS_DIR / "04_Marker_gene_analysis",
        run_go=run_go,
    )
    return adata


# 应用项目中的 cluster-to-celltype 表并保存注释对象和展示图。
def step_annotation(project):
    import scanpy as sc

    mapping = pd.read_csv(project.ANNOTATION_FILE, sep="\t", dtype=str).dropna()
    anno_dict = dict(zip(mapping["cluster"], mapping["celltype"]))
    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    adata = celltype_anno(
        adata,
        project.RESOLUTION_KEY,
        anno_dict,
        project.RESULTS_DIR / "05_Cell_type_annotation",
        project.GROUP_COLORS,
    )
    adata.write(project.ADATA_ANNOTATED)
    return adata


# 读取注释对象，逐样本计算空间中心性并输出结果。
def step_spatial_statistics(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return spatial_statistics(
        adata,
        project.RESULTS_DIR / "06_Spatial_statistics",
        group_colors=project.GROUP_COLORS,
    )


# 按配置尝试不同 CN 数量，并保存带邻域标签的中间对象。
def step_cell_neighborhood(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    for n_clusters in project.CN_CLUSTERS:
        adata = cell_neighbors(
            adata,
            "celltype",
            k=project.CN_NEIGHBORS,
            n_clusters=n_clusters,
            out_dir=project.RESULTS_DIR / "07_Cell_neighborhood",
            group_colors=project.GROUP_COLORS,
        )
    adata.write(project.PROCESSED_DIR / "05_adata_with_neighborhoods.h5ad")
    return adata


# 读取注释对象，按配置的半径和环带计算样本级邻近关系。
def step_proximity(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return proximity_analysis(
        adata,
        project.RESULTS_DIR / "08_Proximity_analysis",
        radius=project.PROXIMITY_RADIUS,
        distance_bins=project.PROXIMITY_BINS,
        group_colors=project.GROUP_COLORS,
    )


# 读取注释对象，输出最近邻距离及样本级组间比较。
def step_nearest_neighbor(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return nearest_neighbor(
        adata,
        project.RESULTS_DIR / "09_Nearest_neighbor",
        group_colors=project.GROUP_COLORS,
    )
