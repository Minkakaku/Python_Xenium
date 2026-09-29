import pandas as pd

from .analysis import celltype_anno, cluster_plot, marker_gene
from .io import load_xenium_data
from .qc import cell_filter, data_qc, filter_parameter_record
from .runtime import prepare_project_dirs
from .spatial import cell_neighbors, nearest_neighbor, proximity_analysis, spatial_statistics


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


def step_filter(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_QC)
    adata = cell_filter(adata, project.QC_THRESHOLDS, project.RESULTS_DIR / "02_Cell_filter")
    adata.write(project.ADATA_FILTERED)
    return adata


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


def step_spatial_statistics(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return spatial_statistics(
        adata,
        project.RESULTS_DIR / "06_Spatial_statistics",
        group_colors=project.GROUP_COLORS,
    )


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


def step_nearest_neighbor(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return nearest_neighbor(
        adata,
        project.RESULTS_DIR / "09_Nearest_neighbor",
        group_colors=project.GROUP_COLORS,
    )
