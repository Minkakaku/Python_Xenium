from pathlib import Path

import matplotlib.pyplot as plt


REPORT_STEPS = (
    "01_Data_QC",
    "02_Cell_filter",
    "03_Dimension_reduction_cluster",
    "04_Marker_gene_analysis",
    "05_Cell_type_annotation",
    "06_Spatial_statistics",
    "07_Cell_neighborhood",
    "08_Proximity_analysis",
    "09_Nearest_neighbor",
)


def prepare_project_dirs(project):
    project.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    project.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    for name in REPORT_STEPS:
        (project.RESULTS_DIR / name).mkdir(parents=True, exist_ok=True)


def save_figure(path_no_ext, fig=None, dpi=300, close=True):
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig = fig or plt.gcf()
    fig.savefig(f"{path_no_ext}.png", dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", facecolor="white", bbox_inches="tight")
    if close:
        plt.close(fig)


def group_palette(groups, configured=None):
    configured = configured or {}
    fallback = plt.get_cmap("tab10")
    return {
        group: configured.get(group, fallback(i % 10))
        for i, group in enumerate(groups)
    }
