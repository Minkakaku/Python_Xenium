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


# 功能：构建 运行与输出 所需对象或目录（定义：prepare_project_dirs）。
def prepare_project_dirs(project):
    for path in (
        project.RAW_DIR,
        project.PROCESSED_DIR,
        project.REFERENCE_DIR,
        project.RESULTS_DIR,
        project.LOGS_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
    for name in REPORT_STEPS:
        (project.RESULTS_DIR / name).mkdir(parents=True, exist_ok=True)


# 功能：保存 运行与输出 产生的对象、表格或图形（定义：save_figure）。
def save_figure(path_no_ext, fig=None, dpi=300, close=True):
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig = fig or plt.gcf()
    fig.savefig(f"{path_no_ext}.png", dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", facecolor="white", bbox_inches="tight")
    if close:
        plt.close(fig)


# 功能：为 CON、EXP 及其他分组生成稳定配色（定义：group_palette）。
def group_palette(groups, configured=None):
    configured = configured or {}
    fallback = plt.get_cmap("tab10")
    return {
        group: configured.get(group, fallback(i % 10))
        for i, group in enumerate(groups)
    }
