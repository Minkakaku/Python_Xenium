# 调用独立 GPU 聚类入口；输入为过滤后对象，输出为聚类 H5AD。

from pathlib import Path

from xenium_workflow.cluster_gpu import run_gpu_clustering
from xenium_workflow.project import load_project

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
project = load_project(Path(__file__).resolve().parents[1])
run_gpu_clustering(
    input_path=project.ADATA_FILTERED,
    output_path=project.ADATA_CLUSTERED,
    resolutions=project.CLUSTER_RESOLUTIONS,
    device=project.GPU_DEVICE,
)
