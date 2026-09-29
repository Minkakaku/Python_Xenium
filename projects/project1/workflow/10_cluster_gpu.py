from pathlib import Path
from xenium_workflow.project import load_project
from xenium_workflow.cluster_gpu import run_gpu_clustering

project = load_project(Path(__file__).resolve().parents[1])
run_gpu_clustering(project.ADATA_FILTERED, project.ADATA_CLUSTERED,
                   resolutions=project.CLUSTER_RESOLUTIONS, device=project.GPU_DEVICE)
