"""公共 YAML + 项目 YAML，项目路径独立于启动目录。"""
from pathlib import Path
from types import SimpleNamespace
import yaml

def merge_config(base, override):
    result = dict(base)
    for key, value in override.items():
        result[key] = merge_config(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else value
    return result

def load_project(project_dir):
    project_dir = Path(project_dir).resolve()
    workspace = project_dir.parents[1]
    settings = {}
    for path in (workspace / "config/default.yaml", workspace / "config/reference.yaml", project_dir / "config/config.yaml"):
        settings = merge_config(settings, yaml.safe_load(path.read_text()) or {})
    project = SimpleNamespace(**{key.upper(): value for key, value in settings.items()})
    project.WORKSPACE_DIR = workspace
    project.PROJECT_DIR = project_dir
    project.CONFIG_DIR = project.PROJECT_DIR / "config"
    project.DATA_DIR = project.PROJECT_DIR / "data"
    project.RAW_DIR = project.DATA_DIR / "raw"
    project.PROCESSED_DIR = project.DATA_DIR / "processed"
    project.RESULTS_DIR = project.PROJECT_DIR / "results"
    project.SAMPLE_INFO = project.CONFIG_DIR / "samples.tsv"
    project.ROI_DIR = project.CONFIG_DIR / "ROI"
    project.QC_THRESHOLDS = project.CONFIG_DIR / "qc_thresholds.tsv"
    project.ANNOTATION_FILE = project.CONFIG_DIR / "celltype_annotation.tsv"
    project.ADATA_RAW = project.PROCESSED_DIR / "00_adata_raw.h5ad"
    project.ADATA_QC = project.PROCESSED_DIR / "01_adata_qc.h5ad"
    project.ADATA_FILTERED = project.PROCESSED_DIR / "02_adata_filtered.h5ad"
    project.ADATA_CLUSTERED = project.PROCESSED_DIR / "03_adata_clustered.h5ad"
    project.ADATA_ANNOTATED = project.PROCESSED_DIR / "04_adata_annotated.h5ad"
    project.LOGS_DIR = project_dir / "logs"
    return project
