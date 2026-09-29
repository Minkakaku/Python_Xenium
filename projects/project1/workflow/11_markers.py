from pathlib import Path
from xenium_workflow.project import load_project
from xenium_workflow.steps import step_cluster_plots_and_markers

project = load_project(Path(__file__).resolve().parents[1])
step_cluster_plots_and_markers(project)
