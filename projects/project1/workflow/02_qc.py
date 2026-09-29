from pathlib import Path
from xenium_workflow.project import load_project
from xenium_workflow.steps import step_qc

project = load_project(Path(__file__).resolve().parents[1])
step_qc(project)
