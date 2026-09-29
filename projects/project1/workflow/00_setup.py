from pathlib import Path
from xenium_workflow.project import load_project

project = load_project(Path(__file__).resolve().parents[1])
for path in (project.PROCESSED_DIR, project.RESULTS_DIR, project.LOGS_DIR):
    path.mkdir(parents=True, exist_ok=True)
print(f"Project ready: {project.PROJECT_DIR}")
