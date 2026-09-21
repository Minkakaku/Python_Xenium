# 读取已聚类对象，输出聚类展示图和 marker 分析结果。

from pathlib import Path

from xenium_workflow.project import load_project
from xenium_workflow.steps import step_cluster_plots_and_markers

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
project = load_project(Path(__file__).resolve().parents[1])
step_cluster_plots_and_markers(project, run_go=project.RUN_GO)
