# 新项目分析顺序

先在项目根目录执行：

```bash
conda activate xenium-workflow
pip install -e .
```

每个实际运行步骤都在 `workflow/` 新建一个 `.py` 脚本，例如 `step01_load.py`、`step02_qc.py`。已运行的脚本不覆盖，不使用 notebook；在项目根目录用 `python -m workflow.step01_load` 这种方式运行。

## Step 01：读取 Xenium 数据

```python
from config import project
from xenium_workflow.steps import step_load

step_load(project)
```

## Step 02：QC 和阈值表

```python
from config import project
from xenium_workflow.steps import step_qc

step_qc(project)
```

检查并修改 `config/qc_thresholds.tsv`。

## Step 03：过滤

```python
from config import project
from xenium_workflow.steps import step_filter

step_filter(project)
```

## Step 04：GPU 聚类

切换到 `xenium-gpu` 环境并完成一次可编辑安装：

```bash
conda activate xenium-gpu
pip install -e .
```

新建 `workflow/step04_cluster_gpu.py`：

```python
from config import project
from cluster_gpu import run_gpu_clustering

run_gpu_clustering(
    input_path=project.ADATA_FILTERED,
    output_path=project.ADATA_CLUSTERED,
    resolutions=project.CLUSTER_RESOLUTIONS,
    device=project.GPU_DEVICE,
)
```

在项目根目录运行：

```bash
python -m workflow.step04_cluster_gpu
```

脚本保留原来的 RAPIDS、CuPy 和 RMM 路线。`GPU_DEVICE` 可指定显卡；GPU 聚类完成后切回 CPU 环境。

## Step 05：聚类图和 marker

```python
from config import project
from xenium_workflow.steps import step_cluster_plots_and_markers

step_cluster_plots_and_markers(project, run_go=True)
```

若当前机器不能访问 Enrichr，可先用 `run_go=False` 完成 marker。

## Step 06：人工注释

填写 `config/celltype_annotation.tsv` 后运行：

```python
from config import project
from xenium_workflow.steps import step_annotation

step_annotation(project)
```

## Step 07–10：空间分析

分别新建 workflow 调用：

```python
from config import project
from xenium_workflow.steps import (
    step_spatial_statistics,
    step_cell_neighborhood,
    step_proximity,
    step_nearest_neighbor,
)

# 每个 workflow 文件只保留其中一个调用
step_spatial_statistics(project)
step_cell_neighborhood(project)
step_proximity(project)
step_nearest_neighbor(project)
```

所有邻域计算都按 `sample` 分开构图；CON/EXP 检验使用每个样本的汇总值。
