# Xenium 空间转录组分析模板

这是从原始 Xenium 分析 notebook 中整理出的可复用项目模板。原始代码完整保存在 `legacy/original/`；`src/` 尽量沿用原函数名称和分析逻辑，只修复明确错误、项目硬编码及会影响统计结论的问题。

## 代码约定

- `src/`：可跨项目复用的 Xenium 函数与分析方法。
- `config/`：当前项目路径、样本信息、QC 参数和人工注释。
- `workflow/`：新项目开始时为空；每实际运行一个分析点，新建一个 `stepXX_name.py`。
- `legacy/original/`：原始 notebook 和脚本，只用于追溯，不直接修改。
- `data/`、`results/`、`logs/`：输入、中间对象和输出，不提交 Git。

## 环境

CPU 分析环境：

```bash
conda env create -f environment.yml
conda activate xenium-workflow
pip install -e .
```

GPU 聚类单独使用：

```bash
conda env create -f environment-gpu.yml
conda activate xenium-gpu
pip install -e .
```

GPU 环境要求可用的 NVIDIA 驱动和与环境匹配的 CUDA。
`cluster_gpu.py` 是聚类阶段单独调用的 GPU 加速入口，默认使用第 0 张 GPU；多卡机器可通过 `--device` 指定设备。

正式分析只使用 Python 脚本，不使用 notebook。`legacy/original/` 中的 notebook 仅为原始代码档案。

## 新项目开始

1. 复制 `config/sample_info.example.tsv` 为 `config/sample_info.tsv`。
2. 修改 `config/project.py`。
3. 将 Xenium Ranger 输出放入 `data/raw/<raw>/`，其中 `<raw>` 与样本表第一列一致。
4. 从下面的标准方法中选择一个，在 `workflow/` 新建本次实际运行的脚本。
步骤脚本使用下划线命名，并在项目根目录用模块方式运行，例如 `python -m workflow.step01_load`。

标准方法：

```python
from config import project
from xenium_workflow.steps import step_load, step_qc, step_filter

step_load(project)
step_qc(project)
# 人工检查 config/qc_thresholds.tsv 后再运行：
step_filter(project)
```

QC 完成后，在 GPU 环境新建 `workflow/step04_cluster_gpu.py`：

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

之后回到 CPU 环境运行 marker、人工注释和空间分析。完整顺序见 [新项目指南](docs/new_project.md)。

## 重要修正

- GPU 脚本读取 QC 后对象，并与后续默认 `res_0.3` 保持一致。
- 删除原 GPU 脚本末尾固定 `P1_CRC` 和裁剪坐标的 Visium 专属绘图。
- 所有空间邻域只在同一样本内部计算，避免不同 FOV 坐标重叠造成假邻居。
- 组间检验改为样本级汇总，不把每个细胞当作独立生物学重复。
- Linux 服务器路径、项目编号、客户信息和 `cp -r` 报告模板不进入通用代码。

详细审查记录见 [原代码审查](docs/original_code_review_CN.md)。
