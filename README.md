# Python-xenium-workflow-template

这是从原始 Xenium notebook 整理出的分层 Python 主流程。公共方法、公共默认和具体项目运行记录分开管理；原始代码保存在 `legacy/original/`，正式分析只使用 Python 脚本。

## 项目结构

```text
├── src/python/xenium_workflow/   # 跨项目复用方法与标准步骤
├── src/python/xenium_workflow/cluster_gpu.py  # 独立 RAPIDS GPU 聚类入口
├── config/                       # 公共默认与参考配置
├── projects/project1/            # 项目配置、workflow、数据、结果和日志
├── envs/                         # CPU/GPU 公共环境
├── tests/
├── source_catalog/               # 其他人员代码参考档案
└── legacy/original/              # 未修改的原始代码
```

项目 `config/config.yaml` 按同名键覆盖根目录 `config/default.yaml` 和 `config/reference.yaml`。样本、ROI、QC 阈值和人工注释均属于具体项目。

## 安装

```bash
conda env create -f envs/environment.yaml
conda activate xenium-workflow
pip install -e .
```

GPU 聚类使用 `envs/environment-gpu.yaml`；`python -m xenium_workflow.cluster_gpu` 是独立入口，默认使用第 0 张 GPU，可通过项目配置修改。

## 新项目开始

1. 复制 `projects/project1/` 为新的项目目录。
2. 修改项目 `config/samples.tsv` 和 `config/config.yaml`。
3. 将 Xenium Ranger 输出放入项目 `data/raw/<raw>/`。
4. 从 `workflow/00_setup.py` 开始逐步运行。
5. QC 后确认 `config/qc_thresholds.tsv`；marker 后填写 `config/celltype_annotation.tsv`。

```bash
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
python projects/project1/workflow/02_qc.py
```

已运行的 workflow 不覆盖，改参或重跑时新增下一编号脚本。完整顺序见 [新项目指南](docs/new_project.md)。

## 重要边界

- GPU 聚类输入 QC 后对象，输出项目 `data/processed/03_adata_clustered.h5ad`。
- 所有空间邻域只在同一样本内计算，避免不同 FOV 坐标重叠造成假邻居。
- 组间检验使用样本级汇总，不把细胞当作独立生物学重复。
- Linux 内部路径、项目编号、客户信息和固定 ROI 不进入复用代码。
- 正式图输出同名 PNG/PDF，绘图输入表保存 TSV/CSV。

`source_catalog/` 只作为分割、邻域和空间统计的方法参考，不作为正式运行入口。
