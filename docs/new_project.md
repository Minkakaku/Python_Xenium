# 新项目使用指南

## 1. 建立项目

复制 `projects/project1/` 为新项目，修改项目 `config/samples.tsv` 和 `config/config.yaml`。样本表 `raw` 列对应项目 `data/raw/<raw>/` 下的 Xenium Ranger 输出。

若一个 FOV 通过 ROI 拆成多个样本，在项目 `config/ROI/` 放置对应 CSV，并把项目 `apply_roi` 设为 `true`。ROI 决策不得写入复用源码。

## 2. 逐步运行

```bash
python projects/my_project/workflow/00_setup.py
python projects/my_project/workflow/01_load.py
python projects/my_project/workflow/02_qc.py
```

检查并修改项目 `config/qc_thresholds.tsv` 后运行 `03_filter.py`。

## 3. GPU 聚类

```bash
conda activate xenium-gpu
pip install -e .
python projects/my_project/workflow/04_cluster_gpu.py
```

GPU 聚类完成后切回 CPU 环境，依次运行 marker、人工注释、空间统计、CN、proximity 和 nearest-neighbor。Enrichr 不可用时可在项目配置把 `run_go` 设为 `false`。

每个 workflow 文件只调用一个分析点；已运行脚本不覆盖。所有邻域按 `sample` 分开构图，CON/EXP 检验使用样本汇总。
