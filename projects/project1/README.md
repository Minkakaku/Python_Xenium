# project1

这是可复制的 Xenium 项目入口。样本、ROI、QC 阈值和人工注释放在 `config/`；短入口放在 `workflow/`；中间 H5AD 写入 `data/processed/`；最终图表写入 `results/`。

```bash
conda env create -f envs/environment.yaml
conda activate xenium-workflow
pip install -e .
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
```

QC 后人工确认 `config/qc_thresholds.tsv` 再过滤。GPU 聚类使用 `envs/environment-gpu.yaml` 和 `04_cluster_gpu.py`。marker 完成后填写 `config/celltype_annotation.tsv`，再运行注释与空间分析。
