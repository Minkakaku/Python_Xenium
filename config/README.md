# config

新项目开始时：

1. 复制 `sample_info.example.tsv` 为 `sample_info.tsv`；
2. 修改 `project.py` 中的分组顺序、QC、聚类和空间分析参数；
3. marker 完成后复制 `celltype_annotation.example.tsv` 为 `celltype_annotation.tsv`，填写真实生物学注释。

`sample_info.tsv` 四列含义：

- `raw`：`data/raw/` 下的 Xenium Ranger 输出目录名；
- `fov`：视野名称；
- `sample`：样本名；
- `group`：实验分组。

如果一个 FOV 通过 ROI 拆成多个样本，沿用原代码格式：`roi1.csv_sample1,roi2.csv_sample2`，并在 `config/ROI/` 放置对应 CSV，同时把 `project.py` 中 `APPLY_ROI` 设为 `True`。
