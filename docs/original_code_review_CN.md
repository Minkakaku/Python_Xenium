# 原始 Xenium 代码审查记录

## 原始文件

原文件未直接修改，完整保存在 `legacy/original/`：

- `xenium5k_report.ipynb`；
- `cluster_gpu.ipynb`；
- `cluster_gpu.py`。

它们与 `E:/NextDATA/backup/VH2604022693/` 中对应文件的 SHA-256 一致。

## 保留内容

- `spatialdata_io.xenium()` 读取 Xenium 数据；
- 四列表型 `raw / fov / sample / group`；
- 一个 FOV 通过 ROI CSV 拆分多个样本的规则；
- `data_qc()`、`filter_parameter_record()`、`cell_filter()`；
- RAPIDS + CuPy + RMM 的 GPU 聚类路线；
- `cluster_plot()`、`marker_gene()`、`celltype_anno()`；
- centrality、cellular neighborhood、proximity 和 nearest-neighbor 分析框架；
- 原始 QC 分位数、复杂度阈值、邻域半径和 CN 数量作为默认配置。

## 已修复的明确问题

1. 原 GPU 脚本读取固定 `adata_raw.h5ad`，主 notebook 实际在 QC 后生成 `adata_filtered.h5ad`。现改为必填 `--input/--output`。
2. 原 GPU 脚本只生成 `res_0.6/0.8/1`，后续却读取 `res_0.3`。现默认同时生成 `0.3/0.6/0.8/1.0`。
3. 原 GPU 脚本末尾是固定 `P1_CRC`、`hires` 和裁剪坐标的 Visium 绘图，不属于 Xenium GPU 聚类，已从新脚本删除。
4. 原目录创建依赖 `output/` 不存在；若 `output/` 已存在而子目录缺失会失败。现所有目录始终按需创建。
5. 原代码调用 `load_xenium_data(..., sample_info='sample_info.txt')`，忽略前面定义的变量。现所有路径来自配置。
6. 原 `.str.contains('.csv_')` 把点号当正则通配符。现使用字面匹配。
7. ROI 合并后直接替换 `AnnData.obs`，可能发生行顺序和索引错位。现以 `fov:cell_id` 映射并保持 AnnData 索引。
8. 原空间邻域在合并后的全部样本坐标上建图。不同 FOV 坐标可重叠，会产生跨样本假邻居。现每个样本独立建图。
9. 原 proximity 和 nearest-neighbor 组间检验以细胞为重复，形成统计伪重复。现先汇总到样本，再做 Mann–Whitney 检验，并保留 `ns`。
10. 原环形距离统计只记录出现过的邻居类型，缺失的零计数被遗漏。现零计数也进入样本级中位数。
11. 原 nearest-neighbor 图的 Y 轴写成 `median cell counts`，实际是距离，并固定截断到 15 µm。现改为样本级中位最近邻距离且使用动态范围。
12. 原报告生成写死服务器模板路径并使用 Linux `cp -r`。这部分属于外部报告系统，没有放入通用分析代码。

## 没有自动替用户决定的内容

- QC 阈值是否适合具体组织；
- 最终使用哪个 resolution；
- cluster 到 cell type 的人工注释；
- 50 µm 是否为适合当前组织的生物学距离；
- CN 聚类数 10/15 是否为最终选择；
- GO 背景库和物种是否符合当前项目。

这些参数均集中在项目 `config/config.yaml` 或人工配置表中。
