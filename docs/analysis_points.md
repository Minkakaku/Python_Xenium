# Python/Xenium 主流程分析点

| 顺序 | 模块 | 分析点 | 主要输出 |
|---|---|---|---|
| 01 | io.py、qc.py | Xenium cell-feature matrix、边界和 metadata 读取；QC | QC 表和图 |
| 02 | qc.py | 细胞过滤 | 过滤 h5ad |
| 03 | analysis.py、cluster_gpu.py | PCA、邻居图、Leiden、UMAP；可选 RAPIDS GPU 聚类 | 聚类 h5ad 与 UMAP |
| 04 | analysis.py | marker gene | 完整 marker 表 |
| 05 | steps.py | 人工 cell-type annotation | 注释 h5ad 和空间图 |
| 06 | spatial.py | 空间自相关和空间统计 | Moran 等结果 |
| 07 | spatial.py | 同一 FOV 内细胞邻域 | 邻域组成表 |
| 08 | spatial.py | 细胞类型 proximity | 样本级距离统计 |
| 09 | spatial.py | nearest-neighbor | 最近邻类型和距离表 |

空间邻居和距离必须在每个样本/FOV 内独立计算，不能跨切片连接。各项目在自己的 `workflow/` 中保留短步骤入口；`xenium_workflow.cluster_gpu` 是独立 GPU 聚类入口。
