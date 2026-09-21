# 方法网页与论文依据

- 10x Xenium Ranger Analysis：https://www.10xgenomics.com/support/software/xenium-ranger/latest/analysis
- Xenium data overview：https://www.10xgenomics.com/support/in-situ-gene-expression/documentation
- Scanpy clustering：https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
- Squidpy spatial graph/statistics：https://squidpy.readthedocs.io/
- Squidpy 论文：PMID 35102346

代码决策：

1. Xenium 表示分割细胞，不与 Visium HD bins 混用。
2. cluster_gpu.py 仅负责 RAPIDS/CuPy/RMM 加速聚类，CPU 主流程仍可独立运行。
3. 空间邻居按样本/FOV 分开建立，避免跨切片假邻居。
4. 组间空间统计先汇总到样本，再比较 CON/EXP。
5. 注释依赖 marker 与空间位置共同复核，聚类编号不是最终细胞类型。