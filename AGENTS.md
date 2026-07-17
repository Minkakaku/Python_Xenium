# Xenium 项目代码约定

- 所有可复用函数和分析方法放在 `src/`。
- 项目参数、样本表、QC 阈值和人工注释放在 `config/`。
- `workflow/` 初始为空；每运行一个新分析点，新建 `stepXX_name.py`。
- 已运行的 workflow 不覆盖；参数变化或重跑使用新的步骤文件。
- workflow 只加载配置并调用本次方法，避免复制大段函数。
- 原始 notebook 和脚本放在 `legacy/original/`，不得直接修改。
- `cluster_gpu.py` 是独立 GPU 聚类入口，保留 RAPIDS、CuPy 和 RMM 路线。
- 不在 `src/` 中写项目编号、样本名、绝对路径、固定 cluster 注释或 ROI。
- 空间邻居必须按样本/FOV 分开计算；组间统计使用样本级汇总。
- 最终图同时保存 PNG 和 PDF，绘图输入表保存 TSV/CSV。
- 优先保留原始分析意图，代码保持简洁，不添加过度防御性封装。
