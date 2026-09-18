# Xenium 项目代码约定

- 跨项目复用函数放在 `src/python/`；公共默认放在根目录 `config/`；项目参数和人工决策放在 `projects/<项目>/config/`。
- 项目 `workflow/` 只放编号短入口；已运行脚本不覆盖，改参或重跑时新增下一编号脚本。
- 原始 notebook 和脚本放在 `legacy/original/`，不得修改；`source_catalog/` 只作为方法参考。
- `xenium_workflow.cluster_gpu` 是独立 RAPIDS/CuPy/RMM 聚类入口。
- 不在复用源码中写项目编号、样本名、绝对路径、固定 cluster 注释或 ROI。
- 空间邻居必须按样本/FOV 分开计算；组间统计使用样本级汇总。
- 中间 H5AD 写入项目 `data/processed/`；结果写入项目编号 `results/` 目录；日志写入项目 `logs/`。
- 正式图同时保存 PNG 和 PDF，绘图输入表保存 TSV/CSV。
- 正式入口只使用 Python 脚本，不使用 notebook。
- 优先保留原始分析意图，代码保持简洁，不添加过度防御性封装。
