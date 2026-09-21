# 公共配置

- `default.yaml`：读取、QC、GPU 聚类和空间分析的公共默认参数。
- `reference.yaml`：物种和参考相关默认值。

项目参数放在 `projects/<项目>/config/config.yaml`，同名顶级键覆盖公共值。样本表、ROI、QC 阈值和人工 cell-type 映射属于项目配置。
