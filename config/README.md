# 公共配置

- `default.yaml`：QC、聚类和空间分析的公共默认参数。
- `reference.yaml`：物种和参考相关默认值。

项目参数放在 `projects/<项目>/config/config.yaml`，同名嵌套映射逐键递归合并，列表和标量由项目值整体替换。样本表、QC 阈值、cluster 注释和 Space Ranger 参数也属于项目配置，不放在这里。
