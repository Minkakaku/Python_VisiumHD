# 原始代码审查与迁移说明

## 保留的分析意图

- Space Ranger v4 count 参数；
- Visium HD 8 µm bins 与 v4 segmentation 读取；
- 逐样本 QC、过滤、UMAP、Leiden、marker 和人工注释；
- 独立 RAPIDS GPU 聚类；
- 空间中心性、cellular neighborhood、半径 proximity 和距离环带。

## 已修正的问题

1. 删除项目编号、固定 Linux 路径、固定肝脏 cluster 字典和客户样本名。
2. 不再把 `tissue_positions.parquet` 写回 Space Ranger 目录伪装成 CSV。
3. barcode 改为 `sample:barcode`，避免多样本拼接后身份不稳定。
4. 分割 polygon 质心按 instance ID 对齐，避免 GeoDataFrame 行顺序与表达表不一致。
5. `log10GenesPerUMI` 对 0/1 UMI 安全计算。
6. 删除错误的线粒体比例下限，只保留可配置上限。
7. marker 默认由 t-test 改为 Wilcoxon，并从筛选 HVG 前保存完整 raw。
8. CN、空间图和距离全部按样本分开；原代码直接在合并坐标上建图会产生跨样本假邻居。
9. CN 使用归一化邻域组成，并以 `log2(CN proportion / tissue proportion)` 表示富集。
10. proximity 组间检验由逐 bin 伪重复改为样本级汇总。
11. 所有正式图保存 PNG/PDF，正式表改为真实 TSV 扩展名。
12. Enrichr 属于联网可选分析，不再作为主流程成功与否的默认条件。

## 暂未自动化的判断

- QC 阈值必须看组织与实验质量人工确认；
- PCA 自动拐点是参考，建议比较相邻 PC 数的聚类稳定性；
- cluster-to-celltype 映射必须由 marker、组织学与参考知识共同判断；
- 8 µm bin 是否适合写成“cell type”取决于样本和验证，报告中更稳妥的表述是 bin annotation 或 dominant program。
