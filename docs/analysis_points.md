# Python Visium HD 分析点说明

## 主流程边界

本仓库是多分析点主流程，输入为 Space Ranger Visium HD 输出，输出一个逐步更新的 H5AD。CellChat、CellPhoneDB、轨迹等接收 H5AD/RDS 的独立方法不并入这里。

| 分析点 | 输入 | 主要方法 | 解释边界 |
|---|---|---|---|
| Space Ranger | FASTQ、组织图像、参考 | `spaceranger count` | 属于上游定量，不在 Python 中重实现 |
| 数据读取 | `square_008um` 或 v4 segmentation | 10x H5 + parquet；SpatialData segmentation | 8 µm bin 与分割细胞必须在 metadata 中区分 |
| QC/过滤 | counts | nGene、UMI、MT%、复杂度、逐样本阈值 | 不存在适合所有组织的固定阈值 |
| 聚类 | filtered H5AD | normalize、HVG、PCA knee、Harmony、Leiden | 聚类不是细胞类型注释 |
| marker | clustered H5AD | Wilcoxon one-vs-rest | 不是样本间差异表达 |
| annotation | marker + 映射表 | 人工 cluster-to-label | bin 标签代表主要表达特征 |
| 空间统计 | annotated H5AD | kNN centrality、Moran's I | Moran's I 不是 CON/EXP 差异 |
| CN | annotated H5AD | 局部标签组成 + MiniBatchKMeans | CN 是邻域状态，不是新细胞类型 |
| proximity | annotated H5AD | 半径和距离环带 | 接近不等于分子通讯 |
| nearest neighbor | annotated H5AD | 来源到目标标签最近距离 | 先汇总样本再比较组别 |

## 原始文件与正式入口对应

- `run_spacerangerv4.sh` → `src/run_spaceranger.py` + `config/spaceranger.example.yaml`
- `HD_V4.ipynb` → `io.py`、`qc.py`、`preprocessing.py`、`analysis.py`
- `HDV4_cluster_gpu.ipynb` → `src/cluster_gpu.py`
- `visium_hd_report.ipynb` → `steps.py` 和 01–09 结果目录
- `CN.py` → `spatial.py` 中按样本隔离的 CN 方法

## 为什么没有正式 notebook

notebook 包含运行状态、绝对路径和项目注释，难以复用和审计。它们原样保存在 `legacy/original/`；正式代码只提供函数与 `.py` 入口，实际项目在空的 `workflow/` 中追加步骤脚本。
