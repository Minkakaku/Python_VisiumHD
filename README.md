# Python-visium-hd-workflow-template

这是从原始 Visium HD notebook、`CN.py` 和 Space Ranger v4 命令整理出的可复用 Python 主流程。正式入口全部使用 `.py`；原始文件完整保存在 `legacy/original/`。

默认分析单位是 Space Ranger `square_008um` 的 8 µm bin，也支持 Space Ranger v4 分割结果。两者不能混为一谈：bin 注释表示主要表达特征，只有明确使用 `input_mode=segmented` 时才是分割单元。

## 正式分析点

| 顺序 | 结果目录 | 方法入口 | 作用 |
|---|---|---|---|
| 01 | `results/01_Data_QC/` | `step_qc()` | 全样本/逐样本四联小提琴、counts-feature 散点、分布/密度、候选保留率和阈值建议 |
| 02 | `results/02_Bin_filter/` | `step_filter()` | 按人工确认的逐样本阈值过滤 |
| 03 | `results/03_Dimension_reduction_cluster/` | `step_cluster_cpu()` 或 `cluster_gpu.py` | HVG、PCA 自动拐点、可选 Harmony、UMAP、Leiden |
| 04 | `results/04_Marker_gene_analysis/` | `marker_genes()` | 全基因 log-normalized raw 上的 Wilcoxon marker |
| 05 | `results/05_Bin_annotation/` | `step_annotation()` | cluster 显式映射、空间图和样本比例 |
| 06 | `results/06_Spatial_statistics/` | `step_spatial_statistics()` | 逐切片空间图中心性和可选 Moran's I |
| 07 | `results/07_Cell_neighborhood/` | `step_cell_neighborhood()` | 逐切片 kNN 组成、CN 和样本比例 |
| 08 | `results/08_Proximity_analysis/` | `step_proximity()` | 半径邻近与距离环带 |
| 09 | `results/09_Nearest_neighbor/` | `step_nearest_neighbor()` | 不同注释之间的最近距离 |

## 安装

CPU 环境：

```bash
conda env create -f environment.yml
conda activate visium-hd-workflow
pip install -e .
```

GPU 聚类环境：

```bash
conda env create -f environment-gpu.yml
conda activate visium-hd-gpu
pip install -e .
```

GPU 环境要求 NVIDIA 驱动与匹配的 CUDA。`cluster_gpu.py` 默认使用第 0 张 GPU，可用 `--device` 修改。

## 新项目开始

1. 复制并修改 `config/sample_info.example.tsv`。
2. 修改 `config/project.py` 的分析参数。
3. 在 `workflow/` 新建 `step01_load.py`，只加载配置并调用 `step_load(project)`。
4. 每运行一个新分析点增加一个新的 `stepXX_name.py`；不要覆盖已经运行过的步骤。
5. QC 后先检查 `config/qc_thresholds.tsv`，再运行过滤。

模板中的 `workflow/` 只有 `.gitkeep`，不会替你预先决定实际项目运行了哪些步骤。完整示例见 [新项目指南](docs/new_project.md)。

## GPU 聚类示例

```bash
python src/cluster_gpu.py \
  --input data/processed/02_adata_filtered.h5ad \
  --output data/processed/03_adata_clustered.h5ad \
  --result-dir results/03_Dimension_reduction_cluster \
  --resolution 0.4 0.5 0.6 0.7
```

不传 `--n-pcs` 时，程序使用 `kneed.KneeLocator` 自动识别 PCA 拐点。该结果是参考值；复杂组织和稀有群仍应比较邻近 PC 数下的稳定性。

## 关键修正

- 多样本 barcode 改为 `sample:barcode`，不再重置为不稳定的数字索引。
- 坐标表和分割 polygon 统一按 ID 显式对齐，不依赖文件当前行顺序。
- 所有空间图、CN、proximity 和 nearest-neighbor 均逐样本计算。
- 原代码错误的 `pct_counts_mt > 下限` 已删除，只保留可配置上限。
- marker 默认由 t-test 改为 Scanpy 常用的 Wilcoxon，并确保 `raw` 保存完整基因。
- 组间统计使用样本级汇总，不把数十万 bins 当作独立生物学重复。
- QC 对每个样本最多抽取 `QC_PLOT_CELLS_PER_SAMPLE` 个单位绘图，完整逐单位指标仍单独保存；分割对象若含 `bin_count` 会自动增加对应分布和密度面板。
- QC 红色虚线是人工复核候选值，不会自动替代 `config/qc_thresholds.tsv`；过滤后会生成同版式 QC 图。
- 正式图均输出 PNG/PDF，绘图输入和统计结果保存 TSV。

当前整理只创建本地模板，没有提交或上传 GitHub。
