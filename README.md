# Python-visium-hd-workflow-template

这是从原始 Visium HD notebook、`CN.py` 和 Space Ranger v4 命令整理出的分层 Python 主流程。公共方法与项目运行记录分开管理，正式入口全部使用 `.py`；原始文件完整保存在 `legacy/original/`。

默认分析单位是 Space Ranger `square_008um` 的 8 µm bin，也支持 Space Ranger v4 分割结果。两者不能混为一谈：bin 注释表示主要表达特征，只有明确使用 `input_mode=segmented` 时才是分割单元。

## 项目结构

```text
├── src/python/visium_hd_workflow/   # 跨项目复用函数与分析步骤
├── src/python/visium_hd_workflow/   # CPU 主流程、GPU 聚类和 Space Ranger 入口
├── config/                          # 公共默认与参考配置
├── projects/project1/               # 项目配置、workflow、数据、结果和日志
├── envs/                            # 公共 CPU/GPU 环境
├── tests/
└── legacy/original/                 # 未修改的原始代码
```

`config/default.yaml` 和 `config/reference.yaml` 提供公共值，`projects/project1/config/config.yaml` 按同名键覆盖。样本表、QC 阈值和人工注释放在项目自己的 `config/`。

## 正式分析点

| 顺序 | 项目结果目录 | 方法入口 | 作用 |
|---|---|---|---|
| 01 | `results/01_Data_QC/` | `step_qc()` | 全样本/逐样本 QC、候选保留率和阈值建议 |
| 02 | `results/02_Bin_filter/` | `step_filter()` | 按人工确认的逐样本阈值过滤 |
| 03 | `results/03_Dimension_reduction_cluster/` | `step_cluster_cpu()` 或 GPU 模块 | HVG、PCA 自动拐点、可选 Harmony、UMAP、Leiden |
| 04 | `results/04_Marker_gene_analysis/` | `marker_genes()` | 全基因 log-normalized raw 上的 Wilcoxon marker |
| 05 | `results/05_Bin_annotation/` | `step_annotation()` | cluster 显式映射、空间图和样本比例 |
| 06 | `results/06_Spatial_statistics/` | `step_spatial_statistics()` | 逐切片空间图中心性和可选 Moran's I |
| 07 | `results/07_Cell_neighborhood/` | `step_cell_neighborhood()` | 逐切片 kNN 组成、CN、样本比例 Welch 检验和组别 alluvial 图 |
| 08 | `results/08_Proximity_analysis/` | `step_proximity()` | 半径邻近与距离环带 |
| 09 | `results/09_Nearest_neighbor/` | `step_nearest_neighbor()` | 不同注释之间的最近距离 |

表中的 `results/` 均位于具体项目下，例如 `projects/project1/results/`。

## 安装

```bash
conda env create -f envs/environment.yaml
conda activate visium-hd-workflow
pip install -e .
```

GPU 环境使用 `envs/environment-gpu.yaml`，要求 NVIDIA 驱动与匹配的 CUDA。

Space Ranger 可先检查生成的命令：

```bash
python -m visium_hd_workflow.run_spaceranger --config projects/project1/config/spaceranger.yaml --dry-run
```

## 新项目开始

1. 复制 `projects/project1/` 为新的项目目录。
2. 修改项目 `config/samples.tsv`、`config/config.yaml` 和可选的 `config/spaceranger.yaml`。
3. 从 `workflow/00_setup.py` 开始逐步运行短入口；已运行脚本不要覆盖，参数变化时新增下一编号脚本。
4. QC 后人工确认 `config/qc_thresholds.tsv`，marker 后人工确认 `config/cluster_annotation.tsv`。

完整说明见 [新项目指南](docs/new_project.md) 和 [project1 README](projects/project1/README.md)。

## GPU 聚类示例

```bash
python -m visium_hd_workflow.cluster_gpu \
  --input projects/project1/data/processed/02_adata_filtered.h5ad \
  --output projects/project1/data/processed/03_adata_clustered.h5ad \
  --result-dir projects/project1/results/03_Dimension_reduction_cluster \
  --resolution 0.4 0.5 0.6 0.7
```

不传 `--n-pcs` 时使用 `kneed.KneeLocator` 自动识别 PCA 拐点。该结果是参考值；复杂组织和稀有群仍应比较邻近 PC 数下的稳定性。

## 关键修正

- 多样本 barcode 使用 `sample:barcode`，空间图、CN、proximity 和 nearest-neighbor 均逐样本计算。
- marker 默认使用 Wilcoxon，并确保 `raw` 保存完整基因。
- 组间统计先汇总到样本层面，不把 bins 当作独立生物学重复。
- CN 每个样本内先计算比例；双组箱线图使用 Welch t-test，组别 alluvial 图展示组内全部单位合并后的 100% 组成。
- 正式图输出同名 PNG/PDF，绘图输入和统计结果保存 TSV。
