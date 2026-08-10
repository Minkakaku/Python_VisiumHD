# 配置文件

新项目中先复制三个示例文件：

```bash
cp config/sample_info.example.tsv config/sample_info.tsv
cp config/qc_thresholds.example.tsv config/qc_thresholds.tsv
cp config/cluster_annotation.example.tsv config/cluster_annotation.tsv
```

`sample_info.tsv` 每行是一个独立切片。`spaceranger_dir` 可指向样本目录或其 `outs/`；`input_mode` 默认用 `bins` 读取 `square_008um`。Space Ranger v4 分割结果可写 `segmented`，并按实际 SpatialData 元素填写 `table_key`、`shape_key`。 分割坐标默认从 Space Ranger `scalefactors_json.json` 读取 `microns_per_pixel`；若坐标系不同，在 `microns_per_unit` 显式填写换算值。

QC 步骤会生成新的 `config/qc_thresholds.tsv` 建议值。人工检查后再运行过滤步骤，不要把模板示例阈值直接用于所有组织。
QC 绘图由 `QC_PLOT_CELLS_PER_SAMPLE` 控制逐样本抽样上限，`QC_LOWER_CANDIDATE_QUANTILES` 控制候选低端分位数，`QC_REVIEW_QUANTILE` 决定图中重点显示哪一个候选。候选线只用于人工判断；正式过滤仍读取人工确认后的 `config/qc_thresholds.tsv`。