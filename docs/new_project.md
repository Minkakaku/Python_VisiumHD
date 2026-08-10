# 新项目使用指南

## 1. 配置

复制示例配置，填写真实样本目录和分组。`spaceranger_dir` 可以是样本目录，也可以直接是 `outs/`。

```bash
cp config/sample_info.example.tsv config/sample_info.tsv
cp config/cluster_annotation.example.tsv config/cluster_annotation.tsv
```

## 2. 逐步创建 workflow

模板不预置已运行步骤。项目开始时建立：

```python
# workflow/step01_load.py
from config import project
from visium_hd_workflow.steps import step_load

step_load(project)
```

在项目根目录运行：

```bash
python -m workflow.step01_load
```

后续依次新建并运行：

```text
step02_qc.py
step03_filter.py
step04_cluster_cpu.py      # 或独立运行 src/cluster_gpu.py
step05_cluster_markers.py
step06_annotation.py
step07_spatial_statistics.py
step08_cell_neighborhood.py
step09_proximity.py
step10_nearest_neighbor.py
```

每个步骤文件只导入 `config.project` 和对应 `step_*()`，不复制 `src/` 中的大段函数。参数变化或重跑时建立新的步骤文件，保留旧运行记录。

## 3. QC 人工检查点

`step02_qc.py` 会生成 `config/qc_thresholds.tsv`。结合每个样本的小提琴图、散点图和组织覆盖检查阈值后，再运行过滤。稀疏 8 µm bins 与分割细胞不应共用同一组固定阈值。

## 4. 注释人工检查点

聚类和 marker 完成后编辑 `config/cluster_annotation.tsv`。正式代码要求每个 cluster 都有映射，避免未注释 cluster 被静默改成缺失值。

## 5. 结果解释

- UMAP 是表达相似性嵌入，不代表组织距离。
- 8 µm bin 可能含多个细胞的 RNA，单标签只是主要特征。
- CN 和 proximity 描述几何组织关系，不证明配体-受体通讯。
- CON/EXP 检验使用样本汇总；只有一个样本/组时仅做描述，不报告可靠组间推断。
