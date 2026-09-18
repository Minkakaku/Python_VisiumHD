# 新项目使用指南

## 1. 建立项目

复制 `projects/project1/`，只改项目目录中的配置和 workflow，不复制 `src/python/`：

```bash
cp -r projects/project1 projects/my_project
```

`config/samples.tsv` 中的相对路径以项目目录为基准。`spaceranger_dir` 可指向样本目录或其 `outs/`；`input_mode` 默认用 `bins` 读取 `square_008um`，Space Ranger v4 分割结果使用 `segmented`。

## 2. 配置覆盖

读取顺序为：

1. 根目录 `config/default.yaml`；
2. 根目录 `config/reference.yaml`；
3. 项目 `config/config.yaml`。

后读取的同名顶级键覆盖前面的值。保持键扁平即可，不需要额外配置框架。

## 3. 逐步运行

在仓库根目录执行：

```bash
python projects/my_project/workflow/00_setup.py
python projects/my_project/workflow/01_load.py
python projects/my_project/workflow/02_qc.py
```

QC 后人工确认项目 `config/qc_thresholds.tsv`，再运行过滤。marker 完成后编辑 `config/cluster_annotation.tsv`，再运行注释和空间分析。每个步骤只加载项目配置并调用一个 `step_*()`；不要把 `src/python/` 中的函数复制进 workflow。

## 4. 结果解释

- UMAP 是表达相似性嵌入，不代表组织距离。
- 8 µm bin 可能含多个细胞的 RNA，单标签只是主要特征。
- CN 和 proximity 描述几何组织关系，不证明配体-受体通讯。
- CON/EXP 检验使用样本汇总；只有一个样本/组时仅做描述。
