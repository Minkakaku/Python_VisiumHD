# 方法网页与文献依据

核对日期：2026-08-04。

1. 10x Genomics, Space Ranger output overview

   https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview

   用于确认 `outs/`、feature-barcode matrix、spatial 输出和 binned output 的职责。

2. spatialdata-io: `spatialdata_io.visium_hd()`

   https://spatialdata.scverse.org/projects/io/en/stable/generated/spatialdata_io.visium_hd.html

   当前接口支持 `bin_size`、filtered/raw counts、`load_segmentations_only`、nucleus segmentation 和多图像读取。正式模板默认直接加载稳定的 8 µm matrix + parquet；Space Ranger v4 segmentation 通过该接口加载。

3. Scanpy: `read_visium()`

   https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html

   该接口面向经典 Visium 目录。原 notebook 为兼容 HD parquet 曾把位置表临时改写成 CSV；正式模板改为直接读取 `filtered_feature_bc_matrix.h5` 和 `tissue_positions.parquet`，不修改 Space Ranger 输出。

4. Scanpy preprocessing and clustering tutorial

   https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html

   正式顺序采用 normalize/log1p、HVG、PCA、neighbors、UMAP、Leiden；marker 使用 Wilcoxon。

5. kneed / Kneedle

   https://kneed.readthedocs.io/

   Satopaa V et al. *Finding a Kneedle in a Haystack: Detecting Knee Points in System Behavior.* ICDCS Workshops, 2011.

   PCA 方差曲线使用 `curve="convex"`、`direction="decreasing"`。无明确 knee 时返回配置的保守 fallback，而不是反复调参强行制造拐点。

6. RAPIDS-singlecell

   https://rapids-singlecell.readthedocs.io/

   `cluster_gpu.py` 保留用户原有 CuPy、RMM 和 RAPIDS 路线，并在 GPU 聚类结束后写回 CPU H5AD。

7. Squidpy spatial graph API

   https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html

   Squidpy 1.8 已将通用 `spatial_neighbors()` 标记为弃用，并提供 `spatial_neighbors_knn()`、`spatial_neighbors_radius()` 等模式专用函数。正式代码使用 `spatial_neighbors_knn()`。

8. Palla G et al. *Squidpy: a scalable framework for spatial omics analysis.* Nature Methods, 2022.

   https://doi.org/10.1038/s41592-021-01358-2

   支持把空间图、邻域富集、共现、中心性和 Moran's I 作为不同分析模块解释。

9. Schürch CM et al. *Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front.* Cell, 2020.

   https://doi.org/10.1016/j.cell.2020.07.005

   CN 采用“局部窗口标签组成再聚类”的常用思想。模板将其用于注释后的 Visium HD 单元，并明确它不是新的细胞类型。
