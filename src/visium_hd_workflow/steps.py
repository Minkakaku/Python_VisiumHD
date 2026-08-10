from .analysis import annotate_bins, cluster_plots, marker_genes, spatial_gene_plots
from .io import load_visium_hd_data
from .preprocessing import cluster_cpu
from .qc import data_qc, filter_bins, write_qc_thresholds
from .runtime import prepare_project_dirs
from .spatial import (
    cell_neighborhoods,
    nearest_neighbor_analysis,
    proximity_analysis,
    spatial_statistics,
)


# 功能：读取样本表中的 Space Ranger 输出，写出带唯一 barcode 和 µm 坐标的原始对象。
def step_load(project):
    prepare_project_dirs(project)
    adata = load_visium_hd_data(project.SAMPLE_INFO, project.BIN_SIZE)
    adata.write(project.ADATA_RAW)
    return adata


# 功能：计算 QC、写出建议阈值；人工确认阈值后再进入过滤步骤。
def step_qc(project):
    import scanpy as sc

    prepare_project_dirs(project)
    adata = sc.read_h5ad(project.ADATA_RAW)
    adata = data_qc(
        adata,
        project.RESULTS_DIR / "01_Data_QC",
        project.SPECIES_MT_PREFIX,
        plot_cells_per_sample=project.QC_PLOT_CELLS_PER_SAMPLE,
        seed=project.RANDOM_SEED,
        lower_candidate_quantiles=project.QC_LOWER_CANDIDATE_QUANTILES,
        candidate_quantile=project.QC_REVIEW_QUANTILE,
        min_genes=project.QC_MIN_GENES,
        max_gene_quantile=project.QC_MAX_GENE_QUANTILE,
        max_pct_counts_mt=project.QC_MAX_MT_PERCENT,
        min_complexity=project.QC_MIN_COMPLEXITY,
    )
    adata.write(project.ADATA_QC)
    write_qc_thresholds(
        adata,
        project.QC_THRESHOLDS,
        project.QC_MIN_GENES,
        project.QC_MIN_GENE_QUANTILE,
        project.QC_MAX_GENE_QUANTILE,
        project.QC_MAX_MT_PERCENT,
        project.QC_MIN_COMPLEXITY,
    )
    return adata


# 功能：使用人工确认后的逐样本阈值过滤，并写出中间 H5AD。
def step_filter(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_QC)
    adata = filter_bins(
        adata,
        project.QC_THRESHOLDS,
        project.RESULTS_DIR / "02_Bin_filter",
        plot_cells_per_sample=project.QC_PLOT_CELLS_PER_SAMPLE,
        seed=project.RANDOM_SEED,
    )
    adata.write(project.ADATA_FILTERED)
    return adata


# 功能：在 CPU 环境执行聚类；大对象可改用独立 src/cluster_gpu.py。
def step_cluster_cpu(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_FILTERED)
    adata = cluster_cpu(
        adata,
        project.RESULTS_DIR / "03_Dimension_reduction_cluster",
        n_top_genes=project.N_TOP_GENES,
        n_neighbors=project.N_NEIGHBORS,
        n_pcs=project.N_PCS,
        max_pcs=project.MAX_PCS,
        resolutions=project.CLUSTER_RESOLUTIONS,
        use_harmony=project.USE_HARMONY,
        seed=project.RANDOM_SEED,
    )
    adata.write(project.ADATA_CLUSTERED)
    return adata


# 功能：绘制聚类结果并计算 marker；marker 默认使用 Wilcoxon。
def step_cluster_plots_and_markers(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    cluster_plots(adata, project.RESOLUTION_KEY, project.RESULTS_DIR / "03_Dimension_reduction_cluster")
    marker_genes(adata, project.RESOLUTION_KEY, project.RESULTS_DIR / "04_Marker_gene_analysis")
    return adata


# 功能：读取人工映射完成 bin 注释，并按配置绘制空间基因表达。
def step_annotation(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_CLUSTERED)
    adata = annotate_bins(
        adata,
        project.RESOLUTION_KEY,
        project.ANNOTATION_FILE,
        project.RESULTS_DIR / "05_Bin_annotation",
    )
    spatial_gene_plots(adata, project.SPATIAL_GENES, project.RESULTS_DIR / "05_Bin_annotation")
    adata.write(project.ADATA_ANNOTATED)
    return adata


# 功能：逐切片运行中心性和可选 Moran's I。
def step_spatial_statistics(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return spatial_statistics(
        adata,
        project.RESULTS_DIR / "06_Spatial_statistics",
        moran_genes=project.MORAN_GENES,
    )


# 功能：运行一个或多个 CN 数量，并将标签写入独立中间对象。
def step_cell_neighborhood(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    for n_clusters in project.CN_CLUSTERS:
        adata = cell_neighborhoods(
            adata,
            project.RESULTS_DIR / "07_Cell_neighborhood",
            k=project.CN_NEIGHBORS,
            n_clusters=n_clusters,
            group_colors=project.GROUP_COLORS,
        )
    adata.write(project.PROCESSED_DIR / "05_adata_with_neighborhoods.h5ad")
    return adata


# 功能：按 µm 距离运行半径和环带邻近分析。
def step_proximity(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return proximity_analysis(
        adata,
        project.RESULTS_DIR / "08_Proximity_analysis",
        radius=project.PROXIMITY_RADIUS_UM,
        distance_bins=project.PROXIMITY_BINS_UM,
        group_colors=project.GROUP_COLORS,
    )


# 功能：逐切片运行不同注释之间的最近距离分析。
def step_nearest_neighbor(project):
    import scanpy as sc

    adata = sc.read_h5ad(project.ADATA_ANNOTATED)
    return nearest_neighbor_analysis(
        adata,
        project.RESULTS_DIR / "09_Nearest_neighbor",
        project.GROUP_COLORS,
    )
