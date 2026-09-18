# 调用独立 GPU 聚类入口；输入为过滤后对象，输出为聚类 H5AD。

from pathlib import Path

from visium_hd_workflow.cluster_gpu import run_gpu_clustering
from visium_hd_workflow.project import load_project

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
project = load_project(Path(__file__).resolve().parents[1])
run_gpu_clustering(
    project.ADATA_FILTERED,
    project.ADATA_CLUSTERED,
    project.RESULTS_DIR / "03_Dimension_reduction_cluster",
    resolutions=project.CLUSTER_RESOLUTIONS,
    n_top_genes=project.N_TOP_GENES,
    n_neighbors=project.N_NEIGHBORS,
    n_pcs=project.N_PCS,
    max_pcs=project.MAX_PCS,
    use_harmony=project.USE_HARMONY,
    device=project.GPU_DEVICE,
    seed=project.RANDOM_SEED,
)
