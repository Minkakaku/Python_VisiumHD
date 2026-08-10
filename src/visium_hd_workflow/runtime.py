from pathlib import Path

import matplotlib.pyplot as plt


REPORT_STEPS = (
    "01_Data_QC",
    "02_Bin_filter",
    "03_Dimension_reduction_cluster",
    "04_Marker_gene_analysis",
    "05_Bin_annotation",
    "06_Spatial_statistics",
    "07_Cell_neighborhood",
    "08_Proximity_analysis",
    "09_Nearest_neighbor",
)


# 功能：建立模板规定的中间对象、日志和结果目录。
def prepare_project_dirs(project):
    project.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    project.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    project.LOGS_DIR.mkdir(parents=True, exist_ok=True)
    for name in REPORT_STEPS:
        (project.RESULTS_DIR / name).mkdir(parents=True, exist_ok=True)


# 功能：把同一张正式图同时保存为 PNG 和 PDF。
def save_figure(path_no_ext, fig=None, dpi=300, close=True):
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig = fig or plt.gcf()
    fig.savefig(f"{path_no_ext}.png", dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", facecolor="white", bbox_inches="tight")
    if close:
        plt.close(fig)


# 功能：为分析结果目录写入简洁中文说明，说明输入、方法和解释边界。
def write_result_readme(out_dir, title, text):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    content = f"# {title}\n\n{text.strip()}\n"
    (out_dir / "README_results_CN.md").write_text(content, encoding="utf-8")


# 功能：为 CON、EXP 及其他组生成稳定配色。
def group_palette(groups, configured=None):
    configured = configured or {}
    fallback = plt.get_cmap("tab10")
    return {
        group: configured.get(group, fallback(i % 10))
        for i, group in enumerate(groups)
    }


# 功能：计算所有分组共享的 UMAP 坐标范围，避免分面图视觉尺度不一致。
def embedding_limits(adata, key="X_umap"):
    coords = adata.obsm[key]
    return (coords[:, 0].min(), coords[:, 0].max()), (coords[:, 1].min(), coords[:, 1].max())


# 功能：统一空间散点图；Visium 图像坐标的 y 轴向下，因此绘图时反转 y 轴。
def spatial_scatter(adata, color, ax=None, coord_key="spatial", size=2, title=None, palette=None):
    import numpy as np
    import pandas as pd

    ax = ax or plt.subplots(figsize=(7, 7))[1]
    coords = np.asarray(adata.obsm[coord_key])
    values = adata.obs[color]
    if pd.api.types.is_numeric_dtype(values):
        points = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=size, cmap="viridis", linewidths=0)
        plt.colorbar(points, ax=ax, fraction=0.04, pad=0.02)
    else:
        categories = values.astype("category").cat.categories
        colors = palette or dict(zip(categories, plt.get_cmap("tab20").colors))
        for category in categories:
            keep = values.astype(str).to_numpy() == str(category)
            ax.scatter(
                coords[keep, 0],
                coords[keep, 1],
                s=size,
                color=colors.get(category, colors.get(str(category))),
                label=str(category),
                linewidths=0,
            )
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=3, frameon=False)
    ax.set_title(title or color)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    return ax
