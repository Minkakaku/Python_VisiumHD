from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, to_hex


NPG_COLORS = (
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
)
NPG_SEQUENTIAL = ("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019")
NPG_DIVERGING = ("#3C5488", "#F7F7F7", "#E64B35")

# 未指定颜色的折线与连续图也使用 NPG，显式绘图参数仍可覆盖。
for name, colors in (("npg_sequential", NPG_SEQUENTIAL), ("npg_diverging", NPG_DIVERGING)):
    if name not in colormaps:
        colormaps.register(LinearSegmentedColormap.from_list(name, colors))
plt.rcParams.update({"axes.prop_cycle": plt.cycler(color=NPG_COLORS), "image.cmap": "npg_sequential"})


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
def save_figure(path_no_ext, fig=None, dpi=600, close=True):
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


# 功能：按固定 NPG 顺序生成稳定分类配色；超过 10 类时从 NPG 锚点扩展。
def categorical_palette(labels, configured=None):
    labels = list(labels)
    configured = configured or {}
    colors = list(NPG_COLORS)
    if len(labels) > len(colors):
        cmap = LinearSegmentedColormap.from_list("npg_extended", NPG_COLORS)
        colors.extend(to_hex(cmap(x)) for x in np.linspace(0, 1, len(labels) - len(colors) + 2)[1:-1])
    return {
        label: configured.get(label, colors[i])
        for i, label in enumerate(labels)
    }


# 功能：为 CON、EXP 及其他组生成稳定 NPG 配色。
def group_palette(groups, configured=None):
    return categorical_palette(groups, configured)


# 功能：把 NPG 颜色写入 AnnData，供 Scanpy 图和后续切片复用。
def set_adata_palette(adata, key, configured=None):
    if not isinstance(adata.obs[key].dtype, pd.CategoricalDtype):
        adata.obs[key] = adata.obs[key].astype("category")
    labels = adata.obs[key].cat.categories
    palette = categorical_palette(labels, configured)
    adata.uns[f"{key}_colors"] = [palette[label] for label in labels]
    return palette


# 功能：把 Scanpy DotPlot 调整为紧凑白底、纵向 marker 标签和水平辅助线。
def style_dotplot(fig=None):
    fig = fig or plt.gcf()
    for ax in fig.axes:
        if not ax.get_xticklabels() or not ax.get_yticklabels():
            continue
        ax.set_facecolor("white")
        ax.tick_params(axis="x", rotation=90)
        ax.set_axisbelow(True)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True, color="#E6E6E6", linewidth=0.5)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(0.6)
    return fig


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
        cmap = LinearSegmentedColormap.from_list("npg_sequential", NPG_SEQUENTIAL)
        points = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=size, cmap=cmap, linewidths=0)
        plt.colorbar(points, ax=ax, fraction=0.04, pad=0.02)
    else:
        categories = values.astype("category").cat.categories
        colors = palette or categorical_palette(categories)
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
    if title is not None:
        ax.set_title(title)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    return ax
