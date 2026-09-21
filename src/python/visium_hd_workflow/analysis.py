from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from .runtime import (
    NPG_SEQUENTIAL,
    categorical_palette,
    embedding_limits,
    save_figure,
    set_adata_palette,
    spatial_scatter,
    style_dotplot,
    write_result_readme,
)


# 功能：固定所有分组 UMAP 的坐标范围和纵横比。
def _fix_embedding_axes(axes, xlim, ylim):
    for ax in np.atleast_1d(axes).ravel():
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")


# 功能：绘制全对象、分组和逐样本空间 cluster 图，并输出样本级比例。
def cluster_plots(adata, cluster_key, out_dir):
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xlim, ylim = embedding_limits(adata)
    set_adata_palette(adata, cluster_key)

    sc.pl.umap(
        adata,
        color=cluster_key,
        title="Visium HD cluster",
        legend_loc="on data",
        legend_fontoutline=2,
        show=False,
    )
    _fix_embedding_axes(plt.gca(), xlim, ylim)
    save_figure(out_dir / "umap_cluster")

    groups = adata.obs["group"].drop_duplicates().tolist()
    if len(groups) > 1:
        fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 5), squeeze=False)
        for ax, group in zip(axes.ravel(), groups):
            sc.pl.umap(
                adata[adata.obs["group"] == group],
                color=cluster_key,
                title=str(group),
                legend_loc="on data",
                legend_fontoutline=2,
                ax=ax,
                show=False,
            )
        _fix_embedding_axes(axes, xlim, ylim)
        fig.tight_layout()
        save_figure(out_dir / "umap_cluster_by_group", fig)

    cluster_levels = adata.obs[cluster_key].astype("category").cat.categories
    palette = categorical_palette(cluster_levels)
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        fig, ax = plt.subplots(figsize=(7, 7))
        spatial_scatter(tmp, cluster_key, ax=ax, size=2, title=str(sample), palette=palette)
        save_figure(out_dir / f"spatial_{sample}_cluster", fig)

    counts = pd.crosstab(adata.obs[cluster_key], adata.obs["sample"])
    proportions = pd.crosstab(adata.obs["sample"], adata.obs[cluster_key], normalize="index")
    counts.to_csv(out_dir / "cluster_counts_by_sample.tsv", sep="\t")
    proportions.to_csv(out_dir / "cluster_proportions_by_sample.tsv", sep="\t")
    pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"]).to_csv(
        out_dir / "umap_coordinates.tsv", sep="\t"
    )
    ax = proportions.plot(
        kind="bar",
        stacked=True,
        figsize=(max(6, len(proportions) * 0.8), 4.5),
        width=0.9,
        color=[palette.get(x) for x in proportions.columns],
    )
    ax.set(xlabel="", ylabel="Bin proportion", title="Cluster proportions by sample")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(False)
    save_figure(out_dir / "cluster_proportions_by_sample", ax.figure)
    write_result_readme(
        out_dir,
        "03 Dimension reduction and clustering",
        "在过滤后对象上完成标准化、HVG、PCA 自动拐点、可选 Harmony、邻接图、UMAP 与 Leiden。空间图逐切片绘制；比例以样本汇总。",
    )
    return adata


# 功能：用 Wilcoxon 计算每个 cluster 的正向 marker，并保存 top marker 图与完整表。
def marker_genes(adata, cluster_key, out_dir, method="wilcoxon", top_n=5):
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        use_raw=adata.raw is not None,
        method=method,
        pts=True,
    )
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=cluster_key,
        n_genes=top_n,
        standard_scale="var",
        cmap=LinearSegmentedColormap.from_list("npg_dotplot", NPG_SEQUENTIAL),
        smallest_dot=18,
        show=False,
    )
    style_dotplot()
    save_figure(out_dir / f"dotplot_top{top_n}_markers")

    categories = adata.obs[cluster_key].astype("category").cat.categories
    tables = []
    for cluster in categories:
        table = sc.get.rank_genes_groups_df(adata, group=cluster)
        table = table[(table["logfoldchanges"] > 0) & (table["pvals_adj"] < 0.05)].copy()
        table["cluster"] = str(cluster)
        table.to_csv(out_dir / f"cluster_{cluster}_positive_markers.tsv", sep="\t", index=False)
        tables.append(table)
    all_markers = pd.concat(tables, ignore_index=True)
    all_markers.to_csv(out_dir / "all_cluster_positive_markers.tsv", sep="\t", index=False)
    write_result_readme(
        out_dir,
        "04 Marker gene analysis",
        "使用 log-normalized 全基因 raw 表达和 Wilcoxon 检验计算 cluster marker。这里是 bin/分割单元层面的描述性 marker，不代表样本间差异表达。",
    )
    return all_markers


# 功能：根据显式 cluster-to-label 表进行人工注释，并输出 UMAP、空间图和样本比例。
def annotate_bins(adata, cluster_key, annotation_file, out_dir):
    import scanpy as sc

    mapping = pd.read_csv(annotation_file, sep="\t", dtype=str).dropna()
    annotation = dict(zip(mapping["cluster"], mapping["celltype"]))
    clusters = adata.obs[cluster_key].astype(str)
    missing = sorted(set(clusters) - set(annotation))
    if missing:
        raise ValueError(f"注释表缺少 cluster: {', '.join(missing)}")
    adata.obs["celltype"] = pd.Categorical(clusters.map(annotation))

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    xlim, ylim = embedding_limits(adata)
    set_adata_palette(adata, "celltype")
    sc.pl.umap(adata, color="celltype", title="Bin annotation", legend_loc="right margin", show=False)
    _fix_embedding_axes(plt.gca(), xlim, ylim)
    save_figure(out_dir / "umap_annotation")

    palette = categorical_palette(adata.obs["celltype"].cat.categories)
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        fig, ax = plt.subplots(figsize=(7, 7))
        spatial_scatter(tmp, "celltype", ax=ax, size=2, title=str(sample), palette=palette)
        save_figure(out_dir / f"spatial_{sample}_annotation", fig)

    proportions = pd.crosstab(adata.obs["sample"], adata.obs["celltype"], normalize="index")
    proportions.to_csv(out_dir / "annotation_proportions_by_sample.tsv", sep="\t")
    adata.obs.to_csv(out_dir / "all_bin_annotations.tsv", sep="\t")
    ax = proportions.plot(
        kind="bar",
        stacked=True,
        figsize=(max(6, len(proportions) * 0.9), 4.5),
        width=0.9,
        color=[palette.get(x) for x in proportions.columns],
    )
    ax.set(xlabel="", ylabel="Bin proportion", title="Annotation proportions by sample")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(False)
    save_figure(out_dir / "annotation_proportions_by_sample", ax.figure)
    write_result_readme(
        out_dir,
        "05 Bin annotation",
        "依据 marker 和组织知识用 config/cluster_annotation.tsv 显式映射。8 µm bin 标签表示主要表达特征，不应自动解释成单个纯细胞。",
    )
    return adata


# 功能：逐样本绘制配置指定基因的空间表达，避免不同切片坐标互相混合。
def spatial_gene_plots(adata, genes, out_dir):
    out_dir = Path(out_dir)
    available = [gene for gene in genes if gene in adata.var_names or (adata.raw is not None and gene in adata.raw.var_names)]
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        for gene in available:
            values = tmp.raw[:, gene].X if tmp.raw is not None and gene in tmp.raw.var_names else tmp[:, gene].X
            values = values.toarray().ravel() if hasattr(values, "toarray") else np.asarray(values).ravel()
            tmp.obs[f"__plot_{gene}"] = values
            fig, ax = plt.subplots(figsize=(7, 7))
            spatial_scatter(tmp, f"__plot_{gene}", ax=ax, size=2, title=f"{sample}: {gene}")
            save_figure(out_dir / f"spatial_{sample}_{gene}", fig)
