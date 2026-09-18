from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, ttest_ind

from .runtime import group_palette, save_figure, spatial_scatter, write_result_readme


NPG15 = (
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#D3B3E5", "#E6A176", "#5F9F8B", "#6C79A8", "#C85C74",
)


# 功能：稳定取得分类顺序，保证不同样本的邻域组成列完全一致。
def categories(values):
    return values.astype("category").cat.categories.astype(str).tolist()


# 功能：对样本级汇总执行双组检验，而不是把每个 bin 当重复。
def two_group_tests(table, value_col, key_cols, group_col="group", method="mannwhitney"):
    groups = table[group_col].drop_duplicates().tolist()
    rows = []
    if len(groups) != 2:
        return pd.DataFrame(columns=[*key_cols, "group_1", "group_2", "test_method", "p_value", "p_signif"])
    for keys, data in table.groupby(key_cols, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        x = data.loc[data[group_col] == groups[0], value_col].dropna()
        y = data.loc[data[group_col] == groups[1], value_col].dropna()
        if method == "welch_ttest":
            p_value = ttest_ind(x, y, equal_var=False).pvalue if len(x) >= 2 and len(y) >= 2 else np.nan
        elif method == "mannwhitney":
            p_value = mannwhitneyu(x, y).pvalue if len(x) and len(y) else np.nan
        else:
            raise ValueError(f"Unsupported two-group test: {method}")
        p_signif = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
        rows.append(
            dict(zip(key_cols, keys))
            | {"group_1": groups[0], "group_2": groups[1], "test_method": method,
               "p_value": p_value, "p_signif": p_signif}
        )
    return pd.DataFrame(rows)


# 功能：按 lyh CN 流程绘制样本比例箱线图和组内合并比例 alluvial 图。
def cn_group_comparison(sample_table, out_dir, n_clusters, group_colors=None):
    import seaborn as sns
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    from matplotlib.ticker import PercentFormatter

    out_dir = Path(out_dir)
    cn_order = [f"CN{x}" for x in range(n_clusters)]
    group_order = sample_table["group"].drop_duplicates().tolist()
    palette = group_palette(group_order, group_colors)
    tests = two_group_tests(sample_table, "proportion", ["CN"], method="welch_ttest")

    fig, ax = plt.subplots(figsize=(max(7, n_clusters * 0.72), 4.6))
    sns.boxplot(data=sample_table, x="CN", y="proportion", hue="group", order=cn_order,
                hue_order=group_order, palette=palette, showfliers=False, width=0.55,
                linewidth=0.75, ax=ax)
    sns.stripplot(data=sample_table, x="CN", y="proportion", hue="group", order=cn_order,
                  hue_order=group_order, palette=palette, dodge=True, legend=False,
                  size=4.2, linewidth=0.4, edgecolor="white", ax=ax)
    pad = max(0.005, sample_table["proportion"].max() * 0.06)
    top = sample_table["proportion"].max()
    for i, cn_id in enumerate(cn_order):
        label = tests.loc[tests["CN"].astype(str) == cn_id, "p_signif"]
        if len(label) and label.iloc[0]:
            y = sample_table.loc[sample_table["CN"].astype(str) == cn_id, "proportion"].max() + pad
            ax.text(i, y, label.iloc[0], ha="center", va="bottom", fontsize=10)
            top = max(top, y + pad)
    ax.set(xlabel="", ylabel="CN proportion", ylim=(0, top + pad))
    ax.tick_params(axis="x", rotation=45, length=0)
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(False)
    sns.despine(ax=ax)
    save_figure(out_dir / f"CNs_{n_clusters}_proportions_by_sample", fig, dpi=600)

    pooled = sample_table.pivot_table(index="group", columns="CN", values="n_cells", aggfunc="sum", fill_value=0)
    pooled = pooled.reindex(index=group_order, columns=cn_order, fill_value=0)
    pooled = pooled.div(pooled.sum(axis=1), axis=0)
    colors = dict(zip(cn_order, NPG15[:n_clusters]))
    x = np.arange(len(group_order))
    bottoms = np.zeros(len(group_order))
    fig, ax = plt.subplots(figsize=(4.8, 5.2))
    for cn_id in cn_order:
        values = pooled[cn_id].to_numpy()
        tops = bottoms + values
        for j in range(len(group_order) - 1):
            left, right = x[j] + 0.18, x[j + 1] - 0.18
            control_1, control_2 = left + (right - left) * 0.45, left + (right - left) * 0.55
            path = MplPath(
                [(left, bottoms[j]), (control_1, bottoms[j]), (control_2, bottoms[j + 1]),
                 (right, bottoms[j + 1]), (right, tops[j + 1]),
                 (control_2, tops[j + 1]), (control_1, tops[j]), (left, tops[j]),
                 (left, bottoms[j])],
                [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.CLOSEPOLY],
            )
            ax.add_patch(PathPatch(path, facecolor=colors[cn_id], edgecolor="none", alpha=0.68))
        ax.bar(x, values, bottom=bottoms, width=0.36, color=colors[cn_id],
               edgecolor="black", linewidth=0.55, label=cn_id, zorder=2)
        bottoms = tops
    ax.set(xlabel="", ylabel="Percentage", ylim=(0, 1), xticks=x, xticklabels=group_order)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.tick_params(axis="x", rotation=45, length=0)
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(False)
    sns.despine(ax=ax)
    save_figure(out_dir / f"CNs_{n_clusters}_group_pooled_stacked_proportion", fig, dpi=600)
    return tests


# 功能：逐切片建立 kNN 空间图，计算注释网络中心性和可选 Moran's I。
def spatial_statistics(adata, out_dir, label_col="celltype", moran_genes=(), n_neighbors=10):
    import squidpy as sq

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    score_tables = []
    moran_tables = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample].copy()
        if tmp.n_obs < 2:
            continue
        tmp.obsm["spatial"] = np.asarray(tmp.obsm["spatial_um"])
        sq.gr.spatial_neighbors_knn(tmp, n_neighs=min(n_neighbors, tmp.n_obs - 1))
        sq.gr.centrality_scores(tmp, cluster_key=label_col)
        scores = tmp.uns[f"{label_col}_centrality_scores"].copy().reset_index()
        scores = scores.rename(columns={scores.columns[0]: label_col})
        scores["sample"] = str(sample)
        scores["group"] = str(tmp.obs["group"].iloc[0])
        score_tables.append(scores)

        sq.pl.centrality_scores(tmp, cluster_key=label_col, figsize=(15, 5))
        save_figure(out_dir / f"{sample}_centrality_scores")

        available = [gene for gene in moran_genes if gene in tmp.var_names]
        if available:
            sq.gr.spatial_autocorr(tmp, mode="moran", genes=available)
            moran = tmp.uns["moranI"].reset_index().rename(columns={"index": "gene"})
            moran["sample"] = str(sample)
            moran["group"] = str(tmp.obs["group"].iloc[0])
            moran_tables.append(moran)

    pd.concat(score_tables, ignore_index=True).to_csv(
        out_dir / "centrality_scores_by_sample.tsv", sep="\t", index=False
    )
    if moran_tables:
        pd.concat(moran_tables, ignore_index=True).to_csv(
            out_dir / "moran_by_sample.tsv", sep="\t", index=False
        )
    write_result_readme(
        out_dir,
        "06 Spatial statistics",
        "每个切片独立建立 kNN 空间图，再计算注释网络中心性；配置的基因可额外计算 Moran's I。Moran's I 衡量空间自相关，不是组间差异表达。",
    )


# 功能：逐样本计算每个 bin 周围 k 个邻居的标签组成，防止不同切片坐标重叠造成假邻居。
def neighborhood_composition(adata, label_col="celltype", k=10, coord_key="spatial_um"):
    labels = adata.obs[label_col].astype(str)
    levels = sorted(labels.unique())
    composition = pd.DataFrame(0.0, index=adata.obs_names, columns=levels)

    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        coords = np.asarray(tmp.obsm[coord_key])
        n_neighbors = min(int(k) + 1, tmp.n_obs)
        _, indices = cKDTree(coords).query(coords, k=n_neighbors)
        if indices.ndim == 1:
            indices = indices[:, None]
        local_labels = tmp.obs[label_col].astype(str).to_numpy()
        for row, neighbors in enumerate(indices):
            neighbors = neighbors[neighbors != row][:k]
            counts = pd.Series(local_labels[neighbors]).value_counts(normalize=True)
            composition.loc[tmp.obs_names[row], counts.index] = counts.to_numpy()
    return composition


# 功能：基于局部标签组成聚类 cellular neighborhoods，并输出富集、样本比例和空间图。
def cell_neighborhoods(
    adata,
    out_dir,
    label_col="celltype",
    k=10,
    n_clusters=10,
    group_colors=None,
    seed=0,
):
    import seaborn as sns
    from sklearn.cluster import MiniBatchKMeans

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    composition = neighborhood_composition(adata, label_col, k)
    model = MiniBatchKMeans(n_clusters=n_clusters, random_state=seed, n_init="auto")
    cn = pd.Categorical([f"CN{x}" for x in model.fit_predict(composition)], ordered=False)
    cn_col = f"CNs_{n_clusters}"
    adata.obs[cn_col] = cn

    centers = pd.DataFrame(
        model.cluster_centers_,
        columns=composition.columns,
        index=[f"CN{x}" for x in range(n_clusters)],
    )
    tissue = pd.get_dummies(adata.obs[label_col].astype(str)).reindex(columns=composition.columns, fill_value=0).mean()
    enrichment = np.log2((centers + 1e-6).div(tissue + 1e-6, axis=1))
    composition.assign(**{cn_col: cn}).to_csv(
        out_dir / f"CNs_{n_clusters}_neighbor_composition.tsv", sep="\t"
    )
    enrichment.to_csv(out_dir / f"CNs_{n_clusters}_log2_enrichment.tsv", sep="\t")
    grid = sns.clustermap(
        enrichment,
        cmap="coolwarm",
        center=0,
        vmin=-2,
        vmax=2,
        row_cluster=False,
        col_cluster=False,
        figsize=(max(7, 0.55 * enrichment.shape[1] + 3), 6),
    )
    save_figure(out_dir / f"CNs_{n_clusters}_log2_enrichment", grid.fig)

    counts = pd.crosstab(adata.obs["sample"], adata.obs[cn_col])
    sample_groups = adata.obs[["sample", "group"]].drop_duplicates("sample")
    long = counts.reset_index().melt(id_vars="sample", var_name="CN", value_name="n_cells")
    long = long.merge(sample_groups, on="sample", how="left")
    long["proportion"] = long["n_cells"] / long.groupby("sample")["n_cells"].transform("sum")
    long.to_csv(out_dir / f"CNs_{n_clusters}_proportions_by_sample.tsv", sep="\t", index=False)
    tests = cn_group_comparison(long, out_dir, n_clusters, group_colors)
    tests.to_csv(out_dir / f"CNs_{n_clusters}_group_tests.tsv", sep="\t", index=False)

    cn_palette = dict(zip(adata.obs[cn_col].cat.categories, NPG15[:n_clusters]))
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        fig, ax = plt.subplots(figsize=(7, 7))
        spatial_scatter(tmp, cn_col, ax=ax, size=2, title=f"{sample}: {cn_col}", palette=cn_palette)
        save_figure(out_dir / f"{sample}_spatial_{cn_col}", fig)
    write_result_readme(
        out_dir,
        "07 Cell neighborhood",
        "对每个切片分别计算 kNN 标签组成，再用 MiniBatchKMeans 得到 CN。热图是相对全对象标签背景的 log2 富集；按样本计算CN比例，双组比较使用Welch t-test，另输出组内合并比例alluvial图。",
    )
    return adata


# 功能：逐样本汇总半径邻近比例和距离环带计数，所有标签组合都保留零值。
def proximity_sample_summary(
    adata,
    radius=50,
    distance_bins=(0, 50, 100, 150, 200),
    label_col="celltype",
    coord_key="spatial_um",
):
    levels = sorted(adata.obs[label_col].astype(str).unique())
    radius_rows = []
    annulus_rows = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        coords = np.asarray(tmp.obsm[coord_key])
        labels = tmp.obs[label_col].astype(str).to_numpy()
        group = str(tmp.obs["group"].iloc[0])
        tree = cKDTree(coords)
        radius_neighbors = tree.query_ball_point(coords, radius)
        max_neighbors = tree.query_ball_point(coords, max(distance_bins))

        for source in levels:
            source_indices = np.flatnonzero(labels == source)
            if not len(source_indices):
                continue
            radius_matrix = np.zeros((len(source_indices), len(levels)))
            annulus = {(low, high): np.zeros((len(source_indices), len(levels))) for low, high in zip(distance_bins[:-1], distance_bins[1:])}
            for row, source_idx in enumerate(source_indices):
                neighbors = np.array([x for x in radius_neighbors[source_idx] if x != source_idx], dtype=int)
                if len(neighbors):
                    counts = pd.Series(labels[neighbors]).value_counts(normalize=True)
                    radius_matrix[row] = [counts.get(target, 0) for target in levels]
                candidates = np.array([x for x in max_neighbors[source_idx] if x != source_idx], dtype=int)
                if not len(candidates):
                    continue
                distances = np.linalg.norm(coords[candidates] - coords[source_idx], axis=1)
                for low, high in zip(distance_bins[:-1], distance_bins[1:]):
                    selected = candidates[(distances > low) & (distances <= high)]
                    counts = pd.Series(labels[selected]).value_counts()
                    annulus[(low, high)][row] = [counts.get(target, 0) for target in levels]
            for target, value in zip(levels, radius_matrix.mean(axis=0)):
                radius_rows.append({"sample": sample, "group": group, "source": source, "target": target, "proportion": value})
            for (low, high), matrix in annulus.items():
                for target, value in zip(levels, np.median(matrix, axis=0)):
                    annulus_rows.append({"sample": sample, "group": group, "source": source, "target": target, "distance_range": f"{low}-{high}", "distance_mid": (low + high) / 2, "median_count": value})
    return pd.DataFrame(radius_rows), pd.DataFrame(annulus_rows)


# 功能：绘制样本级邻近比例和距离环带曲线，并保存双组检验表。
def proximity_analysis(adata, out_dir, radius=50, distance_bins=(0, 50, 100, 150, 200), group_colors=None):
    import seaborn as sns

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    radius_table, annulus_table = proximity_sample_summary(adata, radius, distance_bins)
    tests = two_group_tests(radius_table, "proportion", ["source", "target"])
    radius_table.to_csv(out_dir / f"proximity_{radius}um_by_sample.tsv", sep="\t", index=False)
    annulus_table.to_csv(out_dir / "proximity_distance_bins_by_sample.tsv", sep="\t", index=False)
    tests.to_csv(out_dir / f"proximity_{radius}um_group_tests.tsv", sep="\t", index=False)
    palette = group_palette(radius_table["group"].drop_duplicates(), group_colors)

    for source in radius_table["source"].drop_duplicates():
        source_dir = out_dir / str(source)
        source_dir.mkdir(exist_ok=True)
        data = radius_table[radius_table["source"] == source]
        fig, ax = plt.subplots(figsize=(max(9, data["target"].nunique() * 0.65), 5))
        sns.boxplot(data=data, x="target", y="proportion", hue="group", palette=palette, showfliers=False, ax=ax)
        sns.stripplot(data=data, x="target", y="proportion", hue="group", palette=palette, dodge=True, legend=False, ax=ax)
        ax.tick_params(axis="x", rotation=75)
        ax.set(xlabel="", ylabel=f"Neighbor proportion within {radius} µm")
        ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left")
        ax.grid(False)
        save_figure(source_dir / f"proximity_{radius}um", fig)

        annulus = annulus_table[annulus_table["source"] == source]
        grid = sns.relplot(data=annulus, x="distance_mid", y="median_count", hue="group", col="target", col_wrap=4, kind="line", marker="o", errorbar="se", palette=palette, facet_kws={"sharey": False})
        grid.set_axis_labels("Distance from source (µm)", "Median neighboring-bin count")
        save_figure(source_dir / "proximity_distance_ranges", grid.fig)
    write_result_readme(
        out_dir,
        "08 Proximity analysis",
        "在每张切片的 µm 坐标中计算固定半径邻近比例和多个距离环带。原代码的逐 bin 显著性检验改为先按样本汇总。",
    )
    return radius_table, annulus_table, tests


# 功能：计算每个来源 bin 到每一种目标标签的最近距离，且不跨样本寻找邻居。
def nearest_neighbor_table(adata, label_col="celltype", coord_key="spatial_um"):
    rows = []
    levels = sorted(adata.obs[label_col].astype(str).unique())
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        coords = np.asarray(tmp.obsm[coord_key])
        labels = tmp.obs[label_col].astype(str).to_numpy()
        group = str(tmp.obs["group"].iloc[0])
        trees = {target: cKDTree(coords[labels == target]) for target in levels if np.any(labels == target)}
        for idx, source in enumerate(labels):
            for target, tree in trees.items():
                if source == target and np.sum(labels == target) == 1:
                    continue
                k = 2 if source == target else 1
                distances, _ = tree.query(coords[idx], k=k)
                distance = distances[1] if k == 2 else distances
                rows.append({"bin_id": tmp.obs_names[idx], "sample": sample, "group": group, "source": source, "target": target, "distance_um": float(distance)})
    return pd.DataFrame(rows)


# 功能：把最近距离先汇总到样本中位数，再完成双组比较和分面图。
def nearest_neighbor_analysis(adata, out_dir, group_colors=None):
    import seaborn as sns

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    distances = nearest_neighbor_table(adata)
    summary = distances.groupby(["sample", "group", "source", "target"], observed=True)["distance_um"].median().reset_index()
    tests = two_group_tests(summary, "distance_um", ["source", "target"])
    distances.to_csv(out_dir / "nearest_neighbor_distances_by_bin.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "nearest_neighbor_median_by_sample.tsv", sep="\t", index=False)
    tests.to_csv(out_dir / "nearest_neighbor_group_tests.tsv", sep="\t", index=False)
    palette = group_palette(summary["group"].drop_duplicates(), group_colors)
    for source in summary["source"].drop_duplicates():
        data = summary[summary["source"] == source]
        grid = sns.catplot(data=data, x="group", y="distance_um", hue="group", col="target", col_wrap=4, kind="box", palette=palette, sharey=False, legend=False)
        grid.set_axis_labels("", "Median nearest distance (µm)")
        save_figure(out_dir / f"{source}_nearest_neighbor_distance", grid.fig)
    write_result_readme(
        out_dir,
        "09 Nearest-neighbor analysis",
        "每个来源 bin 到每一种目标标签的最近距离均在同一切片内计算，随后用样本中位数进行组间比较。",
    )
    return distances, summary, tests
