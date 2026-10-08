from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import save_figure


# 功能：用 Kneedle 在递减凸形 PCA 方差曲线上自动选择拐点。
def select_pcs_by_knee(variance, sensitivity=1.0, fallback=15):
    from kneed import KneeLocator

    values = np.asarray(variance, dtype=float)
    x = np.arange(1, len(values) + 1)
    locator = KneeLocator(
        x,
        values,
        curve="convex",
        direction="decreasing",
        S=float(sensitivity),
    )
    return int(locator.knee) if locator.knee is not None else min(int(fallback), len(values))


# 功能：完成 Normalize/log1p、HVG、PCA、可选 Harmony、邻接图、UMAP 和多分辨率 Leiden。
def cluster_cpu(
    adata,
    out_dir,
    n_top_genes=3000,
    n_neighbors=20,
    n_pcs=None,
    max_pcs=50,
    resolutions=(0.4, 0.5, 0.6, 0.7),
    batch_key="sample",
    use_harmony=True,
    seed=10,
):
    import matplotlib.pyplot as plt
    import scanpy as sc

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    hvg_args = {"n_top_genes": n_top_genes, "flavor": "seurat", "subset": False}
    if batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        hvg_args["batch_key"] = batch_key
    sc.pp.highly_variable_genes(adata, **hvg_args)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    n_comps = min(int(max_pcs), adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps, random_state=seed)

    variance = adata.uns["pca"]["variance_ratio"]
    selected = int(n_pcs) if n_pcs else select_pcs_by_knee(variance)
    selected = min(selected, n_comps)
    table = pd.DataFrame({"PC": np.arange(1, len(variance) + 1), "variance_ratio": variance})
    table["selected_pc"] = selected
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / "pca_variance_and_selected_pc.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(table["PC"], table["variance_ratio"], marker="o", markersize=3)
    ax.axvline(selected, linestyle="--", color="#E64B35", label=f"PC = {selected}")
    ax.set(xlabel="Principal component", ylabel="Variance ratio", title="Automatic PCA knee")
    ax.legend(frameon=False)
    ax.grid(False)
    save_figure(out_dir / "pca_automatic_knee", fig)

    use_rep = "X_pca"
    if use_harmony and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_pca_harmony")
        use_rep = "X_pca_harmony"
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=selected, use_rep=use_rep, random_state=seed)
    sc.tl.umap(adata, min_dist=0.2, random_state=seed)
    for resolution in resolutions:
        sc.tl.leiden(adata, resolution=float(resolution), key_added=f"res_{resolution:g}", random_state=seed)
    adata.uns["selected_n_pcs"] = selected
    return adata
