#!/usr/bin/env python
"""使用 RAPIDS/CuPy/RMM 加速 Visium HD 聚类，输入与输出均为 H5AD。"""

import argparse
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")


# 功能：把 resolution 转成与下游配置一致的 metadata 列名。
def resolution_key(value):
    return f"res_{value:g}"


# 功能：使用 RAPIDS 完成标准化、HVG、PCA、可选 Harmony、UMAP 与 Leiden。
def run_gpu_clustering(
    input_path,
    output_path,
    result_dir=None,
    resolutions=(0.4, 0.5, 0.6, 0.7),
    n_top_genes=3000,
    n_neighbors=20,
    n_pcs=None,
    max_pcs=50,
    batch_key="sample",
    use_harmony=True,
    regress_keys=(),
    min_dist=0.2,
    device=0,
    seed=10,
):
    import cupy as cp
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import rapids_singlecell as rsc
    import rmm
    import scanpy as sc
    from rmm.allocators.cupy import rmm_cupy_allocator
    from visium_hd_workflow.preprocessing import select_pcs_by_knee
    from visium_hd_workflow.runtime import save_figure, write_result_readme

    cp.cuda.Device(device).use()
    rmm.reinitialize(managed_memory=True, pool_allocator=False, devices=device)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    adata = sc.read_h5ad(input_path)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    rsc.get.anndata_to_GPU(adata, convert_all=True)
    rsc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    rsc.pp.log1p(adata)

    # raw 必须在筛选 HVG 前保存，marker 才能使用完整基因集。
    rsc.get.anndata_to_CPU(adata, convert_all=True)
    adata.raw = adata.copy()
    rsc.get.anndata_to_GPU(adata, convert_all=True)

    hvg_args = {"n_top_genes": n_top_genes, "flavor": "seurat"}
    if batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        hvg_args["batch_key"] = batch_key
    rsc.pp.highly_variable_genes(adata, **hvg_args)
    adata = adata[:, adata.var["highly_variable"].values].copy()
    valid_regress = [key for key in regress_keys if key in adata.obs]
    if valid_regress:
        rsc.pp.regress_out(adata, valid_regress)
    rsc.pp.scale(adata, max_value=10)

    n_comps = min(int(max_pcs), adata.n_obs - 1, adata.n_vars - 1)
    rsc.pp.pca(adata, n_comps=n_comps, random_state=seed)
    rsc.get.anndata_to_CPU(adata, convert_all=True)
    variance = np.asarray(adata.uns["pca"]["variance_ratio"])
    selected = int(n_pcs) if n_pcs else select_pcs_by_knee(variance)
    selected = min(selected, n_comps)

    result_dir = Path(result_dir or Path(output_path).parent)
    result_dir.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame({"PC": np.arange(1, len(variance) + 1), "variance_ratio": variance})
    table["selected_pc"] = selected
    table.to_csv(result_dir / "pca_variance_and_selected_pc.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(table["PC"], table["variance_ratio"], marker="o", markersize=3)
    ax.axvline(selected, linestyle="--", color="#E64B35", label=f"PC = {selected}")
    ax.set(xlabel="Principal component", ylabel="Variance ratio", title="Automatic PCA knee (GPU workflow)")
    ax.legend(frameon=False)
    ax.grid(False)
    save_figure(result_dir / "pca_automatic_knee", fig)

    rsc.get.anndata_to_GPU(adata, convert_all=True)
    use_rep = "X_pca"
    if use_harmony and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        rsc.pp.harmony_integrate(adata, key=batch_key)
        use_rep = "X_pca_harmony"
    rsc.pp.neighbors(
        adata,
        n_neighbors=min(int(n_neighbors), adata.n_obs - 1),
        n_pcs=selected,
        use_rep=use_rep,
        metric="cosine",
        random_state=seed,
    )
    rsc.tl.umap(adata, min_dist=min_dist, random_state=seed)
    for resolution in resolutions:
        rsc.tl.leiden(
            adata,
            resolution=float(resolution),
            key_added=resolution_key(resolution),
            random_state=seed,
        )

    rsc.get.anndata_to_CPU(adata, convert_all=True)
    adata.uns["selected_n_pcs"] = selected
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path)
    write_result_readme(
        result_dir,
        "03 Dimension reduction and clustering",
        "本次使用 RAPIDS GPU 路线完成聚类。PCA 数量由 kneed 自动拐点或命令行参数决定；多样本可启用 Harmony。",
    )
    return adata


# 功能：解析独立 GPU 聚类入口参数。
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="QC-filtered H5AD")
    parser.add_argument("--output", required=True, help="clustered H5AD")
    parser.add_argument("--result-dir", help="PCA plot/table output directory")
    parser.add_argument("--resolution", nargs="+", type=float, default=[0.4, 0.5, 0.6, 0.7])
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--n-pcs", type=int, help="omit to use automatic knee")
    parser.add_argument("--max-pcs", type=int, default=50)
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--no-harmony", action="store_true")
    parser.add_argument("--regress-key", nargs="*", default=[])
    parser.add_argument("--min-dist", type=float, default=0.2)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--seed", type=int, default=10)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_gpu_clustering(
        input_path=args.input,
        output_path=args.output,
        result_dir=args.result_dir,
        resolutions=args.resolution,
        n_top_genes=args.n_top_genes,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        max_pcs=args.max_pcs,
        batch_key=args.batch_key,
        use_harmony=not args.no_harmony,
        regress_keys=args.regress_key,
        min_dist=args.min_dist,
        device=args.device,
        seed=args.seed,
    )
