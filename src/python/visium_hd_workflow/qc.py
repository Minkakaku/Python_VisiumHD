from pathlib import Path

import numpy as np
import pandas as pd

from .qc_plots import plot_qc_review, qc_frame, threshold_candidates
from .runtime import write_result_readme


# 功能：计算基因数、UMI、线粒体/核糖体/HB 比例和文库复杂度。
def add_qc_metrics(adata, mt_prefix="MT-"):
    import scanpy as sc

    names = adata.var_names.str.upper()
    adata.var["mt"] = names.str.startswith(mt_prefix.upper())
    adata.var["ribo"] = names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = names.str.match(r"^HB[^P]")
    percent_top = tuple(x for x in (10, 50, 100) if x < adata.n_vars) or None
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=percent_top,
        log1p=False,
        inplace=True,
    )
    genes = adata.obs["n_genes_by_counts"].to_numpy(float)
    counts = adata.obs["total_counts"].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        adata.obs["log10GenesPerUMI"] = np.where(
            (genes > 0) & (counts > 1),
            np.log10(genes) / np.log10(counts),
            np.nan,
        )
    return adata


# 功能：输出全样本、逐样本、分布和密度 QC 图；保留旧调用接口。
def plot_qc(
    adata,
    out_dir,
    prefix,
    plot_cells_per_sample=50000,
    seed=10,
    references=None,
    candidate_quantile=0.10,
):
    return plot_qc_review(
        adata,
        out_dir,
        prefix,
        plot_cells_per_sample=plot_cells_per_sample,
        seed=seed,
        references=references,
        candidate_quantile=candidate_quantile,
    )


# 功能：完成原始对象 QC，并保存指标、候选阈值、绘图输入与样本汇总。
def data_qc(
    adata,
    out_dir,
    mt_prefix="MT-",
    plot_cells_per_sample=50000,
    seed=10,
    lower_candidate_quantiles=(0.01, 0.025, 0.05, 0.10, 0.15),
    candidate_quantile=0.10,
    min_genes=20,
    max_gene_quantile=0.995,
    max_pct_counts_mt=20,
    min_complexity=0.70,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = add_qc_metrics(adata, mt_prefix)
    frame = qc_frame(adata)
    candidates = threshold_candidates(
        frame,
        lower_candidate_quantiles,
        min_genes,
        max_gene_quantile,
        max_pct_counts_mt,
        min_complexity,
    )
    candidates.to_csv(out_dir / "qc_threshold_candidates.tsv", sep="\t", index=False)
    plot_qc(
        adata,
        out_dir,
        "raw",
        plot_cells_per_sample,
        seed,
        candidates,
        candidate_quantile,
    )
    adata.obs.to_csv(out_dir / "all_bin_features_raw.tsv", sep="\t")
    summary = adata.obs.groupby(["sample", "group"], observed=True).agg(
        bins=("total_counts", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median"),
        median_complexity=("log10GenesPerUMI", "median"),
    )
    summary.to_csv(out_dir / "qc_summary_by_sample.tsv", sep="\t")
    write_result_readme(
        out_dir,
        "01 Data QC",
        "计算每个 8 µm bin 或分割单元的表达量、基因数、线粒体/核糖体比例和文库复杂度。"
        "同时输出全样本与逐样本四联小提琴图、counts-feature 散点图、分布图、密度图、候选阈值保留率图及精确绘图输入。"
        "红色虚线仅为候选阈值参考；过滤前必须人工确认 config/qc_thresholds.tsv。bin 或分割单元不是生物学重复。",
    )
    return adata


# 功能：按每个样本的基因数分位数生成可人工修改的 QC 阈值表。
def write_qc_thresholds(
    adata,
    output_file,
    min_genes=20,
    min_gene_quantile=0.01,
    max_gene_quantile=0.995,
    max_pct_counts_mt=20,
    min_complexity=0.70,
):
    rows = []
    for sample in adata.obs["sample"].drop_duplicates():
        tmp = adata[adata.obs["sample"] == sample]
        rows.append(
            {
                "sample": sample,
                "min_genes": max(
                    min_genes,
                    int(np.quantile(tmp.obs["n_genes_by_counts"], min_gene_quantile)),
                ),
                "max_genes": int(
                    np.quantile(tmp.obs["n_genes_by_counts"], max_gene_quantile)
                ),
                "max_pct_counts_mt": max_pct_counts_mt,
                "min_log10GenesPerUMI": min_complexity,
            }
        )
    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)


# 功能：按样本阈值过滤，并保留过滤前后数量用于审计。
def filter_bins(
    adata_qc,
    threshold_file,
    out_dir,
    plot_cells_per_sample=50000,
    seed=10,
):
    import anndata as ad

    thresholds = pd.read_csv(threshold_file, sep="\t")
    objects = []
    rows = []
    for row in thresholds.itertuples(index=False):
        tmp = adata_qc[adata_qc.obs["sample"].astype(str) == str(row.sample)]
        keep = (
            (tmp.obs["n_genes_by_counts"] >= row.min_genes)
            & (tmp.obs["n_genes_by_counts"] <= row.max_genes)
            & (tmp.obs["pct_counts_mt"] <= row.max_pct_counts_mt)
            & (tmp.obs["log10GenesPerUMI"] >= row.min_log10GenesPerUMI)
        )
        objects.append(tmp[keep].copy())
        rows.append({"sample": row.sample, "before": tmp.n_obs, "after": int(keep.sum())})

    adata = ad.concat(objects, join="inner", merge="same", index_unique=None)
    adata.uns["visium_hd_samples"] = adata_qc.uns.get("visium_hd_samples", {})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_qc(
        adata,
        out_dir,
        "filtered",
        plot_cells_per_sample,
        seed,
        references=thresholds,
    )
    adata.obs.to_csv(out_dir / "all_bin_features_filtered.tsv", sep="\t")
    pd.DataFrame(rows).to_csv(
        out_dir / "filter_summary_by_sample.tsv", sep="\t", index=False
    )
    write_result_readme(
        out_dir,
        "02 Bin filter",
        "按 config/qc_thresholds.tsv 对每个样本分别过滤，并输出过滤后的全样本与逐样本标准 QC 面板。"
        "只设置线粒体比例上限，不再沿用原代码会误删低线粒体 bins 的下限条件。",
    )
    return adata
