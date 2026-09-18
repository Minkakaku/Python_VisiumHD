from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import save_figure


STANDARD_METRICS = (
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "log10GenesPerUMI",
)


def qc_frame(adata) -> pd.DataFrame:
    """提取可审计的逐 bin/分割单元 QC 表。"""
    metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
        "log10GenesPerUMI",
    ]
    if "bin_count" in adata.obs:
        metrics.insert(0, "bin_count")
    missing = [metric for metric in STANDARD_METRICS if metric not in adata.obs]
    if missing:
        raise KeyError(f"QC columns are missing: {', '.join(missing)}")
    frame = adata.obs.loc[:, ["sample", *metrics]].copy()
    frame.insert(0, "unit_id", adata.obs_names.astype(str))
    frame["sample"] = frame["sample"].astype(str)
    return frame


def sample_plot_table(frame: pd.DataFrame, n_per_sample=50000, seed=10) -> pd.DataFrame:
    """逐样本等上限抽样，避免大规模 Visium HD 图形过载。"""
    sampled = []
    for _, sample_frame in frame.groupby("sample", sort=False):
        n = min(len(sample_frame), int(n_per_sample))
        sampled.append(sample_frame.sample(n=n, random_state=seed, replace=False))
    return pd.concat(sampled, ignore_index=True)


def threshold_candidates(
    frame: pd.DataFrame,
    lower_quantiles=(0.01, 0.025, 0.05, 0.10, 0.15),
    min_genes=20,
    max_gene_quantile=0.995,
    max_pct_counts_mt=20,
    min_complexity=0.70,
) -> pd.DataFrame:
    """计算逐样本候选阈值和预计保留率，不修改 AnnData。"""
    rows = []
    for sample, sample_frame in frame.groupby("sample", sort=False):
        max_genes = int(sample_frame["n_genes_by_counts"].quantile(max_gene_quantile))
        for quantile in lower_quantiles:
            quantile = float(quantile)
            min_detected = max(
                int(min_genes),
                int(sample_frame["n_genes_by_counts"].quantile(quantile)),
            )
            expression_pass = (
                (sample_frame["n_genes_by_counts"] >= min_detected)
                & (sample_frame["n_genes_by_counts"] <= max_genes)
                & (sample_frame["pct_counts_mt"] <= float(max_pct_counts_mt))
            )
            complexity_pass = expression_pass & (
                sample_frame["log10GenesPerUMI"] >= float(min_complexity)
            )
            row = {
                "sample": sample,
                "candidate": f"sample_q{quantile * 100:g}",
                "lower_quantile": quantile,
                "min_genes": min_detected,
                "max_genes": max_genes,
                "max_pct_counts_mt": float(max_pct_counts_mt),
                "min_log10GenesPerUMI": float(min_complexity),
                "min_total_counts_reference": float(
                    sample_frame["total_counts"].quantile(quantile)
                ),
                "units_total": int(len(sample_frame)),
                "units_retained_expression_qc": int(expression_pass.sum()),
                "retained_fraction_expression_qc": float(expression_pass.mean()),
                "units_retained_with_complexity": int(complexity_pass.sum()),
                "retained_fraction_with_complexity": float(complexity_pass.mean()),
            }
            if "bin_count" in sample_frame:
                row["min_bin_count_reference"] = float(
                    sample_frame["bin_count"].quantile(quantile)
                )
                row["max_bin_count_reference"] = float(
                    sample_frame["bin_count"].quantile(0.999)
                )
            rows.append(row)
    return pd.DataFrame(rows)


# 优先采用配置的样本色；其余样本按顺序使用 tab10 颜色。
def _palette(samples, configured=None):
    import matplotlib.pyplot as plt

    configured = configured or {}
    fallback = plt.get_cmap("tab10")
    return {
        sample: configured.get(sample, fallback(index % 10))
        for index, sample in enumerate(samples)
    }


# 选择最接近目标分位数的候选阈值，并按样本建立索引。
def _reference_rows(table, quantile):
    if table is None or table.empty:
        return None
    if "lower_quantile" not in table:
        return table.drop_duplicates("sample").set_index("sample")
    available = table["lower_quantile"].dropna().unique()
    if not len(available):
        return table.drop_duplicates("sample").set_index("sample")
    selected = available[np.argmin(np.abs(available.astype(float) - float(quantile)))]
    return table.loc[np.isclose(table["lower_quantile"], selected)].set_index("sample")


# 读取一个样本的阈值数值，缺失样本、列或值时返回 None。
def _value(reference, sample, column):
    if reference is None or sample not in reference.index or column not in reference:
        return None
    value = reference.loc[sample, column]
    if isinstance(value, pd.Series):
        value = value.iloc[0]
    return float(value) if pd.notna(value) else None


# 绘制基因数、计数、线粒体比例和复杂度的小提琴图，并叠加参考阈值。
def _violin(data, samples, colors, reference, title, path, dpi):
    import matplotlib.pyplot as plt

    panels = (
        ("n_genes_by_counts", "Detected genes", ("min_genes", "max_genes")),
        ("total_counts", "Total counts", ("min_total_counts_reference",)),
        ("pct_counts_mt", "Mitochondrial counts (%)", ("max_pct_counts_mt",)),
        ("log10GenesPerUMI", "log10 genes per UMI", ("min_log10GenesPerUMI",)),
    )
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.8))
    for axis, (metric, label, threshold_columns) in zip(axes, panels):
        values = [
            data.loc[data["sample"] == sample, metric].dropna().to_numpy()
            for sample in samples
        ]
        parts = axis.violinplot(
            values,
            positions=np.arange(len(samples)),
            widths=0.78,
            showmedians=True,
            showextrema=False,
        )
        for body, sample in zip(parts["bodies"], samples):
            body.set_facecolor(colors[sample])
            body.set_edgecolor("black")
            body.set_alpha(0.78)
            body.set_linewidth(0.35)
        parts["cmedians"].set_color("black")
        for position, sample in enumerate(samples):
            for column in threshold_columns:
                threshold = _value(reference, sample, column)
                if threshold is not None:
                    axis.hlines(
                        threshold,
                        position - 0.34,
                        position + 0.34,
                        color="#B2182B",
                        linestyle="--",
                        linewidth=1.05,
                    )
        axis.set_xticks(np.arange(len(samples)), samples, rotation=90)
        axis.set(xlabel="", ylabel=label, title=label)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(False)
    fig.suptitle(title, fontsize=14, y=1.01)
    fig.tight_layout()
    save_figure(path, fig, dpi=dpi)


# 绘制计数与基因数、复杂度的关系，按样本配色并叠加参考线。
def _scatter(data, samples, colors, reference, title, path, dpi):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    for sample in samples:
        sample_frame = data.loc[data["sample"] == sample]
        axes[0].scatter(
            sample_frame["total_counts"],
            sample_frame["n_genes_by_counts"],
            s=2.2,
            alpha=0.22,
            linewidths=0,
            rasterized=True,
            color=colors[sample],
            label=sample,
        )
        axes[1].scatter(
            sample_frame["total_counts"],
            sample_frame["log10GenesPerUMI"],
            s=2.2,
            alpha=0.22,
            linewidths=0,
            rasterized=True,
            color=colors[sample],
        )
        min_counts = _value(reference, sample, "min_total_counts_reference")
        if min_counts is not None:
            for axis in axes:
                axis.axvline(min_counts, color=colors[sample], linestyle="--", linewidth=0.9)
        for column in ("min_genes", "max_genes"):
            threshold = _value(reference, sample, column)
            if threshold is not None:
                axes[0].axhline(
                    threshold, color=colors[sample], linestyle="--", linewidth=0.9
                )
        complexity = _value(reference, sample, "min_log10GenesPerUMI")
        if complexity is not None:
            axes[1].axhline(
                complexity, color=colors[sample], linestyle="--", linewidth=0.9
            )
    axes[0].set(xlabel="Total counts", ylabel="Detected genes")
    axes[1].set(xlabel="Total counts", ylabel="log10 genes per UMI")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(False)
    axes[0].legend(frameon=False, markerscale=4)
    fig.suptitle(title, fontsize=14, y=1.01)
    fig.tight_layout()
    save_figure(path, fig, dpi=dpi)


# 展示各样本 QC 指标分布；绘图时截去最高 0.1% 的极端尾部。
def _distributions(data, samples, colors, title, path, dpi):
    import matplotlib.pyplot as plt

    panels = []
    if "bin_count" in data:
        panels.append(("bin_count", "Bins per segmented unit", True))
    panels.extend(
        [
            ("total_counts", "Total counts", True),
            ("n_genes_by_counts", "Detected genes", True),
            ("pct_counts_mt", "Mitochondrial counts (%)", False),
            ("log10GenesPerUMI", "log10 genes per UMI", False),
            ("pct_counts_ribo", "Ribosomal counts (%)", False),
        ]
    )
    if len(panels) < 6 and "pct_counts_hb" in data:
        panels.append(("pct_counts_hb", "Hemoglobin counts (%)", False))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.6))
    for axis, (metric, label, use_log_x) in zip(axes.ravel(), panels):
        for sample in samples:
            values = data.loc[data["sample"] == sample, metric].dropna().to_numpy()
            if not len(values):
                continue
            values = values[values <= np.quantile(values, 0.999)]
            axis.hist(
                values,
                bins=90,
                density=True,
                histtype="step",
                linewidth=1.35,
                color=colors[sample],
                label=sample,
            )
        if use_log_x:
            axis.set_xscale("log")
        axis.set(xlabel=label, ylabel="Density")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(False)
    for axis in axes.ravel()[len(panels) :]:
        axis.set_visible(False)
    axes[0, 0].legend(frameon=False)
    fig.suptitle(title, fontsize=14, y=1.01)
    fig.tight_layout()
    save_figure(path, fig, dpi=dpi)


# 逐样本绘制六边形密度图，并在对应坐标尺度上标记候选阈值。
def _density(data, samples, reference, title, path, dpi):
    import matplotlib.pyplot as plt

    has_bin_count = "bin_count" in data
    ncols = 3 if has_bin_count else 2
    fig, axes = plt.subplots(
        len(samples),
        ncols,
        figsize=(5.1 * ncols, 4.6 * len(samples)),
        squeeze=False,
    )
    for row, sample in enumerate(samples):
        sample_frame = data.loc[data["sample"] == sample]
        panels = [
            (
                np.log10(sample_frame["total_counts"].to_numpy() + 1),
                np.log10(sample_frame["n_genes_by_counts"].to_numpy() + 1),
                "log10(total counts + 1)",
                "log10(detected genes + 1)",
            )
        ]
        if has_bin_count:
            panels.append(
                (
                    sample_frame["bin_count"].to_numpy(),
                    np.log10(sample_frame["total_counts"].to_numpy() + 1),
                    "Bins per segmented unit",
                    "log10(total counts + 1)",
                )
            )
        panels.append(
            (
                np.log10(sample_frame["total_counts"].to_numpy() + 1),
                sample_frame["log10GenesPerUMI"].to_numpy(),
                "log10(total counts + 1)",
                "log10 genes per UMI",
            )
        )
        for column, (x, y, xlabel, ylabel) in enumerate(panels):
            valid = np.isfinite(x) & np.isfinite(y)
            axis = axes[row, column]
            axis.hexbin(x[valid], y[valid], gridsize=80, bins="log", mincnt=1, cmap="Blues")
            axis.set(xlabel=xlabel, ylabel=ylabel, title=f"{sample}: QC density")
            axis.spines[["top", "right"]].set_visible(False)
            axis.grid(False)
        min_counts = _value(reference, sample, "min_total_counts_reference")
        min_genes = _value(reference, sample, "min_genes")
        complexity = _value(reference, sample, "min_log10GenesPerUMI")
        if min_counts is not None:
            axes[row, 0].axvline(np.log10(min_counts + 1), color="#B2182B", linestyle="--")
        if min_genes is not None:
            axes[row, 0].axhline(np.log10(min_genes + 1), color="#B2182B", linestyle="--")
        if complexity is not None:
            axes[row, -1].axhline(complexity, color="#B2182B", linestyle="--")
        if has_bin_count:
            for column, style in (
                ("min_bin_count_reference", "--"),
                ("max_bin_count_reference", ":"),
            ):
                threshold = _value(reference, sample, column)
                if threshold is not None:
                    axes[row, 1].axvline(
                        threshold, color="#B2182B", linestyle=style
                    )
            if min_counts is not None:
                axes[row, 1].axhline(
                    np.log10(min_counts + 1), color="#B2182B", linestyle="--"
                )
    fig.suptitle(title, fontsize=14, y=1.005)
    fig.tight_layout()
    save_figure(path, fig, dpi=dpi)


# 比较不同候选分位数下的保留率，区分是否附加复杂度筛选。
def _retention(candidates, colors, path, dpi):
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8.4, 5.8))
    for sample in candidates["sample"].drop_duplicates():
        sample_frame = candidates.loc[candidates["sample"] == sample].sort_values(
            "lower_quantile"
        )
        x = sample_frame["lower_quantile"].to_numpy() * 100
        axis.plot(
            x,
            sample_frame["retained_fraction_expression_qc"].to_numpy() * 100,
            marker="o",
            color=colors[sample],
            label=f"{sample}: expression/QC",
        )
        axis.plot(
            x,
            sample_frame["retained_fraction_with_complexity"].to_numpy() * 100,
            marker="s",
            linestyle="--",
            color=colors[sample],
            label=f"{sample}: + complexity",
        )
    axis.set(
        xlabel="Sample-specific lower quantile candidate (%)",
        ylabel="Units retained (%)",
        title="QC threshold sensitivity",
        ylim=(0, 102),
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(False)
    axis.legend(frameon=False, ncol=2)
    fig.tight_layout()
    save_figure(path, fig, dpi=dpi)


def plot_qc_review(
    adata,
    out_dir,
    prefix,
    plot_cells_per_sample=50000,
    seed=10,
    references=None,
    candidate_quantile=0.10,
    sample_colors=None,
    dpi=300,
) -> pd.DataFrame:
    """输出全样本、逐样本和密度/分布/保留率 QC 面板。"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = sample_plot_table(qc_frame(adata), plot_cells_per_sample, seed)
    data.to_csv(out_dir / f"{prefix}_qc_plot_input_sampled.tsv", sep="\t", index=False)
    samples = data["sample"].drop_duplicates().tolist()
    colors = _palette(samples, sample_colors)
    reference = _reference_rows(references, candidate_quantile)
    stage = "Unfiltered" if prefix == "raw" else "Filtered"

    _violin(
        data, samples, colors, reference, f"{stage} Visium HD QC",
        out_dir / f"{prefix}_qc_violin", dpi,
    )
    _scatter(
        data, samples, colors, reference,
        f"{stage} Visium HD count-feature relationships",
        out_dir / f"{prefix}_qc_scatter", dpi,
    )
    for sample in samples:
        sample_data = data.loc[data["sample"] == sample]
        sample_reference = (
            reference.loc[[sample]]
            if reference is not None and sample in reference.index
            else None
        )
        safe_sample = re.sub(r"[^0-9A-Za-z._-]+", "_", sample).strip("_") or "sample"
        _violin(
            sample_data, [sample], colors, sample_reference,
            f"{sample}: {stage.lower()} QC",
            out_dir / f"{prefix}_qc_violin_{safe_sample}", dpi,
        )
        _scatter(
            sample_data, [sample], colors, sample_reference,
            f"{sample}: {stage.lower()} count-feature relationships",
            out_dir / f"{prefix}_qc_scatter_{safe_sample}", dpi,
        )
    _distributions(
        data, samples, colors, f"{stage} Visium HD QC distributions",
        out_dir / f"{prefix}_qc_distributions_by_sample", dpi,
    )
    _density(
        data, samples, reference, f"{stage} Visium HD QC density",
        out_dir / f"{prefix}_qc_density_by_sample", dpi,
    )
    if references is not None and "lower_quantile" in references:
        _retention(
            references, colors, out_dir / f"{prefix}_qc_candidate_retention", dpi
        )
    return data
