import numpy as np
import pandas as pd

from visium_hd_workflow.qc_plots import (
    plot_qc_review,
    sample_plot_table,
    threshold_candidates,
)


def _frame():
    rows = []
    for sample_index, sample in enumerate(("sample1", "sample2")):
        for index in range(20):
            rows.append(
                {
                    "unit_id": f"{sample}:{index}",
                    "sample": sample,
                    "bin_count": 2 + index % 6,
                    "total_counts": 30 + sample_index * 10 + index * 7,
                    "n_genes_by_counts": 12 + sample_index * 3 + index * 2,
                    "pct_counts_mt": 1 + index * 0.5,
                    "pct_counts_ribo": 5 + index * 0.2,
                    "pct_counts_hb": index * 0.1,
                    "log10GenesPerUMI": 0.65 + index * 0.01,
                }
            )
    return pd.DataFrame(rows)


def test_sampling_is_per_sample_and_deterministic():
    frame = _frame()
    first = sample_plot_table(frame, n_per_sample=5, seed=1001)
    second = sample_plot_table(frame, n_per_sample=5, seed=1001)
    assert first.groupby("sample").size().to_dict() == {"sample1": 5, "sample2": 5}
    pd.testing.assert_frame_equal(first, second)


def test_candidates_include_retention_and_optional_bin_references():
    candidates = threshold_candidates(
        _frame(),
        lower_quantiles=(0.10,),
        min_genes=10,
        min_complexity=0.8,
    )
    assert len(candidates) == 2
    assert {"min_bin_count_reference", "max_bin_count_reference"}.issubset(candidates)
    assert (
        candidates["retained_fraction_with_complexity"]
        <= candidates["retained_fraction_expression_qc"]
    ).all()
    assert candidates["min_log10GenesPerUMI"].eq(0.8).all()


def test_qc_review_writes_standard_and_per_sample_panels(tmp_path):
    frame = _frame().set_index("unit_id")

    class DummyAdata:
        obs = frame
        obs_names = frame.index

    candidates = threshold_candidates(
        frame.reset_index(),
        lower_quantiles=(0.10,),
        min_genes=10,
        min_complexity=0.8,
    )
    plot_qc_review(
        DummyAdata(),
        tmp_path,
        "raw",
        plot_cells_per_sample=10,
        seed=1001,
        references=candidates,
        dpi=50,
    )
    stems = (
        "raw_qc_violin",
        "raw_qc_scatter",
        "raw_qc_violin_sample1",
        "raw_qc_scatter_sample1",
        "raw_qc_violin_sample2",
        "raw_qc_scatter_sample2",
        "raw_qc_distributions_by_sample",
        "raw_qc_density_by_sample",
        "raw_qc_candidate_retention",
    )
    for stem in stems:
        assert (tmp_path / f"{stem}.png").exists()
        assert (tmp_path / f"{stem}.pdf").exists()
    assert (tmp_path / "raw_qc_plot_input_sampled.tsv").exists()
