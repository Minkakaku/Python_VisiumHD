"""用逐点距离核对邻近统计：重复坐标、边界、零邻居和跨样本隔离。"""
import anndata as ad
import numpy as np
import pandas as pd
from visium_hd_workflow.spatial import proximity_sample_summary_fast


def test_fast_summary_matches_direct_distances():
    # 第二样本重复第一样本坐标；距离为零的其他细胞保留，自身必须排除。
    coordinates = np.array([[0, 0], [0, 0], [50, 0], [100, 0], [10000, 0]] * 2, dtype=float)
    obs = pd.DataFrame({"sample": ["one"] * 5 + ["two"] * 5,
                        "group": ["CON"] * 5 + ["EXP"] * 5,
                        "celltype": pd.Categorical(["A", "B", "A", "B", "A"] * 2, categories=["A", "B", "absent"])},
                       index=[f"c{i}" for i in range(10)])
    obj = ad.AnnData(obs=obs)
    obj.obsm["spatial"] = coordinates
    radius, rings, cells = proximity_sample_summary_fast(obj, distance_bins=(0, 50, 100), return_cell_proportions=True)
    for sample in ("one", "two"):
        mask = obs["sample"].eq(sample).to_numpy()
        ids = np.flatnonzero(mask)
        distance = np.linalg.norm(coordinates[ids, None, :] - coordinates[None, ids, :], axis=2)
        neighbors = distance <= 50
        np.fill_diagonal(neighbors, False)
        labels = obs.iloc[ids]["celltype"].astype(str).to_numpy()
        counts = np.column_stack([neighbors[:, labels == target].sum(axis=1) for target in ["A", "B", "absent"]])
        denominator = counts.sum(axis=1, keepdims=True)
        proportions = np.divide(counts, denominator, out=np.zeros_like(counts, dtype=float), where=denominator > 0)
        np.testing.assert_allclose(cells.loc[obs.index[ids], ["A", "B", "absent"]], proportions)
        np.testing.assert_array_equal(cells.loc[obs.index[ids], "neighbor_count"], denominator.ravel())
        for source in ("A", "B"):
            for j, target in enumerate(["A", "B", "absent"]):
                row = radius.query("sample == @sample and source_celltype == @source and target_celltype == @target").iloc[0]
                assert np.isclose(row.proportion, proportions[labels == source, j].mean())
                for low, high in [(0, 50), (50, 100)]:
                    ring_counts = ((distance > low) & (distance <= high))[:, labels == target].sum(axis=1)
                    row = rings.loc[(rings["sample"] == sample) & (rings["source_celltype"] == source) &
                                    (rings["target_celltype"] == target) & (rings["distance_mid"] == (low + high) / 2)].iloc[0]
                    assert np.isclose(row.median_count, np.median(ring_counts[labels == source]))
