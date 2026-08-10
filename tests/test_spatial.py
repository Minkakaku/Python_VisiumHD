import unittest

import anndata as ad
import numpy as np
import pandas as pd

from visium_hd_workflow.spatial import (
    nearest_neighbor_table,
    neighborhood_composition,
    proximity_sample_summary,
)


class SpatialIsolationTests(unittest.TestCase):
    # 功能：构造两个坐标重叠的切片，验证算法仍按样本隔离。
    def setUp(self):
        obs = pd.DataFrame(
            {
                "sample": ["CON_1", "CON_1", "EXP_1", "EXP_1"],
                "group": ["CON", "CON", "EXP", "EXP"],
                "celltype": pd.Categorical(["T", "B", "T", "M"]),
            },
            index=["a1", "a2", "b1", "b2"],
        )
        self.adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)
        self.adata.obsm["spatial_um"] = np.array(
            [[0.0, 0.0], [8.0, 0.0], [0.0, 0.0], [80.0, 0.0]]
        )

    # 功能：同坐标的不同样本不能互相成为 CN 邻居。
    def test_neighborhood_composition_does_not_cross_samples(self):
        table = neighborhood_composition(self.adata, k=1)
        self.assertEqual(table.loc["a1", "B"], 1.0)
        self.assertEqual(table.loc["b1", "M"], 1.0)

    # 功能：半径邻近按样本汇总，并为不存在的组合保留零值。
    def test_proximity_keeps_sample_level_zeroes(self):
        radius, annulus = proximity_sample_summary(self.adata, radius=10, distance_bins=(0, 10, 100))
        con = radius.query("sample == 'CON_1' and source == 'T' and target == 'B'")["proportion"].iloc[0]
        exp = radius.query("sample == 'EXP_1' and source == 'T' and target == 'B'")["proportion"].iloc[0]
        self.assertEqual(con, 1.0)
        self.assertEqual(exp, 0.0)
        self.assertFalse(annulus.empty)

    # 功能：目标标签最近距离不跨样本计算。
    def test_nearest_neighbor_does_not_cross_samples(self):
        table = nearest_neighbor_table(self.adata)
        con_t_to_b = table.query("sample == 'CON_1' and source == 'T' and target == 'B'")
        exp_t_to_m = table.query("sample == 'EXP_1' and source == 'T' and target == 'M'")
        self.assertEqual(con_t_to_b["distance_um"].iloc[0], 8.0)
        self.assertEqual(exp_t_to_m["distance_um"].iloc[0], 80.0)


if __name__ == "__main__":
    unittest.main()
