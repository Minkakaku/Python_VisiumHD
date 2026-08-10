import unittest

import numpy as np

from visium_hd_workflow.preprocessing import select_pcs_by_knee


class PCAKneeTests(unittest.TestCase):
    # 功能：验证自动 PCA 选择总是返回有效的 1-based PC 数量。
    def test_knee_is_valid_pc_count(self):
        variance = np.array([0.35, 0.20, 0.12, 0.07, 0.04, 0.03, 0.025, 0.02])
        selected = select_pcs_by_knee(variance)
        self.assertGreaterEqual(selected, 1)
        self.assertLessEqual(selected, len(variance))


if __name__ == "__main__":
    unittest.main()
