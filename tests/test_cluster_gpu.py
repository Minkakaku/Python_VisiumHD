import unittest

from cluster_gpu import resolution_key


class GPUClusterNamingTests(unittest.TestCase):
    # 功能：验证 GPU 聚类列名与正式配置一致。
    def test_resolution_key(self):
        self.assertEqual(resolution_key(0.7), "res_0.7")
        self.assertEqual(resolution_key(1.0), "res_1")


if __name__ == "__main__":
    unittest.main()
