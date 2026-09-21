import unittest

from xenium_workflow.cluster_gpu import resolution_key


# 验证 GPU 聚类分辨率对应的列名与下游配置一致，无需运行 GPU。
class GPUClusterNamingTests(unittest.TestCase):
    # 功能：验证 GPU 聚类分辨率列名与下游配置一致。
    def test_resolution_keys_match_downstream_configuration(self):
        self.assertEqual(resolution_key(0.3), "res_0.3")
        self.assertEqual(resolution_key(1.0), "res_1")


if __name__ == "__main__":
    unittest.main()
