import unittest

from cluster_gpu import resolution_key


class GPUClusterNamingTests(unittest.TestCase):
    def test_resolution_keys_match_downstream_configuration(self):
        self.assertEqual(resolution_key(0.3), "res_0.3")
        self.assertEqual(resolution_key(1.0), "res_1")


if __name__ == "__main__":
    unittest.main()
