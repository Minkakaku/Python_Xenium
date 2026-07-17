import unittest

import anndata as ad
import numpy as np
import pandas as pd

from xenium_workflow.spatial import nearest_neighbor_table, proximity_sample_summary


class SpatialIsolationTests(unittest.TestCase):
    def setUp(self):
        obs = pd.DataFrame(
            {
                "sample": ["CON_1", "CON_1", "EXP_1", "EXP_1"],
                "group": ["CON", "CON", "EXP", "EXP"],
                "celltype": pd.Categorical(["T cell", "B cell", "T cell", "Macrophage"]),
            },
            index=["a1", "a2", "b1", "b2"],
        )
        self.adata = ad.AnnData(X=np.zeros((4, 1)), obs=obs)
        self.adata.obsm["spatial"] = np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [100.0, 0.0]]
        )

    def test_nearest_neighbors_do_not_cross_samples(self):
        table = nearest_neighbor_table(self.adata)
        self.assertEqual(len(table), 4)
        self.assertEqual(sorted(table["distance"].tolist()), [1.0, 1.0, 100.0, 100.0])
        self.assertTrue((table.loc[table["sample"] == "CON_1", "group"] == "CON").all())
        self.assertTrue((table.loc[table["sample"] == "EXP_1", "group"] == "EXP").all())

    def test_proximity_uses_sample_level_graphs_and_zero_counts(self):
        radius, annulus = proximity_sample_summary(
            self.adata,
            radius=2,
            distance_bins=(0, 2, 200),
        )
        con_value = radius.query(
            "sample == 'CON_1' and source_celltype == 'T cell' and target_celltype == 'B cell'"
        )["proportion"].iloc[0]
        exp_value = radius.query(
            "sample == 'EXP_1' and source_celltype == 'T cell' and target_celltype == 'B cell'"
        )["proportion"].iloc[0]
        self.assertEqual(con_value, 1.0)
        self.assertEqual(exp_value, 0.0)
        self.assertFalse(annulus.empty)


if __name__ == "__main__":
    unittest.main()
