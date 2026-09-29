#!/usr/bin/env python
"""GPU-accelerated normalization, PCA, UMAP and Leiden clustering for Xenium."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")


def resolution_key(value):
    return f"res_{value:g}"


def run_gpu_clustering(
    input_path,
    output_path,
    resolutions=(0.3, 0.6, 0.8, 1.0),
    min_counts=50,
    min_cells=30,
    n_top_genes=2000,
    n_neighbors=15,
    n_pcs=15,
    batch_key="sample",
    device=0,
    seed=10,
):
    import cupy as cp
    import rapids_singlecell as rsc
    import rmm
    import scanpy as sc
    from rmm.allocators.cupy import rmm_cupy_allocator

    cp.cuda.Device(device).use()
    rmm.reinitialize(managed_memory=True, pool_allocator=False, devices=device)
    cp.cuda.set_allocator(rmm_cupy_allocator)

    adata = sc.read_h5ad(input_path)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.layers["counts"] = adata.X.copy()

    rsc.get.anndata_to_GPU(adata)
    rsc.pp.normalize_total(adata, inplace=True)
    rsc.pp.log1p(adata)
    adata.raw = adata.copy()

    hvg_kwargs = {"n_top_genes": n_top_genes, "flavor": "seurat"}
    if batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        hvg_kwargs["batch_key"] = batch_key
    rsc.pp.highly_variable_genes(adata, **hvg_kwargs)
    adata = adata[:, adata.var["highly_variable"]].copy()
    rsc.pp.scale(adata, max_value=10)
    adata.layers["scaled"] = adata.X.copy()
    rsc.pp.pca(adata)
    rsc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=seed,
    )
    rsc.tl.umap(adata, init_pos="random", random_state=seed)
    for value in resolutions:
        rsc.tl.leiden(
            adata,
            resolution=value,
            key_added=resolution_key(value),
            random_state=seed,
        )

    rsc.get.anndata_to_CPU(adata)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(output_path)
    return adata


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="QC-filtered input h5ad")
    parser.add_argument("--output", required=True, help="Clustered output h5ad")
    parser.add_argument("--resolution", nargs="+", type=float, default=[0.3, 0.6, 0.8, 1.0])
    parser.add_argument("--min-counts", type=int, default=50)
    parser.add_argument("--min-cells", type=int, default=30)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--n-pcs", type=int, default=15)
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--device", type=int, default=0, help="CUDA device index")
    parser.add_argument("--seed", type=int, default=10)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_gpu_clustering(
        input_path=args.input,
        output_path=args.output,
        resolutions=args.resolution,
        min_counts=args.min_counts,
        min_cells=args.min_cells,
        n_top_genes=args.n_top_genes,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        batch_key=args.batch_key,
        device=args.device,
        seed=args.seed,
    )
