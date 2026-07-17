#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')

import os
os.environ['OPENBLAS_NUM_THREADS'] = '16'
import cupy as cp
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import rapids_singlecell as rsc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

import warnings
warnings.filterwarnings("ignore")

plt.style.use('default')  # 使用默认样式
plt.rcParams['figure.facecolor'] = 'white'  # 确保图像背景为白色
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0,            # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)


# In[ ]:


adata = sc.read_h5ad('adata_raw.h5ad')

adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])

sc.pp.filter_cells(adata, min_counts=50)
sc.pp.filter_genes(adata, min_cells=30)
adata = adata[(adata.obs['log10GenesPerUMI'] > 0.9) & (adata.obs['log10GenesPerUMI'] < 0.99) & (adata.obs['pct_counts_mt'] < 20)]

adata.layers['counts'] = adata.X.copy()

rsc.get.anndata_to_GPU(adata)

# Dimension reduce
rsc.pp.normalize_total(adata, inplace=True)
rsc.pp.log1p(adata)
adata.raw = adata.copy()

# HVG
rsc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat",
    batch_key="sample",
)
adata = adata[:, adata.var['highly_variable']]
rsc.pp.scale(adata, max_value=10)
adata.layers["scaled"] = adata.X.copy()

rsc.pp.pca(adata)

# Cluster
rsc.pp.neighbors(adata, n_neighbors=15, n_pcs=15, use_rep='X_pca', random_state=10)
rsc.tl.umap(adata, init_pos='random', random_state=10)
rsc.tl.leiden(adata, resolution=0.6, key_added='res_0.6', random_state=10)
rsc.tl.leiden(adata, resolution=0.8, key_added='res_0.8', random_state=10)
rsc.tl.leiden(adata, resolution=1, key_added='res_1', random_state=10)

with rc_context({'figure.figsize': (6, 6)}):
    sc.pl.umap(
        adata, 
        color=['res_0.6', 'res_0.8', 'res_1'],
        legend_loc='on data', 
        legend_fontsize=9, 
        legend_fontoutline=3,
        frameon=True,
        size=20,
        show=False
    )

rsc.get.anndata_to_CPU(adata)

adata.write('adata_cluster.h5ad')


# In[ ]:


adata_p1 = adata[adata.obs['sample']=='P1_CRC', :]

with rc_context({'figure.figsize': (8, 8)}):
    sc.pl.spatial(
        adata_p1,
        color=['res_0.6'],
        img_key='hires',
        library_id='P1_CRC',
        alpha_img=0.15,
        crop_coord=[15000, 20000, 25000, 20000],
        marker='s'
        # size=1.5
    )

