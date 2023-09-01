# %%
# https://link.springer.com/article/10.1007/s13402-022-00765-7
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation.dataLoading import loadClusteredAnndatas, cleanSurfaceGenes
from optimalSeparation.visualization import plotHists
from optimalSeparation.searchOptimal import searchGeneSeparation
# %%
cellLines = ['MDAMB436']
adatas = {}
for cellLine in cellLines:
    adata = sc.read_h5ad(f'../../data/h5ads/dave{cellLine}-postqc-normalized-clustered.h5ad')
    compute_dimensionality_reductions(adata)
    # sc.tl.leiden(adata, resolution = 0.05)
    # sc.pl.umap(adata, color = 'leiden', title=cellLine)
    adata.var = adata.var.reset_index().set_index('gene_ids').rename_axis(index=None)
    adatas[cellLine] = adata
# %%
cellLineGenes = {'MDAMB231': 'ESAM', 'BT549': 'COL6A2', 'MDAMB436': 'CAV1'}
c = 0

fig, axes = plt.subplots(1,3, figsize=(12, 4))
for cellLine, gene in cellLineGenes.items():
    adata = adatas[cellLine]
    sc.pl.umap(adata, color=gene, use_raw = False, title = f'{cellLine}: {gene} Expression', ax = axes[c], show = 0)
    c += 1

fig.savefig('../figures/runComparison/daveOptimalGenes.png', dpi=500)

# %% Rework BT549
cellLine = 'BT549'
adata = sc.read_h5ad(f'../data/h5ads/dave{cellLine}-postqc-normalized-clustered.h5ad')
adata.var = adata.var.reset_index().set_index('gene_ids').rename_axis(index=None)
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.1)
sc.pl.umap(adata, color = 'leiden', title=cellLine)
# %%
adata = adata[adata.obs['leiden'].isin(['0','1'])]
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.1)
sc.pl.umap(adata, title=cellLine)
# %%
sc.pl.umap(adata, color='COL6A2', use_raw=False)