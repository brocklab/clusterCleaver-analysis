# %%
import scanpy as sc
import pandas as pd
import numpy as np

from scrna.cluster.main import compute_dimensionality_reductions
# %%
adata = sc.read_h5ad('../data/h5ads/231-1KB3-20231013.h5ad')
adata = adata[adata.obs['sample'].isin(['PT', 'C2'])]
adata.obs['subsample'] = 'C2'
# %%
sc.pp.highly_variable_genes(adata)
compute_dimensionality_reductions(adata)
# %%
sc.pl.umap(adata, color = 'sample')
# %%
sc.tl.leiden(adata, resolution=0.025/1.5)
sc.pl.umap(adata, color = ['sample', 'leiden', 'ESAM'])
# %%
obs = adata.obs
obs['pheno'] = 'treated'
obs.loc[obs['leiden'] == '2', 'pheno'] = 'esamPos'
obs.loc[obs['leiden'] == '3', 'pheno'] = 'esamNeg'
adata.obs = obs
adata.obs['pheno'].unique()
# %%
x = adata.obs.groupby(['pheno', 'lineage']).size().sort_values(ascending=False).to_frame().reset_index()