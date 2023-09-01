# %%
import scvi
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal
from optimalSeparation.dataLoading import cleanSurfaceGenes
from optimalSeparation import visualization
# %%
adata = sc.read_h5ad('../../data/h5ads/231-1KB3-final.h5ad')
# adata = sc.read_h5ad('../../data/h5ads/231-AA115-final.h5ad')
adata = adata[adata.obs['sample'].isin(['PT'])]
# %%
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.05/2)
sc.pl.umap(adata, color = ['leiden', 'ESAM'])
# %% Get random sample
randPerm = np.random.permutation(range(0, adata.shape[0]))
midPoint = int(len(randPerm)/2)
rand1 = randPerm[0:midPoint]
rand2 = randPerm[midPoint:]
adata.obs['batch'] = 0
obs = adata.obs
obs['batch'].iloc[rand1]= 1
np.unique(obs['batch'], return_counts=True)
adata.obs = obs
adata.obs['batch'] = adata.obs['batch'].astype('category')
sc.pl.umap(adata, color = 'batch')
# %%
batch_key = 'batch'
sc.pp.filter_genes(adata, min_cells=1)
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[batch_key], wspace=0.25)
# %%
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="cell_ranger", batch_key=batch_key
)
# %%
n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")
# %%
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
# %%
adata_scvi = adata_hvg.copy()
# %%
scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_key)
# %%
model_scvi = scvi.model.SCVI(adata_scvi)
# %%
model_scvi.view_anndata_setup()
# %%
max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
# %%
model_scvi.train()
# %%
adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
# %%
sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
sc.tl.umap(adata_scvi)
# %%
sc.pl.umap(adata_scvi, color=[batch_key, 'leiden'], wspace=0.25)
# %%
