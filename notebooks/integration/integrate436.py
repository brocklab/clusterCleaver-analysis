# %%
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import anndata

from scrna.cluster.main import compute_dimensionality_reductions
# %%
dave = sc.read_h5ad('../data/h5ads/daveMDAMB436-postqc-normalized-clustered.h5ad')
dave = dave[dave.obs['cellLine'].str.startswith('MDAMB436')]
dave.obs['batch'] = 'dave'
dave.var.index = dave.var['gene_ids']
dave.var.index = dave.var.index.astype(str)
dave.var_names_make_unique()
# %%
gambardella = sc.read_h5ad('../data/h5ads/gambardella-postqc-normalized-clustered.h5ad')
gambardella = gambardella[gambardella.obs['cellLine'].str.startswith('MDAMB436')]
gambardella.obs['batch'] = 'gambardella'
# %%
kinker = sc.read_h5ad('../data/h5ads/kinker-postqc-normalized-clustered.h5ad')
kinker.obs['cellLine'] = kinker.obs['Cell_line']
kinker.obs['cellLine'] = kinker.obs['cellLine'].replace('MDAMB436_BREAST', 'MDAMB436')
kinker = kinker[kinker.obs['cellLine'].str.startswith('MDAMB436')]
kinker.obs['batch'] = 'kinker'
# %%
adata = anndata.concat([dave, gambardella, kinker])
# %%
label_key = 'cellLine'
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
sc.pl.umap(adata, color=[label_key, batch_key], wspace=0.25)
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
sc.pl.umap(adata_scvi, color=[label_key, batch_key], wspace=0.25)
# %%
