# %%
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# %%
adata_raw = sc.read_h5ad(
    "../../data/h5ads/adataConcatInner.h5ad"
)
adata_raw.layers["logcounts"] = adata_raw.X
adata_raw

# %%
label_key = "cellLine"
batch_key = "batch"
# %%
adata_raw.obs[batch_key].value_counts()
adata = adata_raw
adata.layers["counts"] = adata.X.copy()

adata.X = adata.layers["counts"].copy()
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
n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")
n_batches
# %%
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
# %%
adata_scvi = adata_hvg.copy()
scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key=batch_key)
model_scvi = scvi.model.SCVI(adata_scvi)
# %%
model_scvi.view_anndata_setup()
# %%
max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
max_epochs_scvi
# %%
model_scvi.train()
# %%
