# %%
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# %%
adata_raw = sc.read_h5ad(
    "../../data/h5ads/adataConcatInner.h5ad"
)

# %%
adataJostner = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adataJostner.obs['batch'] = 'jostner'
adataJostner.obs['cellLine'] = adataJostner.obs['sample']
adata_raw = sc.concat([adata_raw, adataJostner])

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
from scrna.cluster.main import compute_dimensionality_reductions

sc.pp.highly_variable_genes(adata_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata_raw)
# %%
adata_raw.write_h5ad('../../data/h5ads/adata_concatInnerNormalized.h5ad')
# %%
adata = adata_raw
adata.obs['cellLine-batch'] = adata.obs['cellLine'].str[0:] + '-' + adata.obs['batch'].str[0:]
fig = plt.figure(figsize = (8, 8))
gs = gridspec.GridSpec(1, 1)
gs.update(wspace = 0.5, hspace = 0.5)
ax1 = plt.subplot(gs[0, 0])

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('UMAP 1', loc = 'left')
ax1.set_ylabel('UMAP 2', loc = 'bottom')


umappts = adata.obsm['X_umap'].copy()
identity = adata.obs['cellLine-batch']
for cellLine in adata.obs['cellLine-batch'].unique():
        isSample = identity == cellLine
        X = umappts[isSample, 0].copy()
        Y = umappts[isSample, 1].copy()

        ax1.scatter(X, Y, s = 2)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()


ax1.arrow(xmin, ymin, 2.33, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 3, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)

ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)

fig.savefig('../../figures/panCancerUMAP.png', dpi = 500, bbox_inches='tight')

# %%
