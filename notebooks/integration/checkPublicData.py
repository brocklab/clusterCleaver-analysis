# %%
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

from scrna.cluster.main import compute_dimensionality_reductions

# %%
adata = sc.read_h5ad(
    "../../data/h5ads/adataConcatInner.h5ad"
)

# %%
adataJostner = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adataJostner.obs['batch'] = 'jostner'
adataJostner.obs['cellLine'] = adataJostner.obs['sample']
# %%
adata = sc.concat([adata, adataJostner])
# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
# %%
sc.pl.umap(adata, color = 'batch')
# %%
adata.obs['cellLine'] = adata.obs['cellLine'].str.lower()
adata.obs['cellBatch'] = adata.obs['cellLine'].astype('string') + adata.obs['batch'].astype('string')

# %%
batches = ['dave231', 'dave436', 'gambardella', '1kb3', 'AA115', 'jost', 'jostner']
adataUseful = adata[adata.obs['batch'].astype('string').isin(batches)]
# %%
sc.pp.highly_variable_genes(adataUseful, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adataUseful)
# %%
adataUseful.obs['cellBatch'] = adataUseful.obs['cellBatch'].astype('category')
sc.tl.dendrogram(adataUseful, groupby = 'cellBatch')

sc.pl.dendrogram(adataUseful, groupby = 'cellBatch')
# %%
