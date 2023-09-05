# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal
from optimalSeparation import dataLoading
# %%
adataFull = sc.read_h5ad('../data/h5ads/jostner-processed.h5ad')
# %%
samples = adataFull.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
    sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataSub)
    sc.pl.umap(adataSub, title = sample)
    adatas[sample] = adataSub
# %%
cellLineRes = {}
for cellLine, adataSub in adatas.items():
    leidenResolution = 0.5
    nLeiden = 5
    c = 1
    while nLeiden != 2:
        sc.tl.leiden(adataSub, resolution= leidenResolution)
        leidenResolution /= 1.2
        nLeiden = len(adataSub.obs['leiden'].unique())
        c += 1
        if c > 20:
            leidenResolution = 0
            break
    cellLineRes[cellLine] = leidenResolution
# %%
plt.rcParams["figure.figsize"] = (10,10)
fig, axs = plt.subplots(3,2)
axs = axs.ravel()
idx = 0
for cellLine, adataSub in adatas.items():
    sc.tl.leiden(adataSub, resolution= cellLineRes[cellLine])
    sc.pl.umap(adataSub, 
               color = ['total_counts'], 
               title = cellLine,
               ax = axs[idx],
               show = False)
    idx += 1


# %%
adata = adatas['mdamb436'].copy()

sc.tl.leiden(adata, resolution=0.1)
sc.pl.umap(adata, color = 'leiden')

adata = adata[adata.obs['leiden'].isin(['0', '1'])]
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('..')
# %%
optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
# %%
optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)
# %%
