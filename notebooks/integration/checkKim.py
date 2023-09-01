# %%
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation.dataLoading import loadClusteredAnndatas, cleanSurfaceGenes
from optimalSeparation.visualization import plotHists
from optimalSeparation.searchOptimal import searchGeneSeparation
# %%
adataKim = sc.read_h5ad('../data/h5ads/kim-postqc-normalized-clustered.h5ad')
# %%
# compute_dimensionality_reductions(adataKim)
# sc.tl.leiden(adataKim, resolution = 0.1)

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sc.pl.umap(adataKim, color = 'cellLine', title='Lung Cancer Cell Lines', ax = ax)
fig.savefig('../figures/runComparison/kimAllLines.png')
# %%
for cellLine in adataKim.obs['cellLine'].unique():
    adataKimSub = adataKim[adataKim.obs['cellLine'] == cellLine]
    compute_dimensionality_reductions(adataKimSub)
    sc.tl.leiden(adataKimSub, resolution = 0.1)
    sc.pl.umap(adataKimSub, color = 'leiden', title=cellLine)
# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
# adataKimSub = adataKim[adataKim.obs['cellLine'] == 'H1299']
# compute_dimensionality_reductions(adataKimSub)
sc.tl.leiden(adataKimSub, resolution = 0.1)
sc.pl.umap(adataKimSub, color = 'CAPNS1', title='H1299', ax = ax)
fig.savefig('../figures/runComparison/kimOptimalGenes.png')
# %%
sc.pl.umap(adataKimSub, color = 'CAPNS1', title=cellLine)
# %%