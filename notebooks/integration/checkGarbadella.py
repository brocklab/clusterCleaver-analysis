# %%
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation.dataLoading import loadClusteredAnndatas, cleanSurfaceGenes
from optimalSeparation.visualization import plotHists
from optimalSeparation.searchOptimal import searchGeneSeparation
# %%
surfaceGenes = cleanSurfaceGenes('..')
panCancerAdatas = loadClusteredAnndatas('.','.', '../data/h5ads/kinder-leiden2ClustLines.pickle')
adata = sc.read_h5ad('../data/h5ads/gambardella-postqc-normalized-clustered.h5ad')
adata.X = adata.X.toarray()
adataSubPan = panCancerAdatas['MDAMB436_BREAST']
# %%
sc.pl.umap(adata, color='cellLine')
# %%
cellLines = ['BT549', 'MDAMB436']
adatas = {}
for cellLine in cellLines:
    adataSubBreast = adata[adata.obs['cellLine'].isin([cellLine])]
    compute_dimensionality_reductions(adataSubBreast)
    sc.tl.leiden(adataSubBreast, resolution = 0.1)
    sc.pl.umap(adataSubBreast, color = 'leiden', title = cellLine)
    adatas[cellLine] = adataSubBreast
# %%
cellLinesPan = {'BT549_BREAST': 'COL6A2', 'MDAMB436_BREAST': 'CAV1', 'NCIH1299_LUNG': 'CAPNS1'}

c = 0

fig, axes = plt.subplots(1,3, figsize=(12, 4))
for cellLine, gene in cellLinesPan.items():
    adata = panCancerAdatas[cellLine]
    sc.pl.umap(adata, color=gene, use_raw = False, title = f'{cellLine}: {gene} Expression', ax = axes[c], show = 0)
    c += 1
fig.savefig('../figures/runComparison/kinkerOptimalGenes.png', dpi=500)
# %%
optimalGenes1_Breast = searchGeneSeparation(adataSubBreast, surfaceGenes['gene']).sort_values(by='auc', ascending=False)
optimalGenes1_Atlas = searchGeneSeparation(adataSubPan, surfaceGenes['gene']).sort_values(by='auc', ascending=False)
# %%
cellLineGenes = {'MDAMB436': 'CAV1', 'BT549': 'COL6A2'}
c = 0

fig, axes = plt.subplots(1,2, figsize=(7, 4))
for cellLine, gene in cellLineGenes.items():
    adata = adatas[cellLine]
    sc.pl.umap(adata, color=gene, use_raw = False, title = f'{cellLine}: {gene} Expression', ax = axes[c], show = 0)
    c += 1

fig.savefig('../figures/runComparison/garbadellaOptimalGenes.png', dpi=500)
# %%
sc.pl.umap(adataSubBreast, color='CAV1')
# %%
plotHists(adataSubBreast, gene = 'CAV1')
# %%
fig, axes = plt.subplots(1,2)
sc.pl.umap(adataSubBreast, color = 'CAV1', show=False, ax=axes[0], title='Breast Database')
sc.pl.umap(adataSubPan, color = 'CAV1', show=False, ax=axes[1], title='Pan-Cancer Database')
# %%
fig, axes = plt.subplots(1,2)
sc.pl.umap(adataSubBreast, color = 'EPB41L3', show=False, ax=axes[0], title='Breast Database')
sc.pl.umap(adataSubPan, color = 'EPB41L3', show=False, ax=axes[1], title='Pan-Cancer Database')
# %%
fig, axes = plt.subplots(1,2)
sc.pl.umap(adataSubBreast, color = 'EREG', show=False, ax=axes[0], title='Breast Database')
sc.pl.umap(adataSubPan, color = 'EREG', show=False, ax=axes[1], title='Pan-Cancer Database')
# %%