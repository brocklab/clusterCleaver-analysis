# %% [markdown]
"""
This notebook will test a new experimental method using EMD to look for genes that may have overlap between subpopulations
but are still characterized by a high expressing population.
"""
# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import wasserstein_distance
import seaborn as sns
from tqdm import tqdm

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed-regressed-clustered.h5ad')
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')

samples = adataFull.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
    adataSub = adataSub[adataSub.obs['scDblFinder_class'] == 'singlet']
    # sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataSub)
    # sc.pl.umap(adataSub, title = sample)
    adatas[sample] = adataSub
# %%
allEMDGenes, allEMDGenesNo0 = {}, {}
cellLines = ['mdamb231', 'bt474', 'hs578t', 'mdamb453', 'mdamb436']
# cellLines = ['bt474']
for cellLine in cellLines:
    print(f'Searching {cellLine}')
    adata = adatas[cellLine].copy()
    adata.X = adata.layers['log1p_norm']
    emdGenes = searchOptimal.searchExpressionDist1D(adata, surfaceGenes['gene'], modifier = None)

    emdGenes.columns = ['genes', 'scoresOld', 'cluster']
    # emdGenes['cellLine'] = cellLine

    allEMDGenes[cellLine] = emdGenes


    emdGenesNew = searchOptimal.searchExpressionDist1D(adata, surfaceGenes['gene'], modifier = 'remove0')

    emdGenesNew.columns = ['genes', 'scoresNew', 'cluster']
    # emdGenesNew['cellLine'] = cellLine
    allEMDGenesNo0[cellLine] = emdGenesNew
# %%
dfEMDGenesNo0   = pd.concat(allEMDGenesNo0, names = ['cellLine', 'idx']).reset_index().drop(columns = 'idx')
dfEMDGenes0      = pd.concat(allEMDGenes, names = ['cellLine', 'idx']).reset_index().drop(columns = 'idx')
dfEMDGenesAll =dfEMDGenesNo0.merge(dfEMDGenes0, on = ['genes', 'cellLine'])

dfEMDGenesAll.head()

dfEMDGenesAll.to_csv('../../data/optimalGenes/allEMDGenesNewOld.csv')
# %%
dfEMDGenes = dfEMDGenesAll.loc[dfEMDGenesAll['cellLine'] == 'mdamb231']
dfEMDGenes = dfEMDGenes[['genes', 'scoresNew']]
dfEMDGenes.coolumns = ['genes', 'score']
dfEMDGenes.to_csv('../../data/emdGenesRankScore231.csv')

dfEMDGenes = dfEMDGenesAll.loc[dfEMDGenesAll['cellLine'] == 'mdamb436']
dfEMDGenes = dfEMDGenes[['genes', 'scoresNew']]
dfEMDGenes.coolumns = ['genes', 'score']
dfEMDGenes.to_csv('../../data/emdGenesRankScore436.csv')
# %%
adata = adatas['mdamb231']
adata.X = adata.layers['log1p_norm']
sc.tl.rank_genes_groups(adata, groupby = 'leiden')
dfGeneRank = sc.get.rank_genes_groups_df(adata, group = '1')
dfGeneRank.to_csv("../../data/topGenes231.csv")

adata = adatas['mdamb436']
adata.X = adata.layers['log1p_norm']
sc.tl.rank_genes_groups(adata, groupby = 'leiden')
dfGeneRank = sc.get.rank_genes_groups_df(adata, group = '1')
dfGeneRank.to_csv("../../data/topGenes436.csv")
# %%
dfEMDGenesAll.sort_values(by = 'scoresNew', ascending=False)
# %%
visualization.plotHists(adatas['bt474'], gene = 'TNFSF10')
sc.pl.umap(adatas['bt474'], color = ['leiden', 'TNFSF10'])
# %%
cellLine = 'mdamb231'
dfEMDGenes = dfEMDGenesAll.loc[dfEMDGenesAll['cellLine'] == cellLine, :]
print('Scores New:')
print(dfEMDGenes.sort_values(by = 'scoresNew', ascending = True).head())
print('Scores Old:')
print(dfEMDGenes.sort_values(by = 'scoresOld', ascending = True).head())
# %%
gene = 'CLEC2B'
dfEMDGenes.loc[dfEMDGenes['genes'] == gene]
visualization.plotHists(adatas[cellLine], gene = gene)

# %%
gene = 'BST2'
cellLine = 'mdamb436'
visualization.plotHists(adatas[cellLine], gene = gene, saveFig='../../figures/emdExemplar/436BST2.png')
# %%
gene = 'IL13RA2'
cellLine = 'mdamb436'
visualization.plotHists(adatas[cellLine], gene = gene, saveFig='../../figures/emdExemplar/436IL13RA2.png')
# %%
adataNo0 = adatas[cellLine]
is0Clust = adataNo0.obs['leiden'] == '0'
is0Expr = adataNo0[:, gene].X == 0
is0Expr = is0Expr.ravel()

adataNo0 = adataNo0[~(np.logical_and(is0Clust.tolist(), is0Expr.tolist())),:]

visualization.plotHists(adataNo0, gene = gene, saveFig='../../figures/emdExemplar/436IL13RA2No0.png')

# %%
