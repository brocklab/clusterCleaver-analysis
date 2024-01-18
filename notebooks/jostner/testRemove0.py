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

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')

samples = adataFull.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
    adataSub = adataSub[adataSub.obs['scDblFinder_class'] == 'singlet']
    sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataSub)
    # sc.pl.umap(adataSub, title = sample)
    adatas[sample] = adataSub
# %%
cellLineRes = {}
for cellLine, adataSub in adatas.items():
    leidenResolution = 2
    nLeiden = 5
    c = 1
    while nLeiden != 2:
        leidenResolution /= 1.3
        sc.tl.leiden(adataSub, resolution= leidenResolution)
        nLeiden = len(adataSub.obs['leiden'].unique())
        c += 1
        if c > 20:
            leidenResolution = 0
            break
    cellLineRes[cellLine] = leidenResolution
    print(adataSub.obs['leiden'].unique())
    print(f'Cell Line: {cellLine} \t Resolution: {leidenResolution} \t nLeiden: {nLeiden}')
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%
adata = adatas['mdamb231']
allEMDGenes = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'])
allEMDCombos = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], nGenes = 2, topGenes = allEMDGenes['genes'])

# %%
allOptimalGenes, allOptimalCombos = [], []
allEMDGenes, allEMDCombos = [], []
allPareto1, allPareto2 = {}, {}
cellLines = ['mdamb231', 'bt474', 'hs578t', 'mdamb453', 'hcc38', 'mdamb436']
for cellLine in cellLines:
    print(f'Searching {cellLine}')
    adata = adatas[cellLine]
    # optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
    # optimalGenesSurface = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes, nGenes = 1)
    # optimalGenesPareto = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface)
    emdGenes = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'])

    # optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)
    # # optimalGenesSurface2 = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalCombos, nGenes = 2)
    # # optimalGenesPareto2 = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2, nGenes = 2)
    # emdCombos = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], nGenes = 2, topGenes = emdGenes['genes'][0:75].tolist())

    # optimalGenes['cellLine'] = cellLine
    # optimalCombos['cellLine'] = cellLine
    emdGenes['cellLine'] = cellLine
    # emdCombos['cellLine'] = cellLine


    # allOptimalGenes.append(optimalGenes)
    # allOptimalCombos.append(optimalCombos)

    allEMDGenes.append(emdGenes)
    # allEMDCombos.append(emdCombos)
# %%
from optimalSeparation.visualization import plotHists
cellLine = 'mdamb436'
plotHists(adatas[cellLine], gene = 'EZR')
# %%
adatas[cellLine][:, ['BST2', 'ESAM']].X
# %%
lessOptimalGenes = emdGenes['genes'].iloc[10:35].tolist()
for gene in lessOptimalGenes:
    plt.figure()
    plotHists(adatas[cellLine], gene = gene)
