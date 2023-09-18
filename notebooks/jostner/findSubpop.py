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
# %%
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
plt.rcParams["figure.figsize"] = (12,12)
fig, axs = plt.subplots(3,2)
axs = axs.ravel()
idx = 0
for cellLine, adataSub in adatas.items():
    sc.tl.leiden(adataSub, resolution = cellLineRes[cellLine])
    sc.pl.umap(adataSub, 
               color = ['leiden'], 
               title = cellLine,
               ax = axs[idx],
               show = False)
    idx += 1


# %%
sc.tl.rank_genes_groups(adataFull, groupby = 'sample')
# %%
adata = adatas['mdamb436'].copy()
# %%
plt.rcParams["figure.figsize"] = (4,4)
sc.tl.leiden(adata, resolution=0.06)
sc.pl.umap(adata, color = 'leiden', title = 'mdamb436 leiden')

sc.tl.dendrogram(adata, groupby = 'leiden')
sc.pl.dendrogram(adata, groupby = 'leiden')
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%
allOptimalGenes, allOptimalCombos = {}, {}
allPareto1, allPareto2 = {}, {}
cellLines = ['mdamb231', 'bt474', 'hs578t', 'mdamb453', 'hcc38']
for cellLine in cellLines:
    print(f'Searching {cellLine}')
    adata = adatas[cellLine]
    optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
    optimalGenesSurface = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes, nGenes = 1)
    optimalGenesPareto = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface)

    optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)
    optimalGenesSurface2 = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalCombos, nGenes = 2)
    optimalGenesPareto2 = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2, nGenes = 2)

    allOptimalGenes[cellLine] = optimalGenes
    allOptimalCombos[cellLine] = optimalCombos

    allPareto1[cellLine] = optimalGenesPareto
    allPareto2[cellLine] = optimalGenesPareto2
# %%
optimGenesDf = []
optimCombosDf = []
pareto1Df = []
pareto2Df = []
for cellLine in cellLines:
    optimalGenes = allOptimalGenes[cellLine]
    optimalGenes['cellLine'] = cellLine

    optimalCombos = allOptimalCombos[cellLine]
    optimalCombos['cellLine'] = cellLine

    pareto1 = allPareto1[cellLine]
    pareto1['cellLine'] = cellLine

    pareto2 = allPareto2[cellLine]
    pareto2['cellLine'] = cellLine

    assert 0 in optimalGenes.head(10)['cluster']
    assert 1 in optimalGenes.head(10)['cluster']

    optimGenesDf.append(optimalGenes.head(10))
    optimCombosDf.append(optimalCombos.head(100))
    pareto1Df.append(pareto1)
    pareto2Df.append(pareto2)

optimGenesDf = pd.concat(optimGenesDf)
optimCombosDf = pd.concat(optimCombosDf)

pareto1Df = pd.concat(pareto1Df)
pareto2Df = pd.concat(pareto2Df)

optimGenesDf.to_csv('../../data/optimalGenes/jostner/jostnerOptimGenes1.csv', index = None)
optimCombosDf.to_csv('../../data/optimalGenes/jostner/jostnerOptimGenes2.csv', index = None)
pareto1Df.to_csv('../../data/optimalGenes/jostner/jostnerParetoGenes1.csv', index = None)
pareto2Df.to_csv('../../data/optimalGenes/jostner/jostnerParetoGenes2.csv', index = None)

# %%
for cellLine in cellLines:
    print(cellLine)
    print(allOptimalGenes[cellLine].head())
    print(allOptimalCombos[cellLine].head())
    print('')
# %%
plt.rcParams["figure.figsize"] = (6,6)

for cellLine in ['mdamb231', 'bt474', 'hs578t', 'mdamb453', 'hcc38']:

    bestGene1 = allOptimalGenes[cellLine]['gene1'].iloc[0]
    bestGene2 = [allOptimalCombos[cellLine]['gene1'].iloc[0], allOptimalCombos[cellLine]['gene2'].iloc[0]]

    visualization.plotHists(adatas[cellLine], bestGene1)
    plt.title(f'{cellLine} Expression Histogram')
    plt.show()
    visualization.plotExpression(adatas[cellLine], bestGene2)
    plt.title(f'{cellLine} Expression Values')
    plt.show()
# %% ESAM Expression Values for 231s
adata = adatas['mdamb231']
visualization.plotHists(, gene = 'ESAM')
plt.title('MDAMB231 Separation')

sc.pl.umap(adata, color = 'ESAM')

sc.pl.umap(adata, color = 'leiden')


# %%
cellLine = 'mdamb453'
print(f'Searching {cellLine}')
adata = adatas[cellLine]
optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)

allOptimalGenes[cellLine] = optimalGenes
allOptimalCombos[cellLine] = optimalCombos
# %%
sc.pl.umap(adatas['mdamb231'], color = 'ESAM')
# %%
