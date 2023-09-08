# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt

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
cellLines = ['mdamb231', 'bt474', 'hs578t', 'mdamb453', 'hcc38']
for cellLine in cellLines:
    print(f'Searching {cellLine}')
    adata = adatas[cellLine]
    optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
    optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)

    allOptimalGenes[cellLine] = optimalGenes
    allOptimalCombos[cellLine] = optimalCombos
# %%
for cellLine in cellLines:
    print(cellLine)
    print(allOptimalGenes[cellLine].head())
    print(allOptimalCombos[cellLine].head())
    print('')
# %%
for cellLine in ['mdamb453', 'hcc38']:
    bestGene1 = allOptimalGenes[cellLine]['gene1'].iloc[0]
    bestGene2 = [allOptimalCombos[cellLine]['gene1'].iloc[0], allOptimalCombos[cellLine]['gene2'].iloc[0]]

    visualization.plotHists(adatas[cellLine], bestGene1)
    plt.show()
    visualization.plotExpression(adatas[cellLine], bestGene2)
    plt.show()
#%%
surfaceGenes = dataLoading.cleanSurfaceGenes('../../')

pareto1, pareto2 = {}, {}
for cellLine in ['mdamb453', 'hcc38']:
    optim1 = searchOptimal.mergeSurfaceScores(surfaceGenes, allOptimalGenes[cellLine])
    optim2 = searchOptimal.mergeSurfaceScores(surfaceGenes, allOptimalCombos[cellLine], nGenes = 2)

    pareto1[cellLine] = searchOptimal.findOptimalSurfaceMarkers(optim1)
    pareto2[cellLine] = searchOptimal.findOptimalSurfaceMarkers(optim2, nGenes = 2)
# %%
adata = adatas['mdamb453']
genes = ['THBS1', 'TSPAN8']
X = adata[:, genes].X.toarray()

plt.scatter(X[:,0], X[:, 1], c = adata.obs['leiden'].astype(int))
plt.xlabel(genes[0])
plt.ylabel(genes[1])
# %%
dfExpr = pd.DataFrame([X[:,0], X[:, 1], adata.obs['leiden']]).T
dfExpr.columns = [0, 1, 'leiden']