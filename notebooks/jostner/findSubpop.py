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
# %%

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
genes = ["GSTP1",  "MGP",    "SOD2",   "SPARC",  "PI3",    "EGLN3",  "POSTN",  "TFF1",   "PGAP3",  "F2RL1",  "COL1A2", "NR2F1",  "S100A6", "H3C14", "XAGE1A"]
for gene in genes:
    if gene not in adataFull.var_names:
        continue
    plt.rcParams["figure.figsize"] = (12,12)
    fig, axs = plt.subplots(3,2)
    axs = axs.ravel()
    idx = 0

    vmin = 100
    vmax = -100
    for cellLine, adataSub in adatas.items():
        geneIdx = np.where(adataSub.var.index == gene)[0]
        if np.min(adataSub.X[:, geneIdx]) < vmin:
            vmin = np.min(adataSub.X[:, geneIdx])
        if np.max(adataSub.X[:, geneIdx]) > vmax:
            vmax = np.max(adataSub.X[:, geneIdx])
    for cellLine, adataSub in adatas.items():
        # sc.tl.leiden(adataSub, resolution = cellLineRes[cellLine])
        sc.pl.umap(adataSub, 
                color = [gene], 
                title = cellLine,
                ax = axs[idx],
                vmin = vmin,
                vmax = vmax,
                show = False)
        idx += 1
    fig.suptitle(gene)
    plt.show()
    fig.savefig(f'../../figures/cellTypeChecking/{gene}Expression.png', dpi=500)
# %%
genes = ["GSTP1",  "MGP",    "SOD2",   "SPARC",  "PI3",    "EGLN3",  "POSTN",  "TFF1",   "PGAP3",  "F2RL1",  "COL1A2", "NR2F1",  "S100A6", "XAGE1A"]
fig, axs = plt.subplots(4, 4)
axs = axs.ravel()
idx = 0
plt.rcParams["figure.figsize"] = (20, 20adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata))
for gene in genes:
    sc.pl.violin(adataFull, 
                 gene, 
                 groupby = 'sample',
                 ax = axs[idx],
                 show = False,
                 rotation=90
                 )
    idx += 1

# %%
sc.tl.rank_genes_groups(adata, groupby='sample')
allDfs = []
for sample in adataFull.obs['sample'].unique():
    sampleDf = sc.get.rank_genes_groups_df(adataFull, group = sample)
    sampleDf['sample'] = sample
    allDfs.append(sampleDf)

sampleRankConcat = pd.concat(allDfs)