# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.stats import wasserstein_distance
import time

from sklearn.linear_model import RidgeClassifier
from sklearn import metrics

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
adata = adatas['hcc38']
leidenResolution = 2
nLeiden = 5
c = 1
while nLeiden != 3:
    leidenResolution /= 1.3
    sc.tl.leiden(adata, resolution= leidenResolution)
    nLeiden = len(adata.obs['leiden'].unique())
    c += 1
    if c > 20:
        leidenResolution = 0
        break
cellLineRes['hcc38'] = leidenResolution
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%
allOptimalGenes, allEMDGenes, allOptimalCombos, allEMDCombos = [], [], [], []
emdGenesDict = {}
for cluster in adata.obs['leiden'].unique():
    allClusts = adata.obs['leiden'].unique().tolist()

    isClust = adata.obs['leiden'] == cluster
    notClust = adata.obs['leiden'] != cluster

    adata.obs['newLeiden'] = 0
    adata.obs.loc[isClust, 'newLeiden'] = 1
    sc.pl.umap(adata, color = 'newLeiden')

    optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'], label = 'newLeiden')
    emdGenes = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], label = 'newLeiden')

    allClusts.remove(cluster)

    optimalGenes.loc[optimalGenes['cluster'] != 1, 'cluster'] = ', '.join(allClusts)
    emdGenes.loc[emdGenes['cluster'] != '1', 'cluster'] = ', '.join(allClusts)

    optimalGenes.loc[optimalGenes['cluster'] == 1, 'cluster'] = cluster
    emdGenes.loc[emdGenes['cluster'] == '1', 'cluster'] = cluster

    
    optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75, label = 'newLeiden')
    emdCombos = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], nGenes = 2, topGenes = emdGenes['genes'][0:75].tolist(), label = 'newLeiden')
    
    optimalCombos['cluster'] = cluster
    emdCombos['cluster'] = cluster




    emdGenesDict[cluster] = emdGenes
    # optimalGenes = optimalGenes.loc[optimalGenes['cluster'] == cluster,]
    # emdGenes = emdGenes.loc[emdGenes['cluster'] == cluster,]
    # optimalCombos = optimalCombos.loc[optimalCombos['cluster'] == cluster,]
    # emdCombos = emdCombos.loc[emdCombos['cluster'] == cluster,]

    allOptimalGenes.append(optimalGenes)
    allEMDGenes.append(emdGenes)
    allOptimalCombos.append(optimalCombos)
    allEMDCombos.append(emdCombos)
# %%
allOptimalGenesConcat = pd.concat(allOptimalGenes)
allEMDGenesConcat = pd.concat(allEMDGenes)
allOptimalCombosConcat = pd.concat(allOptimalCombos)
allEMDCombosConcat = pd.concat(allEMDCombos)
# %%
allOptimalGenesConcat.to_csv('../../data/optimalGenes/hcc38allOptimalGenes.csv')
allEMDGenesConcat.to_csv('../../data/optimalGenes/hcc38allEMDGenes.csv')
allOptimalCombosConcat.to_csv('../../data/optimalGenes/hcc38allOptimalCombos.csv')
allEMDCombosConcat.to_csv('../../data/optimalGenes/hcc38allEMDCombos.csv')
# %%
allOptimalGenes = pd.read_csv('../../data/optimalGenes/hcc38allOptimalGenes.csv')
allEMDGenes = pd.read_csv('../../data/optimalGenes/hcc38allEMDGenes.csv')
allOptimalCombos = pd.read_csv('../../data/optimalGenes/hcc38allOptimalCombos.csv')
allEMDCombos = pd.read_csv('../../data/optimalGenes/hcc38allEMDCombos.csv')
# %%
# %%
expressions = np.arange(0.1, 2.5, 0.01)
props, yields = [], []
adata = adatas['mdamb231']
X = adata[:, 'ESAM'].X.toarray()
cluster = adata.obs['leiden']
n1 = cluster.astype(int).sum()
for expression in expressions:
    aboveThresh = X>=expression
    aboveThresh = aboveThresh.ravel()
    xthresh = X[aboveThresh]
    leidenThresh = cluster[aboveThresh]
    total1 = leidenThresh.astype(int).sum()
    prop1 = total1/len(leidenThresh)
    percentOfClust = total1/n1
    props.append(prop1)
    yields.append(percentOfClust)
# %%
plt.figure()
visualization.plotHists(adata, 'ESAM')

plt.figure()
plt.plot(expressions, props, label = 'Percentage of Identified Cells')
plt.plot(expressions, yields, label = 'Percentage of Total Cluster')
plt.xlabel('ESAM Expression Threshold')
plt.ylabel('Percentage Cluster 1')
plt.legend()
# %%
