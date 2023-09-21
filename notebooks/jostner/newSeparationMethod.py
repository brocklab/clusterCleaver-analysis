# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

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
cellLine = 'mdamb231'
adata = adatas['mdamb231'].copy()
label = 'leiden'

adata.obs[label] = adata.obs[label].astype("string")
for combo in tqdm([['ESAM']]):
    surfaceIdx = np.where(adata.var.index.isin(combo))[0]
    X = adata.X[:, surfaceIdx]
    y = adata.obs[label].astype('int')
    clf = RidgeClassifier().fit(X,y)
    # The score is the accuracy when predicting cluster identity
    score = clf.score(X, y)
    
    y_score = clf.decision_function(X)
    fpr, tpr, _ = metrics.roc_curve(y, y_score)
    auc = metrics.auc(fpr, tpr)

    precisions, recalls, threshold = metrics.precision_recall_curve(y, y_score)
    # Cluster with highest expression
    cluster = pd.DataFrame(X).set_index(y).groupby('leiden').mean().reset_index().idxmax(0)[0]

    isClust = np.where(y == cluster)[0]
    percentAboveThresh = np.sum(clf.predict(X[isClust]) == cluster)/len(isClust)

    medExpr = np.median(X[isClust])

# %%
plt.plot(threshold, precisions[0:-1], 'b--')
plt.plot(threshold, recalls[0:-1], 'g')
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')

searchOptimal.searchGeneSeparation(adatas['mdamb231'], surfaceGenes['gene'])
# %%
