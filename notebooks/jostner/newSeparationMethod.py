# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.stats import wasserstein_distance

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
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])

# %%
# optimal231_1 = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'])
optimal231_2 = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], nGenes=2)

# %%
optimalScores = searchOptimal.searchGeneSeparation(adatas['mdamb231'], surfaceGenes['gene'])
# %%
nGenes = 10
# inter = set(distDf['gene'].tolist()[0:nGenes]) & set(optimalScores['gene1'].tolist()[0:nGenes])
# %%
label = 'leiden'
genes = ['TM9SF3', 'CAV1', 'THBS1']
adata.obs[label] = adata.obs[label].astype("string")
for combo in tqdm([genes]):
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
plt.rcParams["figure.figsize"] = (8, 4)
plt.subplot(121)
plt.plot(threshold, precisions[0:-1], 'b--', label = 'precision')
plt.plot(threshold, recalls[0:-1], 'g', label = 'recall')
plt.xlabel('Threshold values')
plt.legend()
plt.subplot(122)
if X.shape[1] == 2:
    visualization.plotExpression(adata, genes = genes)
elif X.shape[1] == 1:
    visualization.plotHists(adata, gene = genes)
# %% Sliced wasserstein
%%timeit
def sliced_wasserstein(X, Y, num_proj):
    dim = X.shape[1]
    ests = []
    for _ in range(num_proj):
        # sample uniformly from the unit sphere
        dir = np.random.randn(dim)
        dir /= np.linalg.norm(dir)

        # project the data
        X_proj = X @ dir
        Y_proj = Y @ dir

        # compute 1d wasserstein
        ests.append(wasserstein_distance(X_proj, Y_proj))
    return np.mean(ests)
is0 = np.array(adata.obs['leiden'] == '0').astype(bool)
is1 = np.array(adata.obs['leiden'] == '1').astype(bool)

X0 = X[is0, :]
X1 = X[is1, :]

res = sliced_wasserstein(X0, X1, 10000)

print(res)

# %% Vectorize sliced wasserstein
%%timeit
dim = X0.shape[1]
ests = []
# sample uniformly from the unit sphere
dir1 = np.random.randn(dim, 10000)
dir2 = np.divide(dir1, np.linalg.norm(dir1, axis = 0))

X0_proj = np.matmul(X0, dir2)
X1_proj = np.matmul(X1, dir2)

ests = []
for i in range(10000):
    ests.append(wasserstein_distance(X0_proj[:, i], X1_proj[:, i]))
ests = np.mean(ests)
# %%
