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
cellLine = 'mdamb231'
adata = adatas['mdamb231'].copy()
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])

# %%
def vectorizedSearchSorted(a, b):
    """
    A vectorized implementation of the numpy function `searchsorted`.
    This is a columnwise implementation
    """
    m,n = a.shape
    max_num = np.maximum(a.max() - a.min(), b.max() - b.min()) + 1
    r = max_num*np.arange(a.shape[1])[None,:]
    print("\t\t np.searchsorted")
    p = np.searchsorted((a+r).ravel(order = 'F'), (b+r).ravel(order = 'F'), 'right').reshape(-1, n, order = 'F')
    res = p - m*(np.arange(n))
    return res
    
def vectorizedWasserstein(u_values, v_values):
    """
    Computes the wasserstein distance for two values. This is based heavily
    on the scipy code, however it does not have weights. 

    Note: Wasserstein distances are currently only evaluated over columns, there is no axis value

    Inputs:
    u_values, v_values: Matrices to be computed

    Outputs:
    distance: Computed columnwise distance
    """
    all_values = np.concatenate((u_values, v_values), axis = 0)
    all_values.sort(kind='mergesort', axis = 0)

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values, axis = 0)
    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    print(f'\tSorting')
    u_values.sort(axis = 0)
    v_values.sort(axis = 0)

    print(f'\tsearch sorting')
    u_cdf_indices = vectorizedSearchSorted(u_values, all_values[:-1])
    v_cdf_indices = vectorizedSearchSorted(v_values, all_values[:-1])
    print("\tDone")
    # Calculate the CDFs of u and v using their weights, if specified.
    u_cdf = u_cdf_indices / u_values.shape[0]

    v_cdf = v_cdf_indices / v_values.shape[0]


    wd = np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas), axis = 0)

    return wd


# %%
surfaceCombos = [['ESAM', 'THBS1']]
label = 'leiden'
nGenes = 2
adata.obs[label] = adata.obs[label].astype("string")

is0 = np.array(adata.obs[label] == '0').astype(bool)
is1 = np.array(adata.obs[label] == '1').astype(bool)
for combo in tqdm(surfaceCombos):
    surfaceIdx = np.where(adata.var.index.isin(combo))[0]
    X = adata.X[:, surfaceIdx]
    X0 = X[is0, :]
    X1 = X[is1, :]
    if nGenes == 1:
        X0 = X0.ravel()
        X1 = X1.ravel()
    break
    if nGenes == 1:
        dist = wasserstein_distance(X0, X1)
    elif nGenes > 1:
        dist = sliced_wasserstein(X0, X1, num_proj = 50)
    comboScores.append(dist)
# %%
num_proj = 10000
def vectorizedSlicedWasserstein(X, Y, num_proj):
    dim = X.shape[1]
    dir = np.random.randn(num_proj, dim)
    dirNorm = dir / np.linalg.norm(dir, axis = 1)[:, None]
    print('Projecting')
    # project the data
    X_proj = X @ dirNorm.T
    Y_proj = Y @ dirNorm.T
    print('Wassersteining')
    wass = vectorizedWasserstein(X_proj, Y_proj)
    print('Done')
    return np.mean(wass)

t = time.time()
res = vectorizedSlicedWasserstein(X0, X1, num_proj)
print(f'Took {time.time() - t:0.4f} seconds')
print(res)
# %%
t = time.time()
dim = X0.shape[1]
ests = []
# sample uniformly from the unit sphere
dir1 = np.random.randn(dim, 10000)
dir2 = np.divide(dir1, np.linalg.norm(dir1, axis = 0))

X0_proj = np.matmul(X0, dir2)
X1_proj = np.matmul(X1, dir2)

ests = []
for i in range(num_proj):
    ests.append(vectorizedWasserstein(X0_proj[:, i][:,None], X1_proj[:, i][:,None]))
ests = np.mean(ests)
print(f'Took {time.time() - t:0.4f} seconds')
print(ests)
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
n_proj = 100
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
t = time.time()
res = sliced_wasserstein(X0, X1, n_proj)
print(f'Took {time.time() - t:0.5f} seconds\nres = {res:.04f}')
# %% Vectorize sliced wasserstein
n_proj = 100
t = time.time()
dim = X0.shape[1]
ests = []
# sample uniformly from the unit sphere
dir1 = np.random.randn(dim, n_proj)
dir2 = np.divide(dir1, np.linalg.norm(dir1, axis = 0))

X0_proj = np.matmul(X0, dir2)
X1_proj = np.matmul(X1, dir2)

ests = []
for i in range(n_proj):
    ests.append(wasserstein_distance(X0_proj[:, i], X1_proj[:, i]))
ests = np.mean(ests)
print(f'Took {time.time() - t:0.5f} seconds\nres = {ests:.04f}')
# %%
