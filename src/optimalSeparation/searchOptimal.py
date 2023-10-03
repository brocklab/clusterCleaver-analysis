import numpy as np
from sklearn.linear_model import RidgeClassifier
from sklearn import metrics

from tqdm import tqdm
import pandas as pd
import itertools
import warnings
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance

def searchGeneSeparation(adata, surfaceGenes, label = 'leiden', nGenes = 1, nCombos = 10000):
    """
    Scores genes based on separability of predefined clusters. 

    Inputs:
        - adata: Anndata object with .obs consisting of numeric leiden cluster column
        - surfaceGenes: List of surface genes to compare
        - label: Column name in .obs containing label identities
        - nGenes: Number of genes to search through 
    Outputs:
        - dfScores: Modified surface genes dataframe with a new separation score
    """
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))

    if len(surfaceCombos) > nCombos:
        warnings.warn('The number of combos generated is beyond the set maximum number of combos. Was this intentional?')
    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')
    comboScores = {}
    
    for combo in tqdm(surfaceCombos):
        surfaceIdx = np.where(adata.var.index.isin(combo))[0]
        X = adata.X[:, surfaceIdx]
        y = adata.obs[label].astype('int')
        clf = RidgeClassifier().fit(X,y)
        # The score is the accuracy when predicting cluster identity
        score = clf.score(X, y)
        
        y_score = clf.decision_function(X)
        fpr, tpr, _ = metrics.roc_curve(y, y_score)
        auc = metrics.auc(fpr, tpr)

        # Cluster with highest expression
        cluster = pd.DataFrame(X).set_index(y).groupby('leiden').mean().reset_index().idxmax(0)[0]

        isClust = np.where(y == cluster)[0]
        percentAboveThresh = np.sum(clf.predict(X[isClust]) == cluster)/len(isClust)

        medExpr = np.median(X[isClust])

        
        comboScores[combo] = [score, auc, cluster, medExpr, percentAboveThresh]

    dfScores = pd.DataFrame(comboScores).T.reset_index()
    nGenes = len(surfaceCombos[0])
    dfScores.columns = [f'gene{geneNum+1}' for geneNum in range(nGenes)] + ['accuracy', 'auc', 'cluster', 'medExpr', 'percentAboveThresh']
    
    dfScores = dfScores.sort_values(by = 'auc', ascending=False).reset_index(drop=True)
    return dfScores

def searchSeparation2(adata, dfScores, label = 'leiden', metric = 'auc', nGenes = 50):
    dfScores = dfScores.sort_values(by = metric, ascending=False).reset_index(drop=True)
    clusters = dfScores['cluster'].unique()
    assert len(clusters) == 2
    genes1 = dfScores.loc[dfScores['cluster'] == clusters[0], 'gene1'][0:nGenes]
    genes2 = dfScores.loc[dfScores['cluster'] == clusters[1], 'gene1'][0:nGenes]

    surfaceCombos = list(itertools.product(genes1, genes2))

    maxScore = 0
    comboScores = {}
    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')
    for combo in tqdm(surfaceCombos):
        surfaceIdx = np.where(adata.var.index.isin(combo))[0]
        X = adata.X[:, surfaceIdx]
        y = adata.obs[label].astype('int')
        clf = RidgeClassifier(class_weight='balanced').fit(X,y)
        # The score is the accuracy when predicting cluster identity
        score = clf.score(X, y)
        
        y_score = clf.decision_function(X)
        fpr, tpr, _ = metrics.roc_curve(y, y_score)
        auc = metrics.auc(fpr, tpr)

        cluster = pd.DataFrame(X).set_index(y).groupby('leiden').mean().reset_index().idxmax(0)[0]
        
        medExpr = np.median(X)
        
        comboScores[combo] = [score, auc, cluster, medExpr]
        if auc > maxScore:
            maxScore = auc
            # print(f'New max score: {maxScore:0.2g}')
    dfScores = pd.DataFrame(comboScores).T.reset_index()
    nGenes = len(surfaceCombos[0])
    dfScores.columns = [f'gene{geneNum+1}' for geneNum in range(nGenes)] + ['accuracy', 'auc', 'cluster', 'medianExpression']
    
    dfScores = dfScores.sort_values(by = 'auc', ascending=False).reset_index(drop=True)
    return dfScores

def paretoOptimal(costs):
    """
    Finds pareto-optimal points
    Modified from: https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
    
    Inputs:
        - costs: An array of features
        - is_efficient: Array of indices of pareto optimal points
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs>costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    
    return is_efficient

def mergeSurfaceScores(surfaceGenes, dfScores, nGenes = 1):
    """Merges surface scores with separation score ability"""
    for geneNum in range(1, nGenes + 1):
        dfScores = dfScores.merge(surfaceGenes[['gene', 'finalScore']].add_suffix(geneNum), 
                                                            left_on  = f'gene{geneNum}', 
                                                            right_on = f'gene{geneNum}')
    return dfScores

def findOptimalSurfaceMarkers(dfScores, metric = 'auc', nGenes = 1):
    """
    Finds the pareto-optimal genes for surface separation

    Inputs:
        - dfScores: Dataframe with columns for genes (gene) and their associated scores 'finalScre 
        - metric: Metric used to separate clusters
        - nGenes: The number of genes to be used for sorting/flow

    Outputs:
        - optimalSurfaceMarkers: Dataframe of pareto-optimal gens
    """
    costs = np.array(dfScores[[metric] + [f'finalScore{geneNum}'for geneNum in range(1, nGenes + 1)]])
    paretoIdx = paretoOptimal(costs)

    return dfScores.iloc[paretoIdx]

def sliced_wasserstein(X, Y, num_proj = 1000):
    """
    Computes the average sliced wasserstein distance for two arrays

    Inputs:
    X, Y: Input arrays (must have same number o columns)
    num_proj: Number of samples to compute distances

    Outputs:
    Mean wasserstein distance

    Notes:
    This was originally taken from
    https://stats.stackexchange.com/questions/404775/calculate-earth-movers-distance-for-two-grayscale-images
    based on the python optimal transport (POT) package.
    """
    dim = X.shape[1]
    ests = []
    # sample uniformly from the unit sphere
    dir1 = np.random.randn(dim, num_proj)
    dir2 = np.divide(dir1, np.linalg.norm(dir1, axis = 0))

    X0_proj = np.matmul(X, dir2)
    X1_proj = np.matmul(Y, dir2)

    ests = []
    for i in range(num_proj):
        ests.append(wasserstein_distance(X0_proj[:, i], X1_proj[:, i]))
    return np.mean(ests)

def vectorizedSearchSorted(a, b):
    """
    A vectorized implementation of the numpy function `searchsorted`.
    This is a columnwise implementation
    """
    m,n = a.shape
    max_num = np.maximum(a.max() - a.min(), b.max() - b.min()) + 1
    r = max_num*np.arange(a.shape[1])[None,:]
    p = np.searchsorted((a+r).ravel(order = 'F'), (b+r).ravel(order = 'F'), 'right').reshape(-1, n, order = 'F')
    res = p - m*(np.arange(n))
    return res
    # z = res
    # print(z)

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
    u_sorter = np.argsort(u_values, axis = 0)
    v_sorter = np.argsort(v_values, axis = 0)

    all_values = np.concatenate((u_values, v_values), axis = 0)
    all_values.sort(kind='mergesort', axis = 0)

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values, axis = 0)
    # Get the respective positions of the values of u and v among the values of
    # both distributions.

    u_values.sort(axis = 0)
    v_values.sort(axis = 0)

    u_cdf_indices = vectorizedSearchSorted(u_values, all_values[:-1])
    v_cdf_indices = vectorizedSearchSorted(v_values, all_values[:-1])

    # Calculate the CDFs of u and v using their weights, if specified.
    u_cdf = u_cdf_indices / u_values.shape[0]

    v_cdf = v_cdf_indices / v_values.shape[0]


    wd = np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas), axis = 0)

    return wd

def searchExpressionDist(adata, surfaceGenes, label = 'leiden', nGenes = 1, nTopGenes = 50, maxCombos = 10000, topGenes = []):
    """
    Computes the wasserstein (earth mover's distance) metric on gene expression data

    Inputs:
        - adata: Anndata object with .obs consisting of numeric leiden cluster column
        - surfaceGenes: List of surface genes to compare
        - label: Column name in .obs containing label identities
        - nGenes: Number of genes to search through 
    Outputs:
        - dfScores: Modified surface genes dataframe with a new separation score    
    """
    if nGenes > 1 and len(topGenes) == 0:
        dfScores1 = searchExpressionDist(adata, surfaceGenes, label, nGenes = 1, maxCombos = maxCombos, nCombos = 50)
        clusters = dfScores1['cluster'].unique()
        assert len(clusters) == 2
        genes1 = dfScores1.loc[dfScores1['cluster'] == clusters[0], 'gene1'][0:nTopGenes]
        genes2 = dfScores1.loc[dfScores1['cluster'] == clusters[1], 'gene1'][0:nTopGenes]
        surfaceCombos = list(itertools.product(genes1, genes2))

    elif len(topGenes) > 0:
        surfaceGenes = topGenes
        availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
        surfaceCombos = list(itertools.combinations(availableGenes, nGenes))
        
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))


    if len(surfaceCombos) > maxCombos:
        warnings.warn('The number of combos generated is beyond the set maximum number of combos. Was this intentional?')
    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')
    comboScores = []
    expressedClusters = []
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
            if np.median(X0) > np.median(X1):
                cluster = '0'
            else:
                cluster = '1'
        else:
            cluster = -1
        
        if nGenes == 1:
            dist = wasserstein_distance(X0, X1)
        elif nGenes > 1:
            dist = sliced_wasserstein(X0, X1, num_proj = 50)
        comboScores.append(dist)
        expressedClusters.append(cluster)   
    if nGenes == 1:   
        dfScores = pd.DataFrame({'genes': np.array(surfaceCombos).ravel(), 'scores': comboScores, 'cluster': cluster})
    else:
        geneDict = {f'gene{num+1}': np.array(surfaceCombos)[:, num] for num in range(0, nGenes)}
        geneDict['scores'] = comboScores
        dfScores = pd.DataFrame(geneDict)
    return dfScores.sort_values(by = 'scores', ascending = False)