import numpy as np
from sklearn.linear_model import RidgeClassifier
from sklearn import metrics

from tqdm import tqdm
import pandas as pd
import itertools
import warnings
from scipy.sparse import issparse
from scipy.stats import wasserstein_distance, gaussian_kde, iqr

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
        cluster = pd.DataFrame(X).set_index(y).groupby(label).mean().reset_index().idxmax(0)[0]

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

        cluster = pd.DataFrame(X).set_index(y).groupby(label).mean().reset_index().idxmax(0)[0]
        
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

def searchExpressionDist1D(
                            adata, 
                            surfaceGenes,
                            modifier = 'remove0',
                            metric = 'EMD', 
                            label = 'leiden',
                            minCounts = 100,
                            scale = False,                           
                            maxCombos = 10000, 
                            ):
    """
    Computes a statistical distance metric on gene expression data.

    Inputs:
        - adata: Anndata object with .obs consisting of numeric leiden cluster column
        - surfaceGenes: List of surface genes to compare, must also be in adata.var.index
        - modifier: Will remove 0s for improved EMD score, else run on unmodified gene expression
        - label: Column name in .obs containing label identities
        - nGenes: Number of genes to search through
        - minCounts: Number of counts sufficient for gene to pass when using modifier "remove0"
        - maxCombos: Maximum number of combinations of genes to search for. Becomes relevant with larger numbers. 
    Outputs:
        - dfScores: Modified surface genes dataframe with a new separation score    
    """
    metricDict = {
        'EMD': wasserstein_distance,
        'bhat': bhattacharyyaHist
    }
    
    # Validate input data
    nGenes = 1
    assert modifier in ['remove0', 'no0', None], 'Modifier must be "remove0" or None'
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    assert len(availableGenes) > 0, 'No surface genes match genes in adata.var.index, check both inputs'

    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))

    if len(surfaceCombos) > maxCombos:
        warnings.warn('The number of combos generated is beyond the set maximum number of combos. Was this intentional?')
    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')

    comboScores = []
    expressedClusters = []
    adata.obs[label] = adata.obs[label].astype("string")

    clusterLabels = list(adata.obs[label].unique())
    assert len(clusterLabels) == 2, 'Number of unique labels in adata.obs[label] must be 2'

    is0 = np.array(adata.obs[label] == clusterLabels[0]).astype(bool)
    is1 = np.array(adata.obs[label] == clusterLabels[1]).astype(bool)
    surfaceCombosWrite = []
    for combo in tqdm(surfaceCombos):
        surfaceIdx = np.where(adata.var.index.isin(combo))[0]
        X = adata.X[:, surfaceIdx]

        if sum(X) == 0:
            continue
        if scale:
            X = (X - min(X))/(max(X) - min(X))

        
        X0 = X[is0, :]
        X1 = X[is1, :]

        # print(f"{combo} : {sum(X0)} \t {sum(X1)} \t {surfaceIdx}")
        if sum(X0) == 0 or sum(X1) == 0:
            continue
        if nGenes == 1:
            X0 = X0.ravel()
            X1 = X1.ravel()
            if np.mean(X0) > np.mean(X1):
                cluster = '0'
            else:
                cluster = '1'
        else:
            cluster = -1
        
        if nGenes == 1:
            X0, X1 = modifyEMD(X0, X1, modifier, minCounts=minCounts)
            distFunc = metricDict[metric]
            if len(X0) < minCounts or len(X1) < minCounts:
                continue
            dist = distFunc(X0, X1)
            comboScores.append(dist)
            expressedClusters.append(cluster)
            surfaceCombosWrite.append(combo)

    
    dfScores = pd.DataFrame({'genes': np.array(surfaceCombosWrite).ravel(), 'scores': comboScores, 'cluster': expressedClusters})


    dfScores['genes'] = dfScores['genes'].astype('string')
    return dfScores.sort_values(by = 'scores', ascending = False)

def modifyEMD(X0, X1, modifier = 'remove0', minCounts = 100):
    """
    Selectively removes gene expression with counts of 0s for the higher expressing cluster. 

    Inputs:
    - X0, X1: Numpy arrays of gene expression for two separate clusters
    - modifier: Selects how to remove genes
                Currently only available as "remove0"
    - minCounts: Prevents selection of genes with low cell numbers after removal
    
    Outputs:
    - X0New, X1New: Modified gene expression values
    """
    X0 = X0.copy()
    X1 = X1.copy()
    assert modifier in ['remove0', 'no0', None]
    if modifier not in ['remove0', 'no0']:
        return X0, X1
    # Other checks
    
    if modifier == 'remove0':

        if X1.mean() > X0.mean():
            X1New = X1[X1>0]
            X0New = X0
        elif X0.mean() > X1.mean():
            X0New = X0[X0>0]
            X1New = X1
        else:
            return X0, X1
        if X0New.shape[0] < minCounts or X1New.shape[0] < minCounts:
            return X0, X1
        else:
            return X0New, X1New
    elif modifier == 'no0':
        X1New = X1[X1>0]
        X0New = X0[X0>0]

        # if len(X1New) < minCounts:
        #     X1New = X1
        # if len(X0New) < minCounts:
        #     X0New = X0
        return X0New, X1New

def calculateKDE(y, covarianceFactor = 0.25):
    kernel = gaussian_kde(y)
    kernel.covariance_factor = lambda : covarianceFactor
    kernel._compute_covariance()
    return kernel

def bhattacharyya(p, q):
    return np.sum(np.sqrt(p*q))

def calculateBhattacharyya(X0, X1, ptsEval = 10000):
    """
    Calculates the Bhattacharyya score by:
    - Finding a kernel density estimate
    - Normalizing KDE
    - Computing score

    Inputs:
    - X0: First gene expression vector
    - X1: Second gene expression vector
    - ptsEval: Number of points to evaluate bhjattacharyya score 
        (will decrease speed on increase of value)
    """
    kernel0 = calculateKDE(X0)
    kernel1 = calculateKDE(X1)

    expr = np.linspace(0, max(np.concatenate([X0, X1])), ptsEval)

    X0KDE = kernel0(expr)
    X1KDE = kernel1(expr)

    X0KDE /= sum(X0KDE)
    X1KDE /= sum(X1KDE)

    bScore = bhattacharyya(X0KDE, X1KDE)

    return bScore

freedmanDiaconis = lambda x: 2*iqr(x)/(len(x)**(1/3))
sturgesRule = lambda x: int(np.ceil(np.log2(len(x))+1))
def bhattacharyyaHist(p, q):
    """
    Very fast (vectorized) bhattacharyya coefficient calculator
    Inputs:
    - p: List of observed values
    - q: List of observed values
    Outputs:
    - bScore: Bhattacharyya coefficient
    """
    # Grab relevant information for later calculations
    full = np.concatenate([p, q])
    maxFull = np.max(full)
    minFull = np.min(full)
    # Calculate and normalize counts
    histRange = (minFull, maxFull)
    hist1, _ = np.histogram(p, bins = 'auto', range = histRange)
    hist2, _ = np.histogram(q, bins = 'auto', range = histRange)
    hist1 = hist1/sum(hist1)
    hist2 = hist2/sum(hist2)

    bScore = bhattacharyya(hist1, hist2)    
    return bScore