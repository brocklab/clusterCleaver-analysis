import numpy as np
from sklearn.linear_model import RidgeClassifier
from sklearn import metrics

from tqdm import tqdm
import pandas as pd
import itertools
import warnings

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
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))

    if len(surfaceCombos) > nCombos:
        warnings.warn('The number of combos generated is beyond the set maximum number of combos. Was this intentional?')
    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')
    maxScore = 0
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

        cluster = pd.DataFrame(X).set_index(y).groupby('leiden').mean().reset_index().idxmax(0)[0]
        comboScores[combo] = [score, auc, cluster]
        if auc > maxScore:
            maxScore = auc
            # print(f'New max score: {maxScore:0.2g}')
    dfScores = pd.DataFrame(comboScores).T.reset_index()
    nGenes = len(surfaceCombos[0])
    dfScores.columns = [f'gene{geneNum+1}' for geneNum in range(nGenes)] + ['accuracy', 'auc', 'cluster']
    
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
        clf = RidgeClassifier().fit(X,y)
        # The score is the accuracy when predicting cluster identity
        score = clf.score(X, y)
        
        y_score = clf.decision_function(X)
        fpr, tpr, _ = metrics.roc_curve(y, y_score)
        auc = metrics.auc(fpr, tpr)

        cluster = pd.DataFrame(X).set_index(y).groupby('leiden').mean().reset_index().idxmax(0)[0]
        comboScores[combo] = [score, auc, cluster]
        if auc > maxScore:
            maxScore = auc
            # print(f'New max score: {maxScore:0.2g}')
    dfScores = pd.DataFrame(comboScores).T.reset_index()
    nGenes = len(surfaceCombos[0])
    dfScores.columns = [f'gene{geneNum+1}' for geneNum in range(nGenes)] + ['accuracy', 'auc', 'cluster']
    
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
