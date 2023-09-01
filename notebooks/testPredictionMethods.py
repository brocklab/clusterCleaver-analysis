# %%
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools
from tqdm import tqdm
import pandas as pd

from sklearn import metrics
from sklearn.linear_model import RidgeClassifier, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score
from optimalSeparation.dataLoading import cleanSurfaceGenes
from scrna.cluster.main import compute_dimensionality_reductions
# %%
adata = sc.read_h5ad('../data/h5ads/231-1KB3-final.h5ad')
# adata = sc.read_h5ad('../data/h5ads/231-AA115-final.h5ad')


adata = adata[adata.obs['sample'].isin(['PT']),:]
surfaceGenesDf = cleanSurfaceGenes('..')
# %%
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.05/2)
sc.pl.umap(adata, color = ['leiden', 'ESAM'])
# %%
models = {'ridge': RidgeClassifier, 'logistic': LogisticRegression, 'rf': RandomForestClassifier, 'svc': SVC}
modelRes = {modelName: [] for modelName in models.keys()}
label = 'leiden'
surfaceGenes = surfaceGenesDf['gene']
nGenes = 1
for modelName, model in models.items():
    print(f'Evaluating {modelName}')
    scGenes = np.array(adata.var.index)

    availableGenes = [gene for gene in surfaceGenes if gene in scGenes]
    surfaceCombos = list(itertools.combinations(availableGenes, nGenes))

    print(f'Searching for {len(surfaceCombos)} combinations of {nGenes} gene(s)')
    maxScore = 0
    comboScores = {}
    
    for combo in tqdm(surfaceCombos):
        surfaceIdx = np.where(adata.var.index.isin(combo))[0]
        X = adata.X[:, surfaceIdx]
        y = adata.obs[label].astype('int')
        clf = model().fit(X,y)
        # The score is the accuracy when predicting cluster identity
        scores = cross_val_score(clf, X, y, cv = 3, scoring = 'roc_auc')
        # if modelName != 'rf':
        #     y_score = clf.decision_function(X)
        # else:
        #     y_score = clf.predict_proba(X)[:,1]
        # fpr, tpr, _ = metrics.roc_curve(y, y_score)
        # auc = metrics.auc(fpr, tpr)
        comboScores[combo] = np.mean(scores)
            # print(f'New max score: {maxScore:0.2g}')
    dfScores = pd.DataFrame(comboScores, index = [0]).T.reset_index()
    nGenes = len(surfaceCombos[0])
    dfScores.columns = [f'gene{geneNum+1}' for geneNum in range(nGenes)] + ['auc']
    
    dfScores = dfScores.sort_values(by = 'auc', ascending=False).reset_index(drop=True)
    dfScores['modelName'] = modelName
    modelRes[modelName] = dfScores
modelResDf = pd.concat(modelRes.values())
modelResDf.to_csv('../data/surfacePreds/modelComp.csv')
# %%
