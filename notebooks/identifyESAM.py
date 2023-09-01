# %%
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn import metrics

from scrna.cluster.main import compute_dimensionality_reductions, compute_leiden_clusters
from optimalSeparation.dataLoading import cleanSurfaceGenes, loadClusteredAnndatas

# from
# %%
adataESAM = sc.read_h5ad('../data/h5ads/231-AA115-final.h5ad')
# %%
# adataESAM = adataESAM[adataESAM.obs['sample'] == 'PT']
sc.pp.highly_variable_genes(adataESAM, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adataESAM)
# %%
sc.pl.umap(adataESAM, color = 'sample')
# %%
sc.tl.leiden(adataESAM, resolution=0.05, key_added='leiden')
sc.pl.umap(adataESAM, color="leiden")
sc.pl.umap(adataESAM, color="ESAM")
# Cluster 1 is ESAM (+)
# %% Ensure there are the same genes
adatas = loadClusteredAnndatas('', '')
# %%
cellLines = list(adatas.keys())
panGenes = adatas[cellLines[0]].var.index.tolist()
genes1KB3 = adataESAM.var.index.tolist()

sharedGenes = list(set(panGenes) & set(genes1KB3))

adataESAM = adataESAM[:, sharedGenes]
for cellLine in cellLines:
    adatas[cellLine] = adatas[cellLine][:, sharedGenes]
# %% Predict with logistic regression
X = adataESAM.X
y = adataESAM.obs['leiden'].astype('int').tolist()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1234, stratify=y)
clf = LogisticRegression(random_state=0).fit(X_train, y_train)
# %%
y_score = clf.decision_function(X_test)
fpr, tpr, _ = metrics.roc_curve(y_test, y_score)
auc = metrics.auc(fpr, tpr)
print(f'AUC: {auc:0.3g}')
# %%

# %%
cellLine = cellLines[0]
adata = adatas[cellLine]
X = adata.X
y_score = clf.decision_function(X)
adata.obs['esamScore'] = y_score

sc.pl.umap(adata, color='esamScore')

# %% Check surface markers
plt.rcParams["figure.figsize"] = (20,20)
fig, axs = plt.subplots(4,4)
axs = axs.ravel()
for idx, cellLine in enumerate(cellLines):
    adata = adatas[cellLine]
    X = adata.X
    y_score = clf.decision_function(X)
    adata.obs['esamScore'] = y_score
    pltTitle = f'{cellLine}'
    sc.pl.umap(adatas[cellLine], color='esamScore', ax = axs[idx], show=False, title=pltTitle)

# %%
for idx, cellLine in enumerate(cellLines):
    adata = adatas[cellLine]
    X = adata.X
    y_score = clf.decision_function(X)
    adata.obs['esamScore'] = y_score
    pltTitle = f'{cellLine}'
    sc.pl.umap(adatas[cellLine], color='esamScore', ax = axs[idx], show=False, title=pltTitle)
