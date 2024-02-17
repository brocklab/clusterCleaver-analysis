# %%
import numpy as np
from scipy.stats import gaussian_kde, iqr, wasserstein_distance
from scipy.sparse import issparse
import matplotlib.pyplot as plt

# %% Test KDE on anndata
from optimalSeparation import searchOptimal, dataLoading, visualization
import scanpy as sc

adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adatas = dataLoading.processFullAnndata(adataFull)
# %%
for cellLine, adata in adatas.items():
    adata.write_h5ad(f'../../data/h5ads/final/jostner_{cellLine}.h5ad')
    
# %%
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
bhatGenes = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = 'no0').reset_index(drop = True).sort_values(by = 'scores', ascending = True)
# %%
cellLine = 'mdamb231'
bhatGenesNoMod = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = None,
                                                 minCounts = 100,
                                                 scale = True).reset_index(drop = True)

bhatGenesRemove0 = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = 'remove0',
                                                 minCounts = 100,
                                                 scale = True).reset_index(drop = True)

bhatGenesNo0 = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = 'no0',
                                                 minCounts = 100,
                                                 scale = True).reset_index(drop = True)
# %%
emdGenesNoMod = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)
emdGenesRemove0 = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'remove0',
                                                 scale = False).reset_index(drop = True)

emdGenesNo0 = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'no0',
                                                 minCounts = 1,
                                                 scale = False).reset_index(drop = True)
# %%
emdGenes = searchOptimal.searchExpressionDist(adatas[cellLine], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)

# %%
emdGenes = emdGenes.sort_values(by = 'scores', ascending = False).reset_index(drop = True)
bhatGenesRemove0 = bhatGenesRemove0.sort_values(by = 'scores', ascending = True).reset_index(drop = True)
emdGenesRemove0.head(20)

# %%
cellLine = 'mdamb231'
scaleData = lambda x:(x - min(x))/(max(x) - min(x))
adata = adatas[cellLine]
if issparse(adata.X):
    adata.X = adata.X.toarray()
gene = 'TRAC'
# gene = 'TMEM156'
gene = 'SLCO4A1'
# gene = 'ESAM'
# gene = 'MAL2'
# gene = 'PCDH7'
is0 = adata.obs['leiden'] == '0'
is1 = adata.obs['leiden'] == '1'

expression = adata[:, gene].X

# expression = scaleData(expression)

x0 = expression[is0].ravel()
x1 = expression[is1].ravel()

x0, x1 = searchOptimal.modifyEMD(x0, x1, modifier = 'no0')


bScore = searchOptimal.bhattacharyyaHist(x0, x1)
emdScore = wasserstein_distance(x0, x1)
visualization.plotHists(adata, gene = gene, truncate0 = False)
visualization.plotModifiedHists(x0, x1, gene)

# print('Bscore ranking')
# print(bhatGenes.loc[bhatGenes['genes'] == gene])
# print('emd ranking')
# print(emdGenes.loc[emdGenes['genes'] == gene])
# print(f'emd: {emdScore:0.3f} ')
print(f'x0: {len(x0)} \t x1: {len(x1)}')


# %%
modDict = {None: emdGenesNoMod,
           'remove0': emdGenesRemove0,
           'no0': emdGenesNo0}

modVals = {}
for k, df in modDict.items():
    dfVal = df.loc[df['genes'] == gene]
    if dfVal.shape[0]>0:
        rank = dfVal.index.tolist()[0]
    else:
        rank = -1
    x0New, x1New = searchOptimal.modifyEMD(x0, x1, modifier = k)

    score = wasserstein_distance(x0New, x1New)
    print(f'{k}: {rank} \t {score:0.3f}')
    # visualization.plotModifiedHists(x0New, x1New, gene)

# %%
import pandas as pd
import seaborn as sns
x0 = np.random.normal(0, 1, 1000)
x1 = np.random.normal(0, 1, 100)

labels = np.concatenate([np.repeat('x0', len(x0)), np.repeat('x1', len(x1))])
histdf = pd.DataFrame([np.concatenate([x0, x1]), labels]).T
histdf.columns = ['values', 'labels']

print(wasserstein_distance(x0, x1))
plt.hist(x0)
plt.hist(x1)