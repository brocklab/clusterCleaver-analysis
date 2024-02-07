# %%
import numpy as np
from scipy.stats import gaussian_kde, iqr, wasserstein_distance
from scipy.sparse import issparse
import matplotlib.pyplot as plt

# %%
def calculateKDE(y, covarianceFactor = 0.25):
    kernel = gaussian_kde(y)
    kernel.covariance_factor = lambda : covarianceFactor
    kernel._compute_covariance()
    return kernel

def bhattacharyya(p, q):
    return np.sum(np.sqrt(p*q))

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
    # Calculate number of bins using Freedman-Diaconis rule
    # If IQR is 0, use sturges rule
    # Could use other rules to slightly improve score accuracy
    fdB = freedmanDiaconis(p)
    if fdB > 0:
        nBins = int(np.ceil((maxFull - minFull)/fdB))
    else:
        nBins = sturgesRule(p)
    # Calculate and normalize counts
    histRange = (minFull, maxFull)
    hist1, _ = np.histogram(p, bins = nBins, range = histRange)
    hist2, _ = np.histogram(q, bins = nBins, range = histRange)
    hist1 = hist1/sum(hist1)
    hist2 = hist2/sum(hist2)

    bScore = bhattacharyya(hist1, hist2)   
    return -np.log(bScore)
# %%
y1 = np.random.normal(0, 1, 10000)
y2 = np.random.normal(1, 1, 10000)
# %%
# %%timeit
# y1 = np.random.exponential(1, 10000)
# y2 = np.random.exponential(5, 10000)
# y2 = y1.copy()

x = np.linspace(-15, 15, 1000)
kernel1 = calculateKDE(y1)
kernel2 = calculateKDE(y2)

y1KDE = kernel1(x)
y2KDE = kernel2(x)

y1KDE = y1KDE/sum(y1KDE)
y2KDE = y2KDE/sum(y2KDE)
# plt.plot(x, y1KDE/sum(y1KDE))
# plt.plot(x, y2KDE/sum(y2KDE))


bScore = bhattacharyya(y1KDE, y2KDE)
# %%
# %%timeit
bhattacharyyaHist(y1, y2)
# %%
p = np.zeros(1000)
q = np.zeros(1000)
bhattacharyyaHist(p, q)
# %%
%%timeit
emd = wasserstein_distance(y1, y2)
# %%
x1 = np.zeros(1000)
x2 = np.zeros(1000)
# %% Test KDE on anndata
from optimalSeparation import searchOptimal, dataLoading, visualization
import scanpy as sc

adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adatas = dataLoading.processFullAnndata(adataFull)
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
bhatGenes = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = 'no0').reset_index(drop = True).sort_values(by = 'scores', ascending = True)
# %%
emdGenesNoMod = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)
emdGenesRemove0 = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'remove0',
                                                 scale = False).reset_index(drop = True)

emdGenesNo0 = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'no0',
                                                 minCounts = 1,
                                                 scale = False).reset_index(drop = True)
# %%
emdGenes = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)

# %%
emdGenes = emdGenes.sort_values(by = 'scores', ascending = False).reset_index(drop = True)

emdGenes.head(20)

# %%
scaleData = lambda x:(x - min(x))/(max(x) - min(x))
adata = adatas['mdamb231']
if issparse(adata.X):
    adata.X = adata.X.toarray()
gene = 'TSPAN8'
# gene = 'TMEM156'
# gene = 'CLEC2B'
# gene = 'ESAM'
# gene = 'MAL2'
is0 = adata.obs['leiden'] == '0'
is1 = adata.obs['leiden'] == '1'

expression = adata[:, gene].X

# expression = scaleData(expression)

x0 = expression[is0].ravel()
x1 = expression[is1].ravel()

x0, x1 = searchOptimal.modifyEMD(x0, x1, modifier = 'no0')


bScore = searchOptimal.bhattacharyyaHist(x0, x1)
emdScore = wasserstein_distance(x0, x1)
# visualization.plotHists(adata, gene = gene, truncate0 = False)
# visualization.plotModifiedHists(x0, x1, gene)

# print('Bscore ranking')
# print(bhatGenes.loc[bhatGenes['genes'] == gene])
# print('emd ranking')
# print(emdGenes.loc[emdGenes['genes'] == gene])
# print(f'emd: {emdScore:0.3f} ')

x0 = expression[is0].ravel()
x1 = expression[is1].ravel()

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
    visualization.plotModifiedHists(x0New, x1New, gene)

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