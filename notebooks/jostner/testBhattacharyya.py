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
adata = adatas['mdamb231']
if issparse(adata.X):
    adata.X = adata.X.toarray()
gene = 'TRAC'

is0 = adata.obs['leiden'] == '0'
is1 = adata.obs['leiden'] == '1'

x0 = adata[is0.tolist(), gene].X.ravel()
x1 = adata[is1.tolist(), gene].X.ravel()

x0, x1 = searchOptimal.modifyEMD(x0, x1)

bScore = bhattacharyyaHist(x0, x1)

visualization.plotHists(adata, gene = gene)
plt.title(f'{gene} Bhattacharyya {bScore:0.3f}')
# %%
bhattacharyyaHist(x0, x1)
# %%
p = x0.copy()
q = x1.copy()
full = np.concatenate([p, q])
maxFull = np.max(full)
minFull = np.min(full)
# Calculate number of bins using Freedman-Diaconis rule
fdB = freedmanDiaconis(p)
nBins = int(np.ceil(np.abs(maxFull - minFull)/fdB))

# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
bhatGenes = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'bhat', 
                                                 modifier = 'remove0')
emdGenes = searchOptimal.searchExpressionDist(adatas['mdamb231'], 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'remove0')
# %%
emdGenes = emdGenes.sort_values(by = 'scores', ascending = False).reset_index(drop = True)

emdGenes.head(20)

# %%
