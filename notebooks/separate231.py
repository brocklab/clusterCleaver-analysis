# %%
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal
from optimalSeparation.dataLoading import cleanSurfaceGenes
from optimalSeparation import visualization
# %%
# adata = sc.read_h5ad('../data/h5ads/231-1KB3-final.h5ad')
adata = sc.read_h5ad('../data/h5ads/231-AA115-final.h5ad')

adata = adata[adata.obs['sample'].isin(['PT']),:]
# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.05/2)
sc.pl.umap(adata, color = ['leiden', 'ESAM'])
# %%
surfaceGenes = cleanSurfaceGenes('..')
# %%
optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
optimalGenesSurface = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes, nGenes = 1)
optimalGenesPareto = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface)
# %%
optimalGenes2 = searchOptimal.searchGeneSeparation(adata, optimalGenes['gene1'][0:100], nGenes = 2, nCombos = 10000)
optimalGenesSurface2 = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes2, nGenes = 2)
optimalGenesPareto2 = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2, nGenes = 2)
# %%
visualization.plotParetoOptimal(optimalGenesSurface, optimalGenesPareto, nGenes = 1)
plt.title('231 Pareto Optimal Genes')
plt.savefig('../figures/231ParetoOptimal1.png', dpi=500)

visualization.plotParetoOptimal(optimalGenesSurface2, optimalGenesPareto2, nGenes = 2)
plt.title('231 Pareto Optimal Genes')
plt.savefig('../figures/231ParetoOptimal2.png', dpi=500)

# %%
gene = 'ESAM'
# gene = 'SLC9A3R2'
# gene = 'BST2'
print(optimalGenes.loc[optimalGenes['gene1'].str.startswith(gene)])
labels = adata.obs['leiden'].astype(int).tolist()
visualization.plotHists(adata, gene)
sc.pl.umap(adata, color = gene)
plt.scatter(np.linspace(0, 1, adata.shape[0]), adata[:, gene].X, c = labels)
# %%
gene1 = 'BST2'
gene2 = 'ESAM'
labels = adata.obs['leiden'].astype(int).tolist()
plt.scatter(adata[:, gene1].X, adata[:, gene2].X, c = labels, label = labels, s = 10)
# plt.legend()
# %%
optimalSurfaceMarkers = optimalGenes.merge(surfaceGenes[['gene', 'finalScore']], left_on= 'gene1', right_on = 'gene')
# %%
plt.scatter(optimalSurfaceMarkers['auc'], optimalSurfaceMarkers['finalScore'])
plt.xlabel('AUC')
plt.ylabel('Surface Score')
# %%
costs = np.array(optimalSurfaceMarkers[['auc', 'finalScore']])
paretoIdx = searchOptimal.paretoOptimal(costs)
paretoSurface = optimalSurfaceMarkers.iloc[paretoIdx]
plt.scatter(optimalSurfaceMarkers['auc'], optimalSurfaceMarkers['finalScore'], label = 'All')
plt.scatter(paretoSurface['auc'], paretoSurface['finalScore'], c = 'red', label = 'Pareto\nOptimal')
plt.legend()

plt.xlabel('AUC')
plt.ylabel('Surface Score')
plt.title('MDAMB231 AA115 Pareto Optimal Genes')
# %%
optimalGenes2 = searchOptimal.searchGeneSeparation(adata, optimalGenes['gene1'].tolist()[0:100], nGenes = 2)
# %%
plt.scatter(adata[:, 'ESAM'].X, adata[:, 'CD63'].X, c = adata.obs['leiden'].astype(int).tolist())
# %% Merging multiple gene calculations
nGenes = 2
optimalSurfaceMarkers2 = optimalGenes2
for geneNum in range(1, nGenes + 1):
    optimalSurfaceMarkers2 = optimalSurfaceMarkers2.merge(surfaceGenes[['gene', 'finalScore']].add_suffix(geneNum), 
                                                            left_on  = f'gene{geneNum}', 
                                                            right_on = f'gene{geneNum}')

optimalSurfaceMarkers2
# %%
costs = np.array(optimalSurfaceMarkers2[['auc'] + [f'finalScore{geneNum}'for geneNum in range(1, nGenes + 1)]])
paretoIdx = searchOptimal.paretoOptimal(costs)
paretoSurface2 = optimalSurfaceMarkers2.iloc[paretoIdx]

ax = plt.figure(figsize=(5,5)).add_subplot(projection='3d')
ax.scatter(optimalSurfaceMarkers2['finalScore1'], optimalSurfaceMarkers2['finalScore2'], optimalSurfaceMarkers2['auc'], c = 'blue', alpha = 0.05)
ax.scatter(paretoSurface2['finalScore1'], paretoSurface2['finalScore2'], paretoSurface2['auc'], c = 'red')
ax.dist = 12
ax.set_xlabel('Surface Score 1')
ax.set_ylabel('Surface Score 2')
ax.set_zlabel('Combined AUC')
ax.set_title('2-Gene Pareto Optimal')