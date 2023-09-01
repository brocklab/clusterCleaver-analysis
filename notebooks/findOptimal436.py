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
adata = sc.read_h5ad('../data/h5ads/daveMDAMB436-postqc-normalized-clustered.h5ad')
adata = adata[adata.obs['cellLine'].str.startswith('MDAMB436')]
adata.var.index = adata.var['gene_ids']
adata.var.index = adata.var.index.astype(str)
adata.var_names_make_unique()
adata.X = adata.X.toarray()
# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
sc.tl.leiden(adata, resolution = 0.1)
sc.pl.umap(adata, color = 'leiden')
# %%
surfaceGenes = cleanSurfaceGenes('..')
# %%
optimalGenesDave = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
optimalGenesSurfaceDave = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenesDave, nGenes = 1)
optimalGenesParetoDave = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurfaceDave)

optimalGenes2Dave = searchOptimal.searchGeneSeparation(adata, optimalGenesDave['gene1'][0:100], nGenes = 2, nCombos = 10000)
optimalGenesSurface2Dave = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes2Dave, nGenes = 2)
optimalGenesPareto2Dave = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2Dave, nGenes = 2)
# %%
visualization.plotParetoOptimal(optimalGenesSurfaceDave, optimalGenesParetoDave, nGenes = 1)
plt.title('436 Pareto Optimal Genes')

visualization.plotParetoOptimal(optimalGenesSurface2Dave, optimalGenesPareto2Dave, nGenes = 2)
plt.title('436 Pareto Optimal Genes')

# %%
visualization.plotHists(adata, 'SLC3A2')
# %%
plt.scatter(adata[:, 'SLC3A2'].X, adata[:, 'ATP2B1'].X, c = adata.obs['leiden'].astype(int))
# %%
plt.scatter(np.linspace(0, 1, adata.shape[0]), adata[:, 'SLC3A2'].X, c = adata.obs['leiden'].astype(int))
# %%
sc.pl.umap(adata, color = 'SLC3A2', use_raw = False)
# %%
gambardella = sc.read_h5ad('../data/h5ads/gambardella-postqc-normalized-clustered.h5ad')
gambardella = gambardella[gambardella.obs['cellLine'].str.startswith('MDAMB436')]
# %%
sc.pp.highly_variable_genes(gambardella, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(gambardella)
sc.tl.leiden(gambardella, resolution = 0.1)
sc.pl.umap(gambardella, color = 'leiden')
gambardella.X = gambardella.X.toarray()
# %%
optimalGenesGambardella = searchOptimal.searchGeneSeparation(gambardella, surfaceGenes['gene'])
optimalGenesSurfaceGambardella = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenesGambardella, nGenes = 1)
optimalGenesParetoGambardella = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurfaceGambardella)

optimalGenes2Gambardella = searchOptimal.searchGeneSeparation(gambardella, optimalGenesGambardella['gene1'][0:100], nGenes = 2, nCombos = 10000)
optimalGenesSurface2Gambardella = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes2Gambardella, nGenes = 2)
optimalGenesPareto2Gambardella = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2Gambardella, nGenes = 2)
# %%
visualization.plotParetoOptimal(optimalGenesSurfaceGambardella, optimalGenesParetoGambardella, nGenes = 1)
plt.title('436 (Gambardella) Pareto Optimal Genes')

visualization.plotParetoOptimal(optimalGenesSurface2Gambardella, optimalGenesPareto2Gambardella, nGenes = 2)
plt.title('436 (Gambardella) Pareto Optimal Genes')
# %%
visualization.plotHists(gambardella, 'SLC3A2')
# %%
plt.scatter(gambardella[:, 'SLC3A2'].X.toarray(), gambardella[:, 'EEF1A1'].X.toarray(), c = gambardella.obs['leiden'].astype(int))
# %% 
kinker = sc.read_h5ad('../data/h5ads/kinker-postqc-normalized-clustered.h5ad')
kinker.obs['cellLine'] = kinker.obs['Cell_line']
kinker.obs['cellLine'] = kinker.obs['cellLine'].replace('MDAMB436_BREAST', 'MDAMB436')
kinker = kinker[kinker.obs['cellLine'].str.startswith('MDAMB436')]
kinker.obs['batch'] = 'kinker'
# %%
sc.pp.highly_variable_genes(kinker, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(kinker)
sc.tl.leiden(kinker, resolution = 0.1)
sc.pl.umap(kinker, color = 'leiden')
kinker.X = kinker.X.toarray()
# %%
optimalGenesKinker = searchOptimal.searchGeneSeparation(kinker, surfaceGenes['gene'])
optimalGenesSurfaceKinker = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenesKinker, nGenes = 1)
optimalGenesParetoKinker = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurfaceKinker)

optimalGenes2Kinker = searchOptimal.searchGeneSeparation(kinker, optimalGenesKinker['gene1'][0:100], nGenes = 2, nCombos = 10000)
optimalGenesSurface2Kinker = searchOptimal.mergeSurfaceScores(surfaceGenes, optimalGenes2Kinker, nGenes = 2)
optimalGenesPareto2Kinker = searchOptimal.findOptimalSurfaceMarkers(optimalGenesSurface2Kinker, nGenes = 2)
# %% Look for shared genes
for nGenes in range(10, 100):
    allOptim = list(optimalGenesDave['gene1'][0:nGenes]) + list(optimalGenesGambardella['gene1'][0:nGenes]) + list(optimalGenesKinker['gene1'][0:nGenes])

    sharedGenes = set([x for x in allOptim if allOptim.count(x) == 3])

    if len(sharedGenes) > 0:
        print(sharedGenes)
        break

# %%
dave = adata
topDave = optimalGenesDave['gene1'].iloc[0]
topGambardella = optimalGenesGambardella['gene1'].iloc[0]
topKinker = optimalGenesKinker['gene1'].iloc[0]

plt.rcParams["figure.figsize"] = (20,20)
fig, axs = plt.subplots(3,3)
axs = axs.ravel()
idx = 0
for ad in [dave, gambardella, kinker]:
    for gene in [topDave, topGambardella, topKinker]:
        sc.pl.umap(ad, color = gene, ax = axs[idx], show=False, title='', use_raw = False)
        idx += 1
# %%
