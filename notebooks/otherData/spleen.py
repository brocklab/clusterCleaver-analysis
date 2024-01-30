# %%
import pandas as pd
import anndata
import scanpy as sc
# %%
adata = sc.read_h5ad('../../data/h5ads/hanSpleen-postqc-normalized-clustered.h5ad')
# %%
cellTypes = pd.read_csv('../../data/MCA_CellAssignments.csv', index_col=0)

cellTypes = cellTypes.loc[cellTypes['Tissue'] == 'Spleen']
# %%
cellTypesJoin = cellTypes[['Annotation', 'Cell.name']]
cellTypesJoin = cellTypesJoin.set_index('Cell.name', drop=True)
cellTypesJoin.index.name = None
adata.obs = adata.obs.join(cellTypesJoin)

# %%
adata.obs['isSplenic'] = 0
adata.obs.loc[adata.obs['Annotation'] == 'Marginal zone B cell(Spleen)', 'isSplenic'] = 1

sc.pl.umap(adata, color = 'isSplenic')
# %%
from optimalSeparation import searchOptimal, dataLoading, visualization
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%
emdGenesNew = searchOptimal.searchExpressionDist(adata, adata.var.index, label = 'isSplenic', modifier = 'remove0')
# %%
visualization.plotHists(adata, gene = 'Ly6d', colorCol = 'isSplenic')
# %%
visualization.plotHists(adata, gene = 'Cd74', colorCol = 'isSplenic')
# %%
topGenes = emdGenesNew['genes'][0:75].tolist()
# %%
allEMDCombos = searchOptimal.searchExpressionDist(adata, adata.var.index, nGenes = 2, topGenes = topGenes, modifier = None)
# %%
allEMDCombos = searchOptimal.searchExpressionDist(adata, adata.var.index, nGenes = 2, topGenes = ['Cd19, Ly6d', ''], modifier = None)

# %%
visualization.plotHists(adata, gene = 'Cd79b', colorCol = 'isSplenic')
# %%
visualization.plotExpression(adata, genes = ['Ly6d', 'Cd19'], colorCol = 'isSplenic')
# %%
visualization.plotExpression(adata, genes = ['Ly6d', 'Cd74'], colorCol = 'isSplenic')
