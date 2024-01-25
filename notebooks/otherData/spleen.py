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
adata.obs['isCluster'] = 0
adata.obs.loc[adata.obs['Annotation'] == 'Marginal zone B cell(Spleen)', 'isCluster'] = 1

sc.pl.umap(adata, color = 'isCluster')
# %%
from optimalSeparation import searchOptimal, dataLoading, visualization
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%
emdGenesNew = searchOptimal.searchExpressionDist(adata, adata.var.index, label = 'isCluster', modifier = 'remove0')
# %%
visualization.plotHists(adata, gene = 'Ly6d', colorCol = 'isCluster')
# %%
visualization.plotHists(adata, gene = 'Cd74', colorCol = 'isCluster')

# %%
