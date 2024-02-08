# %%
import scanpy as sc
from optimalSeparation import searchOptimal, dataLoading, visualization
from scrna.cluster.main import compute_dimensionality_reductions


# %%
adata = sc.read_csv('../../data/toSeurat/231CountsAlra.csv')
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')

# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
# %%
sc.tl.leiden(adata, resolution=0.05)
sc.pl.umap(adata, color = ['leiden', 'ESAM'])
# %%
emdGenes = searchOptimal.searchExpressionDist(adata, 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)
# %%
visualization.plotHists(adata, gene = 'RGS2')