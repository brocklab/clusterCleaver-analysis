# %%
import scanpy as sc
from pathlib import Path

from optimalSeparation import searchOptimal, dataLoading, visualization
from scrna.cluster.main import compute_dimensionality_reductions

'../../data/tem'
# %%
adataPath = Path('../../data/h5ads/temp/231CountsAlra.h5ad')
if adataPath.exists():
    adataImpute = sc.read_h5ad(adataPath)
else:
    print('Loading from csv')
    adata = sc.read_csv('../../data/toSeurat/231CountsAlra.csv')
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adata)
    sc.tl.leiden(adata, resolution=0.05)
adata231 = sc.read_h5ad('../../data/h5ads/final/jostner_mdamb231.h5ad')
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%

sc.pl.umap(adataImpute, color = ['leiden', 'ESAM'])
sc.pl.umap(adata231, color = ['leiden', 'ESAM'])

# %%
emdGenesImpute = searchOptimal.searchExpressionDist(adataImpute, 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)

emdGenes = searchOptimal.searchExpressionDist(adata231, 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'remove0',
                                                 scale = False).reset_index(drop = True)
# %%
emdGenesImpute.head(20)
# %%
emdGenes.head(20)
# %%
sc.pl.umap(adataImpute, color = 'RGS2', title = 'Imputed RGS2 Expression')
sc.pl.umap(adata231, color = 'RGS2', title = f'Non-Imputed RGS2 Expression')
# %%
gene = 'RGS2'
visualization.plotHists(adata231,       gene = gene)
visualization.plotHists(adataImpute,    gene = gene)

# %%
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.write_h5ad('../../data/h5ads/temp/231CountsAlra.h5ad')
