# %%
import scanpy as sc
from pathlib import Path
import matplotlib.pyplot as plt

from optimalSeparation import searchOptimal, dataLoading, visualization
from scrna.cluster.main import compute_dimensionality_reductions
# %%
alraPath = Path('../../data/h5ads/temp/231CountsAlra.h5ad')
if alraPath.exists():
    adataImpute = sc.read_h5ad(alraPath)
else:
    print('Loading from csv')
    adataImpute = sc.read_csv('../../data/toSeurat/231CountsAlra.csv')
    sc.pp.highly_variable_genes(adataImpute, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataImpute)
    sc.tl.leiden(adataImpute, resolution=0.05)
    adataImpute.write_h5ad('../../data/h5ads/temp/231CountsAlra.h5ad')
seuratPath = Path('../../data/h5ads/temp/231CountsSeurat.h5ad')

if seuratPath.exists():
    adataSeurat = sc.read_h5ad(seuratPath)
else:
    print('Loading from csv')
    adataSeurat = sc.read_csv('../../data/toSeurat/231CountsSeurat.csv')
    sc.pp.highly_variable_genes(adataSeurat, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataSeurat)
    sc.tl.leiden(adataSeurat, resolution=0.05)
    adataSeurat.write_h5ad('../../data/h5ads/temp/231CountsSeurat.h5ad')

surfaceGenes = dataLoading.cleanSurfaceGenes('../..')
# %%

sc.pl.umap(adataImpute, color = ['leiden', 'ESAM'])
sc.pl.umap(adataSeurat, color = ['leiden', 'ESAM'])

# %%
emdGenesImpute = searchOptimal.searchExpressionDist(adataImpute, 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = None,
                                                 scale = False).reset_index(drop = True)

emdGenes = searchOptimal.searchExpressionDist(adataSeurat, 
                                                 surfaceGenes['gene'],
                                                 metric = 'EMD', 
                                                 modifier = 'remove0',
                                                 scale = False).reset_index(drop = True)
# %%
emdGenesImpute.head(20)
# %%
emdGenes.head(20)

# %%
gene = 'RASD1'
visualization.plotHists(adataSeurat, gene = gene)
visualization.plotHists(adataImpute, gene = gene)

# %%
fig1 = sc.pl.umap(adataImpute, color = 'RGS2', title = 'Imputed RGS2 Expression', return_fig = True)
fig2 = sc.pl.umap(adataSeurat, color = 'RGS2', title = f'Non-Imputed RGS2 Expression', return_fig = True)
fig1.savefig('../../figures/02-08-24_Figures/rgs2UMAP.png', dpi = 500)
fig2.savefig('../../figures/02-08-24_Figures/rgs2UMAPImputed.png', dpi = 500)
# %%
gene = 'RGS2'
visualization.plotHists(adataSeurat, gene = gene, saveFig = '../../figures/02-08-24_Figures/rgs2Expression.png')
visualization.plotHists(adataImpute, gene = gene, saveFig = '../../figures/02-08-24_Figures/rgsImputed2Expression.png')