# %%
import scanpy as sc
import pandas as pd
import numpy as np

from optimalSeparation import searchOptimal, dataLoading, visualization
from scrna.cluster.main import compute_dimensionality_reductions

# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')

adatas = dataLoading.processFullAnndata(adataFull)
# %%
adata = adatas['bt474']
# %%
geneFileLoc = '/stor/work/Brock/Tyler/BT474Project/data/Macosko_cell_cycle_genes.txt'
cc_genes = pd.read_table(geneFileLoc, delimiter='\t')
s_genes = cc_genes['S'].dropna()
g2m_genes = cc_genes['G2.M'].dropna()

s_genes_mm = [gene.lower() for gene in s_genes]
g2m_genes_mm = [gene.lower() for gene in g2m_genes]

s_genes_mm_ens = adata.var_names[np.in1d([i.lower() for i in adata.var_names], s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d([i.lower() for i in adata.var_names], g2m_genes_mm)]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
adata.obs['cc_difference'] = adata.obs['S_score'] - adata.obs['G2M_score']
sc.pp.regress_out(adata, 'cc_difference')
# %%
sc.pl.umap(adata, color = 'cc_difference')
# %%
adata = adatas['bt474']
cell_cycle_genes = [x.strip() for x in open('../../data/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# %%
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
# %%
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
sc.pp.scale(adata)
# %%
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
# %%
compute_dimensionality_reductions(adata)
sc.pl.umap(adata)
# %%
sc.tl.leiden(adata, resolution = 0.2)
sc.pl.pca(adata, color = ['leiden'])

# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')

emdGenes = searchOptimal.searchExpressionDist1D(adata, surfaceGenes['gene'], modifier = None)

# %%
adata = adataFull.copy()
cell_cycle_genes = [x.strip() for x in open('../../data/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
sc.pp.scale(adata)
adata_cc_genes = adata[:, cell_cycle_genes]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
compute_dimensionality_reductions(adata)
sc.pl.umap(adata, color = 'sample')
# %%
sc.pl.pca(adata, color = ['sample', 'leiden'])
# %%
adata2 = adata[(adata.obs['sample'] == 'mdamb231') & (adata.obs['scDblFinder_class'] == 'singlet')]
compute_dimensionality_reductions(adata2)
# %%
sc.pl.umap(adata2, color = 'ESAM')
# %%
adata = adata[adata.obs['scDblFinder_class'] == 'singlet']
compute_dimensionality_reductions(adata)
sc.pl.umap(adata, color = 'sample')

# %%

# %%
samples = adata.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adata[adata.obs['sample'].isin([sample])]
    adataSub = adataSub[adataSub.obs['scDblFinder_class'] == 'singlet']
    # sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adataSub)
    # sc.pl.umap(adataSub, title = sample)
    adatas[sample] = adataSub
# %%
cellLineRes = {}
for cellLine, adataSub in adatas.items():
    leidenResolution = 2
    nLeiden = 5
    c = 1
    if cellLine == 'hcc38':
        initLeiden = 3
    else:
        initLeiden = 2
    while nLeiden != initLeiden:
        leidenResolution /= 1.3
        sc.tl.leiden(adataSub, resolution= leidenResolution)
        nLeiden = len(adataSub.obs['leiden'].unique())
        c += 1
        if c > 20:
            leidenResolution = 0
            break
    cellLineRes[cellLine] = leidenResolution
# %%
import seaborn as sns
import matplotlib.pyplot as plt

leidens = []
for cellLine, adata in adatas.items():
    leidens.append(adatas[cellLine].obs['leiden'])
    
if 'leiden' not in adata.obs.columns:
    adata.obs = adata.obs.join(pd.concat(leidens))
    
adata.obs['cellCluster'] = adata.obs['sample'].str[:]+'-'+adata.obs['leiden'].str[:]

pcaVals = {}
for cellCluster in adata.obs['cellCluster'].unique():
    isCellClust = adata.obs['cellCluster'] == cellCluster
    pcaVals[cellCluster] = adata.obsm['X_pca'][isCellClust,:].mean(axis = 0)
    # pcaVals[cellCluster] = np.array(adataFull.X[isCellClust,:].mean(axis = 0))[0]

sc.pl.pca(adata, color = 'cellCluster')

pcaVals = pd.DataFrame(pcaVals)
pcaCorr = pcaVals.corr('pearson')

sns.heatmap(pcaCorr, annot = True, fmt = '.3').set_title('Pearson Correlation')
plt.show()

# %%
