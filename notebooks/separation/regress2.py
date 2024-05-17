# %% 
import scanpy as sc
import pandas as pd
import numpy as np
# %%
adata = sc.read_csv('../../data/h5ads/nestorawa_forcellcycle_expressionMatrix.txt', delimiter='\t').T

# %%
cell_cycle_genes = [x.strip() for x in open('../../data/regev_lab_cell_cycle_genes.txt')]

# %%
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
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
print(np.min(adata_cc_genes.X))
# %%
s_genes_mm = [gene.lower() for gene in s_genes]
g2m_genes_mm = [gene.lower() for gene in g2m_genes]

s_genes_mm_ens = adata.var_names[np.in1d([i.lower() for i in adata.var_names], s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d([i.lower() for i in adata.var_names], g2m_genes_mm)]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
adata.obs['cc_difference'] = adata.obs['S_score'] - adata.obs['G2M_score']
sc.pp.regress_out(adata, 'cc_difference')
# %%
print(np.min(adata.X))
