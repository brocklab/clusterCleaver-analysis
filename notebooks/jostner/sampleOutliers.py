# %%
# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal
from optimalSeparation import dataLoading
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
# %%
sc.tl.rank_genes_groups(adataFull, groupby = 'sample', key = 'rank_genes_groups')
# %%
topGenes = []
for sample in adataFull.obs['sample'].unique():
    rankDf = sc.get.rank_genes_groups_df(adataFull, group=sample).sort_values(
        by = 'logfoldchanges', ascending=True)
    genes = rankDf['names'].iloc[0:3].tolist()
    topGenes += genes
# %%
sc.pl.matrixplot(adataFull, var_names = topGenes, groupby='sample', dendrogram = True)
# %%
sc.pl.umap(adataFull, color = 'POSTN')
# %%
