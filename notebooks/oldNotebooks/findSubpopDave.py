# %%
import scanpy as sc
import os
from pathlib import Path
from tqdm import tqdm
import pickle
import pandas as pd
import anndata

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal
from optimalSeparation import dataLoading

# %%
dataDir = Path('../data/h5ads')
pickledDave = dataDir.parent / 'dave.pickle'
if pickledDave.exists():
    h5s = pickle.load(open('../data/dave.pickle', 'rb'))
else:
    h5s = {}
    for filePath in tqdm(dataDir.iterdir()):
        fileName = str(filePath.name)
        if fileName.startswith('dave') and fileName.endswith('clustered.h5ad'):
            fileName = fileName[4:]
            cellLine = fileName.split('-')[0]
            adata = sc.read_h5ad(filePath)
            adata.var.index = adata.var['gene_ids']
            adata.var.index = adata.var.index.astype(str)
            adata.var_names_make_unique()
            adata.X = adata.X.toarray()

            h5s[cellLine] = adata
    pickle.dump(h5s, open(pickledDave, 'wb'))
# %%
# %%
if 'MDAMB231' in h5s:
    h5s.pop('MDAMB231')
adata = anndata.concat(list(h5s.values()))
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)
sc.pl.umap(adata, color = 'cellLine')
# %%
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# %%
for cellLine, adata in h5s.items():
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    compute_dimensionality_reductions(adata)
# %%
cellLineRes = {}
for cellLine, adata in h5s.items():
    leidenResolution = 0.5
    nLeiden = 5
    c = 1
    while nLeiden != 2:
        sc.tl.leiden(adata, resolution= leidenResolution)
        leidenResolution /= 1.2
        nLeiden = len(adata.obs['leiden'].unique())
        c += 1
        if c > 20:
            leidenResolution = 0
            break
    cellLineRes[cellLine] = leidenResolution

# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('..')
# %%
optimalGenesConcat = []
optimalCombosConcat = []
for cellLine in h5s.keys():
    adata = h5s[cellLine]

    sc.tl.leiden(adata, resolution=cellLineRes[cellLine])
    try:
        optimalGenes = searchOptimal.searchGeneSeparation(adata, surfaceGenes['gene'])
        optimalCombos = searchOptimal.searchSeparation2(adata, optimalGenes, nGenes = 75)

        optimalGenes['cellLine'] = cellLine
        optimalCombos['cellLine'] = cellLine

        optimalGenesConcat.append(optimalGenes)
        optimalCombosConcat.append(optimalCombos)
    finally:
        print(f'{cellLine} failed to search for optimal genes')
        continue
optimalGenesConcat = pd.concat(optimalGenesConcat)
optimalCombosConcat = pd.concat(optimalCombosConcat)

optimalGenesConcat.to_csv('../data/surfacePreds/daveSubpopulationOptimalGenes.csv')
optimalCombosConcat.to_csv('../data/surfacePreds/daveSubpopulationOptimalGeneCombos.csv')

# %%
optimalGenes = pd.read_csv('../data/surfacePreds/daveSubpopulationOptimalGenes.csv', index_col=0)
# %%
for cellLine in optimalGenes['cellLine'].unique():
    adata = h5s[cellLine]

    optimalGeneCellLine = optimalGenes.loc[optimalGenes['cellLine'].isin([cellLine])].sort_values(by = 'auc', ascending=False)
    topGene = optimalGeneCellLine['gene1'].iloc[0]
    sc.pl.umap(adata, color = ['leiden', topGene], use_raw = False)
# %%
