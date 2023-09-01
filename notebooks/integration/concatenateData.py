# %% 
import scanpy as sc
import pandas as pd
import numpy as np

import scvi
import anndata
from pathlib import Path
# %%
baseDir = Path('../data/h5ads')
# %%
# %%
cellLinesKeep = ['NCIH1299_LUNG', 'MDAMB436_BREAST', 'BT549_BREAST', 'BT474_BREAST']
kinker = sc.read_h5ad(baseDir / 'kinker-postqc-normalized-clustered.h5ad')
kinker = kinker[kinker.obs['Cell_line'].isin(cellLinesKeep)]
kinker.obs['cellLine'] = kinker.obs['Cell_line'].str.split('_').str[0]
kinker.obs['cellLine'] = kinker.obs['cellLine'].replace('NCIH1299', 'H1299')
kinker.obs['batch'] = 'kinker'
kinker.var_names_make_unique()
adata = kinker

# del kinker
# %% Dave
dave549 = sc.read_h5ad(baseDir / 'daveBT549-postqc-normalized-clustered.h5ad')
dave231 = sc.read_h5ad(baseDir / 'daveMDAMB231-postqc-normalized-clustered.h5ad')
dave436 = sc.read_h5ad(baseDir / 'daveMDAMB436-postqc-normalized-clustered.h5ad')

dave549.obs['batch'] = 'dave549'
dave231.obs['batch'] = 'dave231'
dave436.obs['batch'] = 'dave436'

# Set to symbol ids
dave549.var = dave549.var.reset_index().set_index('gene_ids').rename_axis(index=None).astype(str)
dave231.var = dave231.var.reset_index().set_index('gene_ids').rename_axis(index=None).astype(str)
dave436.var = dave436.var.reset_index().set_index('gene_ids').rename_axis(index=None).astype(str)

dave549.var.index = dave549.var.index.astype(str)
dave231.var.index = dave231.var.index.astype(str)
dave436.var.index = dave436.var.index.astype(str)

dave549.var_names_make_unique()
dave231.var_names_make_unique()
dave436.var_names_make_unique()

adata = anndata.concat([adata, dave231, dave436, dave549])

# del dave231
# del dave436
# del dave549
# %%
adataKim = sc.read_h5ad('../data/h5ads/kim-postqc-normalized-clustered.h5ad')
adataKim = adataKim[adataKim.obs['cellLine'].isin(['H1299']),:]
adataKim.obs['batch'] = 'kim'
adata = anndata.concat([adata, adataKim])

# del adataKim
# %%
cellLinesKeep = ['BT474', 'MDAMB436', 'BT549']

adataGambardella = sc.read_h5ad('../data/h5ads/gambardella-postqc-normalized-clustered.h5ad')
adataGambardella = adataGambardella[adataGambardella.obs['cellLine'].isin(cellLinesKeep)]

adataGambardella.obs['batch'] = 'gambardella'

adata = anndata.concat([adata, adataGambardella])
# %%
adataDing = sc.read_h5ad('../data/h5ads/ding-postqc-normalized-clustered.h5ad')
adataDing = adataDing[adataDing.obs['samples'].isin(['5'])]
adataDing.obs['cellLine'] = 'HCC1806'
adataDing.obs['batch'] = 'ding'

adata = anndata.concat([adata, adataDing])
# %%
del dave231
del dave436
del dave549
del adataKim
del kinker
# %%
adataBT474  = sc.read_h5ad('../data/h5ads/BT474LineageAssigned.h5ad')
adataBT474 = adataBT474[adataBT474.obs['sample'].str.startswith('Pre')]
adataBT474.obs['cellLine'] = 'BT474'
adataBT474.obs['batch'] = 'jost'
adata = anndata.concat([adata, adataBT474])

# %%
adataHCC1806 = sc.read_h5ad('../data/h5ads/HCC1806LineageAssigned.h5ad')
adataHCC1806 = adataHCC1806[adataHCC1806.obs['sample'].isin(['PT'])]
adataHCC1806.obs['cellLine'] = 'HCC1806'
adataHCC1806.obs['batch'] = 'desan'
adata = anndata.concat([adata, adataHCC1806])

# %%
adata231_1kb3 = sc.read_h5ad('../data/h5ads/231-1KB3-final.h5ad')
adata231_1kb3 = adata231_1kb3[adata231_1kb3.obs['sample'].isin(['PT'])]
adata231_1kb3.obs['cellLine'] = 'MDAMB231'
adata231_1kb3.obs['batch'] = '1kb3'
adata = anndata.concat([adata, adata231_1kb3])
# %%
adata231_AA115 = sc.read_h5ad('../data/h5ads/231-AA115-final.h5ad')
adata231_AA115 = adata231_AA115[adata231_AA115.obs['sample'].isin(['PT'])]
adata231_AA115.obs['cellLine'] = 'MDAMB231'
adata231_AA115.obs['batch'] = 'AA115'
adata = anndata.concat([adata, adata231_AA115])
# %%
adata.write_h5ad('../data/h5ads/adataConcatInner.h5ad')
# %%
adata = sc.read_h5ad('../../data/h5ads/adataConcatInner.h5ad')
# %%
