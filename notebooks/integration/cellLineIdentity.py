# %%
import scanpy as sc
# %%
adata = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adata.obs['cmo'] = adata.obs.index.str.split('-').str[-1]
# %%
genesmdamb231 = ["C4BPB",   "ADGRF5",  "GPAT2",   "PHYHD1",  "ACSL5",   "CCRL2",   "SPINK4",  "ADGRF1",  "SERPINB9"]
genesmdamb453 = ["LRRC31","DACH1", "LHX1",  "ADCY2", "NUPR2", "GRM4",  "ABCA12"]
genesbt474 = ["MUC6",   "SCGB1D2","ALOX15B","PGR",    "FAM3B",  "SCGB2A2","CEACAM6","IKZF3"]
genesmdamb436 = ["C1QTNF4","IL1A",   "PTPRN",  "IL13RA2","NR0B1",  "TBX18",  "DDX43"]
geneshs578t = ["PRSS21","MYOCD", "ITGA8", "NTM",   "LMOD1"]
geneshcc38 = ["KLK6",  "ABO",   "MMP20", "CT45A5","RSPO4", "OR3A2"]

genes = genesmdamb231 + genesmdamb453 + genesbt474 + genesmdamb436 + geneshs578t + geneshcc38

markers = {'mdamb231':  genesmdamb231,
           'mdamb453':  genesmdamb453,
           'bt474':     genesbt474,
           'mdamb436':  genesmdamb436,
           'hs578t':    geneshs578t,
           'hcc38':     geneshcc38}
adataSub = adata[:, genes]
# %%
sc.pl.matrixplot(adata, 
                 var_names = markers, 
                 groupby = 'sample', 
                 dendrogram=True, 
                 standard_scale='group')
# %%
import pandas as pd
import os
# %%
files = os.listdir('.')
dm = {}
for file in files:
    if file.endswith('.csv'):
        cell = file.split('_')[2]
        dm[cell] = pd.read_csv(file, index_col = 0)
        dm[cell] = dm[cell]['x'].tolist()
# %%
sc.pl.matrixplot(adata, 
                 var_names = dm, 
                 groupby = 'cmo', 
                 dendrogram=True, 
                 standard_scale='var')
# %%
# Replacing CMOS 
sampleDict = {'1': 'bt474',
              '2': 'mdamb453',
              '3': 'hcc38',
              '4': 'mdamb231',
              '5': 'hs578t',
              '6': 'mdamb436'}
adata.obs['sample'] = adata.obs.index.str.split('-').str[1].to_series().astype('string').replace(sampleDict).tolist()
adata.write_h5ad('../../data/h5ads/jostner-processed.h5ad')