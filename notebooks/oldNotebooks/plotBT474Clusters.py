# %%
import scanpy as sc
import pandas as pd
# %%
adata = sc.read_h5ad('../../data/h5ads/BT474LineageAssigned.h5ad')
# %%
assignedClustersPath = '../../../BT474Project/data/unsupGroupPred.csv'
assignedClusters = pd.read_csv(assignedClustersPath, index_col=0)

assignedClusters.set_index('cellProper', inplace=True)
assignedClusters = assignedClusters.rename_axis(index = None)
# %%
if 'unsupGroupPred' not in adata.obs.columns:
    adata.obs = adata.obs.merge(assignedClusters, left_index = True, right_index = True)
adata.obs['unsupGroupPred'] = adata.obs['unsupGroupPred'].astype('string')

# %%
unsupColorDict = {
                '1': '#db5f57',
                '2': '#d3db57',
                '3': '#57db5f',
                '4': '#57d2db',
                '5': '#5457db'
}

