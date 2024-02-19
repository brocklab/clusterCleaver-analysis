# %%
import scanpy as sc
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# %%
adata = sc.read_h5ad('../../data/h5ads/BT474LineageAssigned.h5ad')

adata.obs['sample'] = adata.obs['sample'].str.replace('\d+','', regex = True)

adata = adata[adata.obs['sample'].isin(['PreTreat', 'LP', 'LPD'])]
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

sampleColorDict = {
                 'PreTreat' : '#7F7F01',
                 'D'            : '#00A08B',
                 'DLP'          : '#F032E6',
                 'LP'           : '#FB0D0D',
                 'LPD'          : '#002bff'
}
fig = plt.figure(figsize = (16, 8))
gs = gridspec.GridSpec(1, 2)
gs.update(wspace = 0.5, hspace = 0.5)
ax1 = plt.subplot(gs[1])

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])



umappts = adata.obsm['X_umap'].copy()
groupName = 'unsupGroupPred'
identity = adata.obs[groupName]
groupNames = list(adata.obs[groupName].unique())
groupNames.sort()
for group in groupNames:
        isSample = identity == group
        X = umappts[isSample, 0].copy()
        Y = umappts[isSample, 1].copy()

        ax1.scatter(X, Y, s = 2, label = group, c = unsupColorDict[group])

lgnd = ax1.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()




ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)

ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)
ax1.set_title('Unsupervised Clusters', fontsize = 20)

ax2 = plt.subplot(gs[0])

ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('UMAP 1', loc = 'left')
ax2.set_ylabel('UMAP 2', loc = 'bottom')


umappts = adata.obsm['X_umap'].copy()
groupName = 'sample'
identity = adata.obs[groupName]
groupNames = list(adata.obs[groupName].unique())
groupNames = ['PreTreat', 'LP', 'LPD']
for group in groupNames:
        isSample = identity == group
        X = umappts[isSample, 0].copy()
        Y = umappts[isSample, 1].copy()

        ax2.scatter(X, Y, 
                    s = 2, 
                    label = group, 
                    c = sampleColorDict[group], 
                    alpha = 0.5)

lgnd = ax2.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
    handle.set_alpha([1])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax2.get_xlim() 
ymin, ymax = ax2.get_ylim()


ax2.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax2.arrow(xmin, ymin, 0., 2, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax2.xaxis.set_label_coords(0.05, 0.025)
ax2.yaxis.set_label_coords(0.025, 0.05)

ax2.xaxis.label.set_fontsize(15)
ax2.yaxis.label.set_fontsize(15)
ax2.set_title('Samples', fontsize = 20)

fig.savefig('../../figures/bt474ClusterUMAP.png', dpi = 500, bbox_inches = 'tight')
# %%
