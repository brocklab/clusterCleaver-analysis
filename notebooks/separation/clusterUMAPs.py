# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed-regressed.h5ad')
# %%
adataFull = adataFull[adataFull.obs['scDblFinder_class'] == 'singlet']
# adataFull = adataFull[adataFull.obs['sample'].isin(['mdamb231', 'mdamb436'])]

# sc.pp.highly_variable_genes(adataFull, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adataFull)
sc.pl.umap(adataFull, color = 'sample')
# %%
sns.reset_defaults()

cellDict = {
    'bt474': 'BT474',
    # 'hcc38': 'HCC38',
    'hs578t': 'Hs578T',
    'mdamb231': 'MDA-MB-231',
    'mdamb436': 'MDA-MB-436',
    'mdamb453': 'MDA-MB-453'
}

fig, ax = plt.subplots(1, 1, figsize=(10,10))
umapPts = adataFull.obsm['X_umap']
for sample in cellDict.keys():
    inSample = adataFull.obs['sample'] == sample
    x = umapPts[inSample,0]
    y = umapPts[inSample, 1]

    plt.scatter(x,y, s = 5, label = cellDict[sample])
ax.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
ax.set_xlabel('UMAP 1', loc = 'left', fontsize = 18)
ax.set_ylabel('UMAP 2', loc = 'bottom', fontsize = 18)

ax.arrow(xmin, ymin, 2, 0., fc='k', ec='k', lw = 1.5, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax.arrow(xmin, ymin, 0., 1.4, fc='k', ec='k', lw = 1.5, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax.xaxis.set_label_coords(0.05, 0.025)
ax.yaxis.set_label_coords(0.025, 0.05)

lgnd = ax.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

fig.savefig('../../figures/final/umapConcat.png', dpi = 500, bbox_inches = 'tight')
plt.show()
# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB', '#E87941', '#FE64A3']
sns.set_palette(sns.color_palette(colors))

samples = adataFull.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
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
    print(adataSub.obs['leiden'].unique())
    print(f'Cell Line: {cellLine} \t Resolution: {leidenResolution} \t nLeiden: {nLeiden}')
# %%
cellLineNames = {'bt474': 'BT474',
                 'mdamb231': 'MDA-MB-231',
                 'mdamb453': 'MDA-MB-453',
                 'hcc38'   : 'HCC38',
                 'hs578t': 'Hs578T',
                 'mdamb436': 'MDA-MB-436'}
fig, axs = plt.subplots(3, 2, figsize = (7.5,10/3))
axs = axs.ravel()
for i, (cellLine, adata) in enumerate(adatas.items()):
    ax = axs[i]
    ax.spines[["top", "right", 'left', 'bottom']].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])

    if i == 0:
        ax.set_xlabel('UMAP 1', loc = 'left')
        ax.set_ylabel('UMAP 2', loc = 'bottom')
    ax.set_title(cellLineNames[cellLine])

    umappts = adata.obsm['X_umap']
    identity = adata.obs['leiden']

    for cat in identity.unique():
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]

        # if cat == 'LPDOther':
        #     alpha = 1
        # else:
        #     alpha = 1
        ax.scatter(X, Y, s = 2)

    # for handle in lgnd.legend_handles:
    #     handle.set_sizes([80])
    # lgnd.get_frame().set_linewidth(0.0)
    if i == 0:
        xmin, xmax = ax.get_xlim() 
        ymin, ymax = ax.get_ylim()
        ax.arrow(xmin, ymin, 2, 0., fc='k', ec='k', lw = 1, 
            head_width=0.25, head_length=0.25, overhang = 0.3, 
            length_includes_head= False, clip_on = False) 


        ax.arrow(xmin, ymin, 0., 2, fc='k', ec='k', lw = 1, 
                head_width=0.25, head_length=0.25, overhang = 0.3, 
                length_includes_head= False, clip_on = False) 

        ax.xaxis.set_label_coords(0.05, 0.025)
        ax.yaxis.set_label_coords(0.025, 0.05)

        ax.xaxis.label.set_fontsize(10)
        ax.yaxis.label.set_fontsize(10)

# fig.delaxes(axs[-1])
plt.show()
fig.savefig('../../figures/final/umapAllLinesClustered.png', dpi = 500)

# %% Get PCC
leidens = []
for cellLine, adata in adatas.items():
    leidens.append(adatas[cellLine].obs['leiden'])
    
if 'leiden' not in adataFull.obs.columns:
    adataFull.obs = adataFull.obs.join(pd.concat(leidens))
    
adataFull.obs['cellCluster'] = adataFull.obs['sample'].str[:]+'-'+adataFull.obs['leiden'].str[:]
pcaVals = {}
for cellCluster in adataFull.obs['cellCluster'].unique():
    isCellClust = adataFull.obs['cellCluster'] == cellCluster
    pcaVals[cellCluster] = adataFull.obsm['X_pca'][isCellClust,:].mean(axis = 0)
    # pcaVals[cellCluster] = np.array(adataFull.X[isCellClust,:].mean(axis = 0))[0]

sc.pl.pca(adataFull, color = 'cellCluster')

pcaVals = pd.DataFrame(pcaVals)
pcaCorr = pcaVals.corr('pearson')

import seaborn as sns
sns.heatmap(pcaCorr, annot = True, fmt = '.3').set_title('Pearson Correlation')
plt.show()
# %%
sc.pl.pca(adataFull, color = 'cellCluster', title = 'cellLine - cluster #')
# %%
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr
for sample, adata in adatas.items():
    is0 = adata.obs['leiden'] == '0'
    is1 = adata.obs['leiden'] == '1'

    v0 = adata.obsm['X_pca'][is0,:].mean(axis = 0)
    v1 = adata.obsm['X_pca'][is1,:].mean(axis = 0)

    # v0 = np.array(adata.X[is0,:].mean(axis = 0))[0]
    # v1 = np.array(adata.X[is1,:].mean(axis = 0))[0]
    # corr = pearsonr(v0, v1).statistic
    corr = cosine(v0, v1)
    print(f'{sample} correlation: {corr:0.2f}')

    plt.figure()
    plt.scatter(v0, v1)

    plt.figure()
    sc.pl.pca(adata, color = 'leiden')
    
    plt.show()