# %%
import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
print('hi')
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
# %%
samples = adataFull.obs['sample'].unique()
adatas = {}
for sample in samples:
    adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
    adataSub = adataSub[adataSub.obs['scDblFinder_class'] == 'singlet']
    sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
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
                 'hs578t': 'Hs 578T',
                 'mdamb436': 'MDA-MB-436'}
fig, axs = plt.subplots(3, 2, figsize = (7.5,10))
axs = axs.ravel()
for i, (cellLine, adata) in enumerate(adatas.items()):
    print(cellLine)
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
        ax.arrow(xmin, ymin, 3, 0., fc='k', ec='k', lw = 1, 
            head_width=0.25, head_length=0.25, overhang = 0.3, 
            length_includes_head= False, clip_on = False) 


        ax.arrow(xmin, ymin, 0., 2, fc='k', ec='k', lw = 1, 
                head_width=0.25, head_length=0.25, overhang = 0.3, 
                length_includes_head= False, clip_on = False) 

        ax.xaxis.set_label_coords(0.05, 0.025)
        ax.yaxis.set_label_coords(0.025, 0.05)

        ax.xaxis.label.set_fontsize(10)
        ax.yaxis.label.set_fontsize(10)
fig.savefig('../../figures/umapAllLinesClustered.png', dpi = 500)

# %%
