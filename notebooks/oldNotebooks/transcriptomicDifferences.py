# %%
from scrna.cluster.main import compute_dimensionality_reductions

import scanpy as sc
import plotnine
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
# %%
adata231 = sc.read_h5ad('../../data/h5ads/231-1KB3-20231013.h5ad')
adata231 = adata231[adata231.obs['sample'].isin(['PT', 'C2'])]
adata231.obs['subsample'] = 'C2'
# %%
adataEsam = adata231[adata231.obs['sample'].isin(['PT'])]
sc.pp.highly_variable_genes(adataEsam)
compute_dimensionality_reductions(adataEsam)
# %%
sc.tl.leiden(adataEsam, resolution=0.05)
sc.pl.umap(adataEsam, color = ['leiden', 'ESAM'])
# %%
adata231.obs.loc[adataEsam.obs.loc[adataEsam.obs['leiden'] == '0'].index, 'subsample'] = 'ESAM Neg'
adata231.obs.loc[adataEsam.obs.loc[adataEsam.obs['leiden'] == '1'].index, 'subsample'] = 'ESAM Pos'
# %%
sc.pp.highly_variable_genes(adata231)
compute_dimensionality_reductions(adata231)
# %%
adata436 = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
# %%
adata436 = adata436[adata436.obs['scDblFinder_class'] == 'singlet']
adata436 = adata436[adata436.obs['sample'] == 'mdamb436']
sc.pp.highly_variable_genes(adata436, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata436)
# %%
sc.tl.leiden(adata436, resolution=0.1)
sc.pl.umap(adata436, color = ['leiden', 'BST2'])
adata436.obs['subsample'] = adata436.obs['leiden']
# %%
adata = sc.concat([adata231, adata436])
# %%
# sc.pp.highly_variable_genes(adata)
# compute_dimensionality_reductions(adata)
# %%
sc.pl.umap(adata, color = ['sample', 'subsample'])
# %%
sc.tl.dendrogram(adata, groupby = 'subsample')
sc.pl.dendrogram(adata, groupby = 'subsample')
# %%
sampleDict = { 'C2':    'Treated',
               'LPD7':  'BT474\nLap/Pac/Dox\nTreated',
               'PT':    'Untreated',
               'mdamb436': 'MDA-MB-436'
}

subsampleDict = {'ESAM Neg': '231-Subpop 1',
                 'ESAM Pos': '231-Subpop 2',
                 '0'       : '436-Subpop 1',
                 '1'       : '436-Subpop 2',
                 'LPDOther': 'Other',
                 'Lin1':     'Lineage 1',
                 'treated':  'Treated'}

allLabelDict = subsampleDict.copy()
for k, v in sampleDict.items():
    allLabelDict[k] = v

adata.obs['sampleNames'] = adata.obs['sample'].astype(str).replace(sampleDict).astype('category')
adata.obs['subsampleNames'] = adata.obs['subsample'].astype(str).replace(subsampleDict).astype('category')
sc.pl.umap(adata, color = ['sampleNames', 'subsampleNames'], wspace = .45)
# %%
isC2 =      adata.obs['sample']     == 'C2'
isPT =      adata.obs['sample']     == 'PT'
isESAMPos = adata.obs['subsample']  == 'ESAM Pos'
isESAMNeg = adata.obs['subsample']  == 'ESAM Neg'
isLin1 =    adata.obs['subsample']  == 'Lin1'
isOther =   adata.obs['subsample']  == 'LPDOther'
is0 =       adata.obs['subsample']  == '0'
is1 =       adata.obs['subsample']  == '1'

X   = adata.obsm['X_pca']
# X = adata.X
C2      =   X[isC2, :].mean(axis = 0).tolist()
PT      =   X[isPT, :].mean(axis = 0).tolist()
ESAMPos =   X[isESAMPos, :].mean(axis = 0).tolist()
ESAMNeg =   X[isESAMNeg, :].mean(axis = 0).tolist()
lin1  =     X[isLin1, :].mean(axis = 0).tolist()
lpdOther =  X[isOther, :].mean(axis = 0).tolist()
# %%
# pcaVals = {'C2': C2, 'PT': PT, 'ESAM Pos': ESAMPos, 'ESAM Neg': ESAMNeg, 'Lin1': lin1, 'LPDOther': lpdOther}
# for key in list(pcaVals.keys()):
#     newKey = ''
#     if key in sampleDict.keys():
#         newKey = sampleDict[key]
#     elif key in subsampleDict.keys():
#         newKey = subsampleDict[key]

#     pcaVals[newKey] = pcaVals.pop(key)   
# dfPcaVals = pd.DataFrame(pcaVals)
# pcaCorr = dfPcaVals.corr()

# sns.heatmap(pcaCorr, annot=True, fmt='.3').set_title('Pearson Correlation of\nPopulations')
# %%
# Treated vs untreated MDAMB231s
# ESAM + vs ESAM -
# Lineage 1 vs Other
untreated = 'MDA-MB-231\nUntreated'
treated =   'MDA-MB-231\nDox Treated'
esamPos =   'MDA-MB-231\Subpopulation 1'
esamNeg =   'MDA-MB-231\Subpopulation 2'
lin1 =      'BT474\nLap/Pac/Dox\nLineage 1'
linOther =  'BT474\nLap/Pac/Dox\nOther'

# corrTreat = pcaCorr.loc[treated, untreated]
# corrSubPop = pcaCorr.loc[esamPos, esamNeg]
# corrLin = pcaCorr.loc[lin1, linOther]

# linBar = pd.DataFrame({'Treated vs. Untreated': corrTreat, 
#                        'Subpopulation vs. Subpopulation': corrSubPop,
#                        'Lineage vs. Sample': corrLin}, 
#                        index = [0]).T


# %%
colorDict = {'LPD7':            '#002bff',
            #  'treated':       '#F68511',
             'ESAM Pos':        '#8fc26c',
             'ESAM Neg':        '#e43820',
             'C2':              '#4bc1ee',
             'PT':              '#7d1811',
             '0':               '#2e336a',
             '1':               '#fed504'
             }

# %% Plot samples

umappts = adata.obsm['X_umap'].copy()
identity = adata.obs['sample'].copy()
plt.rcParams.update({'font.size': 18})

# plt.rcParams["axes.spines.right"] = False
# plt.rcParams["axes.spines.top"] = False

# fig, axs = plt.subplots(1, 3, figsize = (20,16))
fig = plt.figure(figsize = (25, 25))
gs = gridspec.GridSpec(4, 4)
gs.update(wspace = 0.5, hspace = 0.5)
ax1 = plt.subplot(gs[:2, :2])
# ax1 = axs[0]

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('UMAP 1', loc = 'left')
ax1.set_ylabel('UMAP 2', loc = 'bottom')
ax1.set_title('MDA-MB-231 Samples')

for cat in identity.unique():
    if cat in ['C2', 'PT']:
        isSample = identity == cat
        X = umappts[isSample, 0].copy()
        Y = umappts[isSample, 1].copy()

        ax1.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])
lgnd = ax1.legend(prop=dict(size=15), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()


ax1.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 3, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)

ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)

identity = adata.obs['subsample']
# ax2 = axs[1]
ax2 = plt.subplot(gs[:2, 2:])

ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('', loc = 'left')
ax2.set_ylabel('', loc = 'bottom')
ax2.set_title('MDA-MB-231 Subpopulations')

for cat in identity.unique():
    if cat in ['C2', 'ESAM Neg', 'ESAM Pos']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]


        ax2.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])

lgnd = ax2.legend(prop=dict(size=15), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

ax3 = plt.subplot(gs[2:4, 1:3])

ax3.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel('', loc = 'left')
ax3.set_ylabel('', loc = 'bottom')
ax3.set_title('MDA-MB-436 Subpopulations')

for cat in identity.unique():
    if cat in ['0']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]

        ax3.scatter(X, Y, c = colorDict[cat], s = 25, label = allLabelDict[cat])
    if cat in ['1']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]

        ax3.scatter(X, Y, c = colorDict[cat], s = 25, label = allLabelDict[cat])
lgnd = ax3.legend(prop=dict(size=15), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

lgnd.get_frame().set_linewidth(0.0)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
# fig.savefig('../../figures/umapConcat.png', dpi = 500, bbox_inches='tight')
# %%
umappts = adata.obsm['X_umap'].copy()
identity = adata.obs['sample'].copy()
plt.rcParams.update({'font.size': 18})

# plt.rcParams["axes.spines.right"] = False
# plt.rcParams["axes.spines.top"] = False

# fig, axs = plt.subplots(1, 3, figsize = (20,16))
fig = plt.figure(figsize = (30, 10))
gs = gridspec.GridSpec(1, 3)
gs.update(wspace = 0.5, hspace = 0.5)
ax1 = plt.subplot(gs[0])
# ax1 = axs[0]

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('UMAP 1', loc = 'left')
ax1.set_ylabel('UMAP 2', loc = 'bottom')
ax1.set_title('MDA-MB-231 Samples')

for cat in identity.unique():
    if cat in ['C2', 'PT']:
        isSample = identity == cat
        X = umappts[isSample, 0].copy()
        Y = umappts[isSample, 1].copy()

        ax1.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])
lgnd = ax1.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()


ax1.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 2, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)

ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)

identity = adata.obs['subsample']
# ax2 = axs[1]
ax2 = plt.subplot(gs[1])

ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('', loc = 'left')
ax2.set_ylabel('', loc = 'bottom')
ax2.set_title('MDA-MB-231 Subpopulations')

for cat in identity.unique():
    if cat in ['C2', 'ESAM Neg', 'ESAM Pos']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]


        ax2.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])

lgnd = ax2.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

ax3 = plt.subplot(gs[2])

ax3.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_xlabel('', loc = 'left')
ax3.set_ylabel('', loc = 'bottom')
ax3.set_title('BT474 Subpopulations')

for cat in identity.unique():
    if cat in ['0']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]

        ax3.scatter(X, Y, c = colorDict[cat], s = 25, label = allLabelDict[cat])
    if cat in ['1']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]

        ax3.scatter(X, Y, c = colorDict[cat], s = 25, label = allLabelDict[cat])
lgnd = ax3.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

lgnd.get_frame().set_linewidth(0.0)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
fig.savefig('../../figures/umapConcatHoriz.png', dpi = 500, bbox_inches='tight')
# %%
umappts = adata.obsm['X_umap'].copy()
identity = adata.obs['sample'].copy()
plt.rcParams.update({'font.size': 18})

# plt.rcParams["axes.spines.right"] = False
# plt.rcParams["axes.spines.top"] = False

# fig, axs = plt.subplots(1, 3, figsize = (20,16))
fig = plt.figure(figsize = (25/3, 10))
gs = gridspec.GridSpec(1, 1)
gs.update(wspace = 0.5, hspace = 0.5)
ax2 = plt.subplot(gs[0])


identity = adata.obs['subsample']
# ax2 = axs[1]

ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('', loc = 'left')
ax2.set_ylabel('', loc = 'bottom')
ax2.set_title('MDA-MB-231 Subpopulations')

esamVal = adataEsam[:, 'ESAM'].X.toarray().ravel()
for cat in identity.unique():
    if cat in ['ESAM Neg', 'ESAM Pos']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]


        ax2.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])

for cat in identity.unique():
    if cat in ['C2']:
        isSample = identity == cat
        X = umappts[isSample, 0]
        Y = umappts[isSample, 1]


        ax2.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat])

lgnd = ax2.legend(prop=dict(size=18), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

# %%
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
ax1.plot(range(0,10), range(0,10))
ax1.set_title('1')
ax2 = plt.subplot(gs[:2, 2:])
ax2.plot(range(0,10), range(0,10))

ax3 = plt.subplot(gs[2:4, 1:3])
ax3.plot(range(0,10), range(0,10))

plt.show()
# %%
colorDict['0'] = colorDict['ESAM Neg']
colorDict['1'] = colorDict['ESAM Pos']
allLabelDict['0'] = 'Subpop 1'
allLabelDict['1'] = 'Subpop 2'

umappts = adataEsam.obsm['X_umap'].copy()
identity = adataEsam.obs['leiden'].copy()
plt.rcParams.update({'font.size': 18})

# plt.rcParams["axes.spines.right"] = False
# plt.rcParams["axes.spines.top"] = False

# fig, axs = plt.subplots(1, 3, figsize = (20,16))
fig = plt.figure(figsize = (8, 8))
gs = gridspec.GridSpec(1,1)
gs.update(wspace = 0.5, hspace = 0.5)
ax1 = plt.subplot(gs[0])
# ax1 = axs[0]

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('UMAP 1', loc = 'left')
ax1.set_ylabel('UMAP 2', loc = 'bottom')
ax1.set_title('MDA-MB-231 Subpopulations')

for cat in identity.unique():
    
    isSample = identity == cat
    X = umappts[isSample, 0].copy()
    Y = umappts[isSample, 1].copy()

    ax1.scatter(X, Y, c = colorDict[cat], s = 10, label = allLabelDict[cat]) 
lgnd = ax1.legend(prop=dict(size=15), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([100])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()


ax1.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 1.3, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)

ax1.xaxis.label.set_fontsize(15)
ax1.yaxis.label.set_fontsize(15)

fig.savefig('../../figures/mdamb231SubpopNormalized.png', dpi = 500,  bbox_inches='tight')