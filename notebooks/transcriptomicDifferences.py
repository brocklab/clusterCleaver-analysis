# %%
from scrna.cluster.main import compute_dimensionality_reductions

import scanpy as sc
import plotnine
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# %%
adata231 = sc.read_h5ad('../data/h5ads/231-1KB3-20231013.h5ad')
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
adatabt474 = sc.read_h5ad('../data/h5ads/BT474LineageAssigned.h5ad')
adatabt474 = adatabt474[adatabt474.obs['sample'] == 'LPD7']
sc.pp.highly_variable_genes(adatabt474)
compute_dimensionality_reductions(adatabt474)
# %%
lin1 = 'GTGAGTCAGTCACACAGTCA'
isLin1 = adatabt474.obs['collapsedLineages'] == lin1
adatabt474.obs['lin1'] = 'lpdOther'
adatabt474.obs.loc[isLin1, 'lin1'] = 'lin1'
adatabt474.obs['lin1'] = adatabt474.obs['lin1'].astype('category')
adatabt474.obs['subsample'] = 'LPDOther'
adatabt474.obs.loc[isLin1, 'subsample'] = 'Lin1'
# %%
adata = sc.concat([adata231, adatabt474])
# %%
sc.pp.highly_variable_genes(adata)
compute_dimensionality_reductions(adata)
# %%
sc.pl.umap(adata, color = ['sample', 'subsample'])
# %%
sc.tl.dendrogram(adata, groupby = 'subsample')
sc.pl.dendrogram(adata, groupby = 'subsample')
# %%
sampleDict = { 'C2':    'MDAMB231\nDox Treated',
               'LPD7':  'BT474\nLap/Pac/Dox\nTreated',
               'PT':    'MDAMB231\nUntreated'
}

subsampleDict = {'ESAM Neg': 'ESAM (-)',
                 'ESAM Pos': 'ESAM (+)',
                 'LPDOther': 'BT474\nLap/Pac/Dox\nOther',
                 'Lin1':     'BT474\nLap/Pac/Dox\nLineage 1',
                 'treated':  'MDAMB231\nDox Treated'}

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

X   = adata.X
C2      =   X[isC2, :].mean(axis = 0).tolist()[0]
PT      =   X[isPT, :].mean(axis = 0).tolist()[0]
ESAMPos =   X[isESAMPos, :].mean(axis = 0).tolist()[0]
ESAMNeg =   X[isESAMNeg, :].mean(axis = 0).tolist()[0]
lin1  =     X[isLin1, :].mean(axis = 0).tolist()[0]
lpdOther =  X[isOther, :].mean(axis = 0).tolist()[0]
# %%
pcaVals = {'C2': C2, 'PT': PT, 'ESAM Pos': ESAMPos, 'ESAM Neg': ESAMNeg, 'Lin1': lin1, 'LPDOther': lpdOther}
for key in list(pcaVals.keys()):
    newKey = ''
    if key in sampleDict.keys():
        newKey = sampleDict[key]
    elif key in subsampleDict.keys():
        newKey = subsampleDict[key]

    pcaVals[newKey] = pcaVals.pop(key)   
pcaVals = pd.DataFrame(pcaVals)
pcaCorr = pcaVals.corr()

sns.heatmap(pcaCorr, annot=True, fmt='.3').set_title('Pearson Correlation of\nPopulations')
sc.pl.umap(adata, color = ['sampleNames', 'subsampleNames'], wspace = .25)
# %%
colorDict = {'LPD7':        '#002bff',
            #  'treated':       '#F68511',
             'ESAM Pos':        '#8fc26c',
             'ESAM Neg':        '#e43820',
             'C2':              '#4bc1ee',
             'PT':              '#7d1811',
             'LPDOther':            '#2e336a',
             'Lin1':            '#fed504'
             }

# %% Plot samples

umappts = adata.obsm['X_umap']
identity = adata.obs['sample']
plt.rcParams.update({'font.size': 18})

# plt.rcParams["axes.spines.right"] = False
# plt.rcParams["axes.spines.top"] = False
fig, axs = plt.subplots(1, 2, figsize = (20,8))
ax1 = axs[0]

ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_xlabel('UMAP 1', loc = 'left')
ax1.set_ylabel('UMAP 2', loc = 'bottom')
ax1.set_title('Samples')

for cat in identity.unique():
    isSample = identity == cat
    X = umappts[isSample, 0]
    Y = umappts[isSample, 1]

    # if cat == 'LPDOther':
    #     alpha = 1
    # else:
    #     alpha = 1
    ax1.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat], alpha = alpha)
lgnd = ax1.legend(prop=dict(size=10), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([80])
lgnd.get_frame().set_linewidth(0.0)

xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()


ax1.arrow(xmin, ymin, 3, 0., fc='k', ec='k', lw = 1, 
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
ax2 = axs[1]
ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('', loc = 'left')
ax2.set_ylabel('', loc = 'bottom')
ax2.set_title('Sub-Samples')

for cat in identity.unique():
    isSample = identity == cat
    X = umappts[isSample, 0]
    Y = umappts[isSample, 1]

    # if cat == 'LPDOther':
    #     alpha = 1
    # else:
    #     alpha = 1
    ax2.scatter(X, Y, c = colorDict[cat], s = 2, label = allLabelDict[cat], alpha = alpha)

lgnd = ax2.legend(prop=dict(size=10), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([80])
lgnd.get_frame().set_linewidth(0.0)

# %%

# %%
adata231 = sc.read_h5ad('../data/h5ads/231-1KB3-20231013.h5ad')
# %%
# sc.pl.umap(adata231[adata231.obs['sample'].isin(['PT', 'C2'])], color = 'sample')
adata231 = adata231[adata231.obs['sample'].isin(['PT', 'C2'])]

# sc.pp.highly_variable_genes(adata231)
# compute_dimensionality_reductions(adata231)
X_pca = adata231.X
samplePCA = {}
for sample in adata231.obs['sample'].unique():
    samplePCA[sample] = X_pca[adata231.obs['sample'] == sample, :].mean(axis = 0).tolist()[0]
pcaVals = pd.DataFrame(samplePCA)
pcaCorr = pcaVals.corr(method = 'spearman')

sns.heatmap(pcaCorr, annot=True, fmt='.3').set_title('Pearson Correlation of\nPopulations')
sc.pl.umap(adata231, color = 'sample')
sc.pl.pca(adata231, color = 'sample')
# %%
adatabt474 = sc.read_h5ad('../data/h5ads/BT474LineageAssigned.h5ad')
adatabt474 = adatabt474[adatabt474.obs['sample'].isin(['LPD7', 'DLP5'])]
sc.pp.highly_variable_genes(adatabt474)
compute_dimensionality_reductions(adatabt474)
# %%
X_pca = adatabt474.obsm['X_pca']
samplePCA = {}
for sample in adatabt474.obs['sample'].unique():
    samplePCA[sample] = X_pca[adatabt474.obs['sample'] == sample, :].mean(axis = 0)
pcaVals = pd.DataFrame(samplePCA)
pcaCorr = pcaVals.corr()
pcaVals = pd.DataFrame(samplePCA)
pcaCorr = pcaVals.corr()

sns.heatmap(pcaCorr, annot=True, fmt='.3').set_title('Pearson Correlation of\nPopulations')
sc.pl.umap(adatabt474, color = 'sample')
sc.pl.pca(adatabt474, color = 'sample')
# %%import seaborn as sns
from sklearn.decomposition import PCA
iris = sns.load_dataset('iris')
iris = iris.loc[iris['species'].isin(['setosa', 'virginica'])]
pca = PCA(n_components = 2)
X = iris.loc[:, iris.columns != 'species'].to_numpy()
pcaVals = pca.fit_transform(X)

samplePCA = {}
for sample in iris['species'].unique():
    samplePCA[sample] = pcaVals[iris['species'] == sample, :].mean(axis = 0)
pcaVals = pd.DataFrame(samplePCA)
pcaCorr = pcaVals.corr()
# %%
