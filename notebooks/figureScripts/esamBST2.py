# %%
import scanpy as sc
import os
import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns 
from scipy.stats import wasserstein_distance, gaussian_kde
from matplotlib.colors import ListedColormap
from scrna.cluster.main import compute_dimensionality_reductions
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
sns.set_palette(sns.color_palette(colors))
cmap = ListedColormap(pd.read_csv('../../figures/customCmap.csv', index_col=0).values)
# %%
adataFull = sc.read_h5ad('../../data/h5ads/jostner-processed-regressed-clustered.h5ad')
# %%
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
def makeKDE(x, vals):
    density = gaussian_kde(vals)
    xPoints = np.linspace(min(x), max(x), 5000)
    counts = density(xPoints)
    return xPoints, counts 

# %%
adata231 = adatas['mdamb231']
vals = adata231[:,'ESAM'].layers['log1p_norm'].toarray().ravel()
is0 = adata231.obs['leiden'] == '0'
is1 = adata231.obs['leiden'] == '1'

# %%
flowImg = plt.imread('../../figures/esamFlow.png')
dfHist = pd.DataFrame(vals, adata231.obs['leiden']).reset_index()
dfHist.columns = ['leiden', 'expression']
x = [np.min(vals), np.max(vals)]

KDE0x, KDE0y = makeKDE(x, vals[is0])
KDE1x, KDE1y = makeKDE(x, vals[is1])

KDE0y = KDE0y
KDE1y = KDE1y
# %%
matplotlib.rcParams.update({'font.size': 15})
umappts = adata231.obsm['X_umap'].copy()
inner = [['histogram'],
         ['stripplot']]
outer = [['umap', inner]]

fig, axd = plt.subplot_mosaic(outer, 
                              layout="constrained", 
                              figsize = (16,6),
                              width_ratios=[1, 1.1])
ax1 = axd['umap']
map = ax1.scatter(umappts[:,0], umappts[:,1], c = vals, cmap = 'viridis')
ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])
xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()
ax1.set_xlabel('UMAP 1', loc = 'left', fontsize = 13)
ax1.set_ylabel('UMAP 2', loc = 'bottom', fontsize = 13)

ax1.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 2, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)
ax2 = axd['histogram']
ax2.plot(KDE0x, KDE0y)
ax2.fill_between(KDE0x, KDE0y, alpha = 0.5, label = 0)
ax2.plot(KDE1x, KDE1y)
ax2.fill_between(KDE1x, KDE1y, alpha = 0.5, label = 1)
ax2.legend(title = 'Leiden Cluster')
# ax2.set_xlabel('Normalized ESAM Expression')
ax2.set_ylabel('Normalized Density')
ax2.set_yticks([])
ax2.set_xticks([])
ax2.spines[['top', 'right']].set_visible(False)

ax3 = axd['stripplot']
sns.stripplot(
        data=dfHist, 
        x='expression', 
        hue='leiden', 
        native_scale=True,
        legend = False,
        jitter = 0.45,
        # jitter = True,
        ax = ax3,
        alpha = 0.75).set(
    xlabel = f'Normalized ESAM Expression'
)
ax3.spines[['top', 'left', 'right']].set_visible(False)

fig.tight_layout()

fig.colorbar(map, ax = ax1)
fig.savefig('../../figures/final/231Expression.png', dpi = 500)
# %%
adata231 = adatas['mdamb436']
vals = adata231[:,'BST2'].layers['log1p_norm'].toarray().ravel()
is0 = adata231.obs['leiden'] == '0'
is1 = adata231.obs['leiden'] == '1'
# %%
umappts = adata231.obsm['X_umap'].copy()
inner = [['histogram'],
         ['stripplot']]
outer = [['umap', inner]]

fig, axd = plt.subplot_mosaic(outer, 
                              layout="constrained", 
                              figsize = (16,6),
                              width_ratios=[1, 1.1])
ax1 = axd['umap']
map = ax1.scatter(umappts[:,0], umappts[:,1], c = vals, cmap = 'viridis')
ax1.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])
xmin, xmax = ax1.get_xlim() 
ymin, ymax = ax1.get_ylim()
ax1.set_xlabel('UMAP 1', loc = 'left', fontsize = 13)
ax1.set_ylabel('UMAP 2', loc = 'bottom', fontsize = 13)

ax1.arrow(xmin, ymin, 1.8, 0., fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 


ax1.arrow(xmin, ymin, 0., 1.6, fc='k', ec='k', lw = 1, 
         head_width=0.25, head_length=0.25, overhang = 0.3, 
         length_includes_head= False, clip_on = False) 

ax1.xaxis.set_label_coords(0.05, 0.025)
ax1.yaxis.set_label_coords(0.025, 0.05)
ax2 = axd['histogram']
ax2.plot(KDE0x, KDE0y)
ax2.fill_between(KDE0x, KDE0y, alpha = 0.5, label = 0)
ax2.plot(KDE1x, KDE1y)
ax2.fill_between(KDE1x, KDE1y, alpha = 0.5, label = 1)
ax2.legend(title = 'Leiden Cluster')
# ax2.set_xlabel('Normalized ESAM Expression')
ax2.set_ylabel('Normalized Density')
ax2.set_yticks([])
ax2.set_xticks([])
ax2.spines[['top', 'right']].set_visible(False)

ax3 = axd['stripplot']
sns.stripplot(
        data=dfHist, 
        x='expression', 
        hue='leiden', 
        native_scale=True,
        legend = False,
        jitter = 0.45,
        # jitter = True,
        ax = ax3,
        alpha = 0.75).set(
    xlabel = f'Normalized BST2 Expression'
)
ax3.spines[['top', 'left', 'right']].set_visible(False)

fig.tight_layout()

fig.colorbar(map, ax = ax1)
fig.savefig('../../figures/final/436Expression.png', dpi = 500)