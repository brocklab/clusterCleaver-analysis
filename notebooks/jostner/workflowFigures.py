# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

from scipy.stats import wasserstein_distance, gaussian_kde

from scrna.cluster.main import compute_dimensionality_reductions

# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
fullPalette = list(colors + sns.color_palette("tab10"))
sns.set_palette(sns.color_palette(fullPalette))
# %%

fig, ax = plt.subplots(1,1, figsize = (5,5))
cov = [[6, -3], [-3, 3.5]]
c1 = np.random.multivariate_normal((1,2), cov, 1000)
cov2 = [[-1, 3], [-10, 1]]
c2 = np.random.multivariate_normal((3,8), cov2, 1000)

ax.scatter(c1[:,0], c1[:,1])
ax.scatter(c2[:,0], c2[:,1])

ax.set_xticks([])
ax.set_yticks([])
ax.spines[["top", "right", 'left', 'bottom']].set_visible(False)
ax.set_xlabel('UMAP 1', loc = 'left')
ax.set_ylabel('UMAP 2', loc = 'bottom')
xmin, xmax = ax.get_xlim() 
ymin, ymax = ax.get_ylim()
ax.arrow(xmin, ymin, 4, 0., fc='k', ec='k', lw = 1, 
    head_width=0.25, head_length=0.25, overhang = 0.3, 
    length_includes_head= False, clip_on = False) 


ax.arrow(xmin, ymin, 0., 3.5, fc='k', ec='k', lw = 1, 
        head_width=0.25, head_length=0.25, overhang = 0.3, 
        length_includes_head= False, clip_on = False)

ax.xaxis.set_label_coords(0.05, 0.025)
ax.yaxis.set_label_coords(0.025, 0.05)

ax.xaxis.label.set_fontsize(10)
ax.yaxis.label.set_fontsize(10)

fig.savefig('../../figures/exampleUMAP.png', dpi = 500)
# %%
def makeKDE(vals):
    density = gaussian_kde(vals)
    xPoints = np.linspace(min(vals)-1, max(vals) + 1, 5000)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    return xPoints, density(xPoints)
# %%
x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(3, 1, 10000)
EMD = wasserstein_distance(x1, x2)

x1, y1 = makeKDE(x1)

x2, y2 = makeKDE(x2)
plt.figure(figsize=(5, 5))
plt.plot(x1,y1)
plt.fill_between(x1, y1, alpha = 0.5)

plt.plot(x2,y2)
plt.fill_between(x2, y2, alpha = 0.5)
plt.xticks([])
plt.yticks([])
plt.xlabel('Gene Expression')

plt.savefig('../../figures/emdWorkflowHistogram.png', dpi = 500)
# %%
# %% [markdown]
"""
This notebook will test a new experimental method using EMD to look for genes that may have overlap between subpopulations
but are still characterized by a high expressing population.
"""
# %%
import scanpy as sc
from optimalSeparation import searchOptimal, dataLoading, visualization
# %%
adata = sc.read_h5ad('../../data/h5ads/jostner-processed.h5ad')
adata = adata[adata.obs['sample'] == 'bt474']
adata = adata[adata.obs['scDblFinder_class'] == 'singlet']
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
compute_dimensionality_reductions(adata)

# %%
sc.tl.leiden(adata, resolution = 0.05)
sc.pl.umap(adata, color = 'leiden')
# %%
surfaceGenes = dataLoading.cleanSurfaceGenes('../..')

emdGenes = searchOptimal.searchExpressionDist(adata, surfaceGenes['gene'], modifier = None)
emdGenes.head()
# %%
gene = 'SYNE2'
is0 = adata.obs['leiden'] == '0'
is1 = adata.obs['leiden'] == '1'
is0 = is0.to_numpy().astype(bool)
is1 = is1.to_numpy().astype(bool)
x = adata[is0, gene].X.ravel()
y = adata[is1, gene].X.ravel()

x1, y1 = makeKDE(x)
x2, y2 = makeKDE(y)

fig, ax = plt.subplots(1, 1, figsize = (5,5))
ax.spines[["top", "right", 'left', 'bottom']].set_visible(False)
plt.plot(x1,y1)
plt.fill_between(x1, y1, alpha = 0.5)

plt.plot(x2,y2)
plt.fill_between(x2, y2, alpha = 0.5)
plt.xticks([])
plt.yticks([])
plt.xlabel('Gene Expression')

fig.savefig('../../figures/exampleHistogram.png', dpi = 500)
# %%