# %%
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
import os
import matplotlib.pyplot as plt
import matplotlib
# %%
adata = sc.read_h5ad('../../data/h5ads/HCC1806LineageAssigned.h5ad')
# %%
sc.pl.umap(adata)
# %%
lineages = adata.obs['collapsedLineages']
lineages = lineages[adata.obs['label'] == 'Invasive']
lineages = lineages[~lineages.isna()]

linCt = lineages.value_counts()
expData = np.cumsum(linCt)/sum(linCt)

plt.figure()
plt.plot(range(0,len(expData)),expData)
plt.vlines(x=100,ymin=0,ymax=1,colors='red')
plt.ylabel('Proportion of All Reads')
plt.xlabel('Number of Lineages')
plt.title('Cumulative Proportion of Lineage Abundance')
plt.grid()
plt.show()
# %%
plt.figure()
plt.plot(range(0,len(linCt)),linCt)
plt.vlines(x=10,ymin=0,ymax=np.max(linCt),colors='red')
plt.xlabel('Number of Lineages')
plt.ylabel('Number of Reads')
plt.title('Proportion of Lineage Abundance')
plt.grid()
# %%
topLins = list(expData[0:10].index)
# %%
alphabet = [x for x in 'abcdefghijklmnopqrstuvwxyz']
abbr2 = alphabet.copy()
for letter in alphabet:
    abbr2.append(letter*2)
topLinAbbr = {k:v for k,v in zip(topLins, abbr2[0:len(topLins)])}
# %% Perform unsupervised clustering on top lineages
adataTop = adata[adata.obs['collapsedLineages'].isin(topLins),:]
adataTop.obs['topLins'] = adataTop.obs['collapsedLineages'].astype('category')
ZLabels = list(adataTop.obs['topLins'].cat.categories)
Zabb = [topLinAbbr[lab] for lab in ZLabels]
# random.shuffle(ZLabels)

adataTop.obs['topLins'] = adataTop.obs['topLins'].cat.reorder_categories(ZLabels)
sc.tl.dendrogram(adataTop, groupby='topLins', linkage_method='complete', cor_method = 'pearson',key_added='lineageDendro')
thresh = 1.1
Z = adataTop.uns['lineageDendro']['linkage']
# %%
unsupColorsVals = ['#db5f57','#d3db57','#57db5f','#57d2db','#5457db','#db57d3']
unsupColors = {str(x):color for x,color in zip(range(1,6), unsupColorsVals)}
matplotlib.rcParams.update({'font.size': 16})
plt.figure(figsize=(6,4))
# ax = plt.gca()
# xlbls = ax.get_xmajorticklabels()
hierarchy.set_link_color_palette(list(unsupColorsVals))
hierarchy.dendrogram(Z, color_threshold=thresh, labels=Zabb, above_threshold_color='grey', leaf_rotation=0)
plt.axhline(y=thresh, linestyle='--',c='red')
plt.title('Unsupervised Lineage Clustering')
plt.xlabel('Lineage')
plt.ylabel('Expression Distance')
plt.xticks([])
plt.yticks([])
clustering = list(hierarchy.fcluster(Z, t=thresh, criterion='distance'))
labels = pd.DataFrame([ZLabels, clustering], index=None).transpose()
labels.columns = ['topLins','group']
labels['group'] = labels['group'].astype('string')
labels['color'] = labels['group'].replace(unsupColors)
labels['abbr'] = labels['topLins'].replace(topLinAbbr)
# %%
linDict = {k:v for k,v in zip(labels['topLins'], labels['group'])}
adata.obs['unsupGroup'] = adata.obs['collapsedLineages']
adata.obs['unsupGroup'] = adata.obs['unsupGroup'].astype('string')
adata.obs['unsupGroup'] = adata.obs['unsupGroup'].replace(linDict, regex=False)
adata.obs['unsupGroup'] = adata.obs['unsupGroup'].astype('string')
# %%
adataGroup = adata[adata.obs['unsupGroup'].isin(set(linDict.values()))]
sc.pl.umap(adataGroup, color = 'unsupGroup')