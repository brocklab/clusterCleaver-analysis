# %%Ztabt474 = sc.read_h5ad('../data/h5ads/BT474LineageAssigned.h5ad')
adatabt474 = adatabt474[adatabt474.obs['sample'].isin(['PreTreat1', 'D5'])]
sc.pp.highly_variable_genes(adatabt474)
compute_dimensionality_reductions(adatabt474)
# %%
allLabelDict = {'PT': 'Pretreatment\n231',
                'PreTreat1': 'Pretreatment\nBT474',
                'D5': 'Doxorubicin\n60 nM',
                'C2': 'Doxorubicin\n550 nM'}
umappts =   adata231.obsm['X_umap']
identity =  adata231.obs['sample']
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
ax1.set_title('MDA-MB-231 Treatment')

for cat in identity.unique():
    isSample = identity == cat
    X = umappts[isSample, 0]
    Y = umappts[isSample, 1]

    # if cat == 'LPDOther':
    #     alpha = 1
    # else:
    #     alpha = 1
    ax1.scatter(X, Y, s = 2, label = allLabelDict[cat])
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

ax2 = axs[1]

umappts =   adatabt474.obsm['X_umap']
identity =  adatabt474.obs['sample']

ax2.spines[["top", "right", 'left', 'bottom']].set_visible(False)

ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xlabel('', loc = 'left')
ax2.set_ylabel('', loc = 'bottom')
ax2.set_title('BT474 Treatment')

for cat in ['PreTreat1', 'D5']:
    isSample = identity == cat
    X = umappts[isSample, 0]
    Y = umappts[isSample, 1]

    # if cat == 'LPDOther':
    #     alpha = 1
    # else:
    #     alpha = 1
    ax2.scatter(X, Y, s = 2, label = allLabelDict[cat])
lgnd = ax2.legend(prop=dict(size=10), bbox_to_anchor=(1, 0.5),
                         loc='center left', borderaxespad=0.)

for handle in lgnd.legend_handles:
    handle.set_sizes([80])
lgnd.get_frame().set_linewidth(0.0)

fig.savefig('../figures/treatmentComparison.png', dpi=500)
# %%
