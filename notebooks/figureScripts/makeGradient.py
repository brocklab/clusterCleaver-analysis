# %%
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# %%
img = plt.imread('../../figures/gradient.png')
plt.imshow(img)
# %%
rows, cols, _ = img.shape

rVals, gVals, bVals = [], [], []
for col in range(cols-10):
    r, g, b, a = list(np.mean(img[:,col], axis = 0))
    rVals.append(r)
    gVals.append(g)
    bVals.append(b)

# %%
x = list(range(0, len(rVals)))
plt.plot(x, rVals, c = 'r')
plt.plot(x, bVals, c = 'b')
plt.plot(x, gVals, c = 'g')
# %%
z = np.polyfit(x, rVals, 4)
pR = np.poly1d(z)

z = np.polyfit(x, gVals, 4)
pG = np.poly1d(z)

z = np.polyfit(x, bVals, 4)
pB = np.poly1d(z)

plt.plot(x, pR(x), c = 'r')
plt.plot(x, pB(x), c = 'b')
plt.plot(x, pG(x), c = 'g')
# %% Make colormap

# %%
def plot_examples(cms):
    """
    helper function to plot two colormaps
    """
    np.random.seed(19680801)
    data = np.random.randn(30, 30)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    for [ax, cmap] in zip(axs, cms):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()

N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(90/256, 1, N)
vals[:, 1] = np.linspace(39/256, 1, N)
vals[:, 2] = np.linspace(41/256, 1, N)
newcmp = ListedColormap(vals)
plot_examples([newcmp])
# %%
xI = np.round(np.linspace(0, 601, N))
vals = np.ones((256, 4))
vals[:, 0] = pR(xI)
vals[:, 1] = pG(xI)
vals[:, 2] = pB(xI)
newcmp = ListedColormap(vals)
plot_examples([newcmp])
# %%
import pandas as pd
pd.DataFrame(vals).to_csv('../../figures/customCmap.csv')