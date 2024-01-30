# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import wasserstein_distance, gaussian_kde
# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
fullPalette = list(colors + sns.color_palette("tab10"))
sns.set_palette(sns.color_palette(fullPalette))
# %%
def makeKDE(vals):
    density = gaussian_kde(vals)
    xPoints = np.linspace(min(vals)-1, max(vals) + 1, 5000)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    return xPoints, density(xPoints)

x1 = np.random.normal(0, 1, 10000)
x, y = makeKDE(x1)

plt.plot(x,y)
plt.fill_between(x, y, alpha = 0.5)
# %%
x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(3, 1, 10000)
EMD = wasserstein_distance(x1, x2)

x1, y1 = makeKDE(x1)

x2, y2 = makeKDE(x1)

plt.plot(x1,y1)
plt.fill_between(x, y, alpha = 0.5)

# %%
plt.figure(figsize = (12,4))

plt.subplot(131)
x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(0.75, 1, 10000)
EMD = wasserstein_distance(x1, x2)


plt.hist(x1, alpha = 0.5)
plt.hist(x2, alpha = 0.5)
plt.title(f'EMD: {EMD:0.2f}')

# plt.xticks([])
plt.yticks([])
plt.xlabel('Gene Expression')

plt.subplot(132)
x1 = np.random.normal(0, 1.5, 10000)
x2 = np.random.normal(2, 1, 10000)
EMD = wasserstein_distance(x1, x2)
plt.hist(x1, alpha = 0.5)
plt.hist(x2, alpha = 0.5)
plt.title(f'EMD: {EMD:0.2f}')

# plt.xticks([])
plt.yticks([])
plt.xlabel('Gene Expression')

plt.subplot(133)
x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(8, 1, 10000)
EMD = wasserstein_distance(x1, x2)
plt.hist(x1, alpha = 0.5)
plt.hist(x2, alpha = 0.5)
plt.title(f'EMD: {EMD:0.2f}')

# plt.xticks([])
plt.yticks([])
plt.xlabel('Gene Expression')


plt.savefig('../../figures/emdExemplar/emdExample.png', dpi = 500)
# %%
