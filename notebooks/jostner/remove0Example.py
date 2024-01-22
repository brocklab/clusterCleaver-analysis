# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import wasserstein_distance
# %%
x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(3, 1, 10000)
EMD = wasserstein_distance(x1, x2)

plt.figure(figsize = (8,4))
plt.subplot(121)
plt.hist(x1, alpha = 0.5)
plt.hist(x2, alpha = 0.5)
plt.title(f'EMD: {EMD:0.2f}')

x1 = np.random.normal(0, 1, 10000)
x2 = np.random.normal(0.75, 1, 10000)
EMD = wasserstein_distance(x1, x2)

plt.subplot(122)
plt.hist(x1, alpha = 0.5)
plt.hist(x2, alpha = 0.5)
plt.title(f'EMD: {EMD:0.2f}')

plt.savefig('../../figures/emdExample.png', dpi = 500)
# %%
