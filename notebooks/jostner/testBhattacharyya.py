# %%
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
# %%
def calculateKDE(y, covarianceFactor = 0.25):
    kernel = gaussian_kde(y)
    kernel.covariance_factor = lambda : covarianceFactor
    kernel._compute_covariance()
    return kernel

def bhattacharyya(p, q):
    return np.sum(np.sqrt(p*q))

y1 = np.random.normal(0, 1, 10000)
y2 = np.random.normal(1, 1, 10000)

y1 = np.random.exponential(1, 10000)
y2 = np.random.exponential(5, 10000)
# y2 = y1.copy()

x = np.linspace(-20, 20, 1000)
kernel1 = calculateKDE(y1)
kernel2 = calculateKDE(y2)

y1KDE = kernel1(x)
y2KDE = kernel2(x)
plt.plot(x, y1KDE/sum(y1KDE))
plt.plot(x, y2KDE/sum(y2KDE))


# %%
print(bhattacharyya(y1KDE/sum(y1KDE), y2KDE/sum(y2KDE)))
# %%
np.arange(0, len(x), 0.17)

kernel1(y1KDE)