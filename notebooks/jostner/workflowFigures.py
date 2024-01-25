# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# %%
cov = [[6, -3], [-3, 3.5]]
c1 = np.random.multivariate_normal((1,2), cov, 1000)
cov2 = [[-1, 3], [-10, 1]]
c2 = np.random.multivariate_normal((3,8), cov2, 1000)

plt.scatter(c1[:,0], c1[:,1], c = 'red')
plt.scatter(c2[:,0], c2[:,1], c = 'blue')
# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
fullPalette = list(colors + sns.color_palette("tab10"))
sns.set_palette(sns.color_palette(fullPalette))
# %%
plt.scatter(0, 1)
plt.scatter(1, 1)
plt.scatter(2, 2)
plt.scatter(3, 3)
plt.scatter(4, 4)
