# %%
import numpy as np
from scipy.stats import wasserstein_distance
# %%
x = np.random.randn(1000, 2)
y = x.copy()

# x = np.random.randn(1000)
# y = np.random.randn(1000) + 1

# x = x.reshape(-1, 1)
# y = y.reshape(-1, 1)
# %%
# print(wasserstein_distance(x, y))
# %%
u_values = x
v_values = y

u_weights = None
v_weights = None
# %%


u_sorter = np.argsort(u_values, axis = 0)
v_sorter = np.argsort(v_values, axis = 0)

all_values = np.concatenate((u_values, v_values), axis = 0)
all_values.sort(kind='mergesort', axis = 0)

# Compute the differences between pairs of successive values of u and v.
deltas = np.diff(all_values, axis = 0)
# %%
# Get the respective positions of the values of u and v among the values of
# both distributions.

u_values.sort(axis = 0)
v_values.sort(axis = 0)
# u_cdf_indices = u_values.searchsorted(all_values[:-1], 'right')
# v_cdf_indices = v_values.searchsorted(all_values[:-1], 'right')

u_cdf_indices = np.array([np.searchsorted(u_values[:, i], all_values[:-1, i], 'right') for i in range(all_values.shape[1])]).T
v_cdf_indices = np.array([np.searchsorted(v_values[:, i], all_values[:-1, i], 'right') for i in range(all_values.shape[1])]).T

u_cdf_indices_2 = 0
v_cdf_indices_2 = 0


# for i in range(0, all_values.shape[1]):
#     u_cdf_indices = np.searchsorted(u_values[:, i][u_sorter[:, i]],all_values[:-1, i],side="r") 

# u_cdf_indices = np.array([np.searchsorted(u_values[:, i][u_sorter[:, i]],all_values[:-1, i],side="r") for i in range(all_values.shape[1])])
# v_cdf_indices = np.array([np.searchsorted(v_values[:, i][v_sorter[:, i]],all_values[:-1, i],side="r") for i in range(all_values.shape[1])])

# %%
# Calculate the CDFs of u and v using their weights, if specified.
if u_weights is None:
    u_cdf = u_cdf_indices / u_values.shape[0]

if v_weights is None:
    v_cdf = v_cdf_indices / v_values.shape[0]

 #%%

wd = np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas), axis = 0)
print(wd)
# %%
x = np.random.randn(1000)
y = np.random.randn(1000) + 1

# print(wasserstein_distance(x, y))
u_values = x
v_values = y

u_weights = None
v_weights = None


u_sorter = np.argsort(u_values, axis = 0)
v_sorter = np.argsort(v_values, axis = 0)

all_values = np.concatenate((u_values, v_values), axis = 0)
all_values.sort(kind='mergesort', axis = 0)

# Compute the differences between pairs of successive values of u and v.
deltas = np.diff(all_values, axis = 0)
# Get the respective positions of the values of u and v among the values of
# both distributions.

u_values.sort(axis = 0)
v_values.sort(axis = 0)
u_cdf_indices = u_values.searchsorted(all_values[:-1], 'right')
v_cdf_indices = v_values.searchsorted(all_values[:-1], 'right')

u_cdf_indices / u_values.size

# %% Fix searchsorted

a = u_values.T.copy()
b = all_values[:-1, :].T.copy()

m,n = a.shape
max_num = np.maximum(a.max() - a.min(), b.max() - b.min()) + 1
r = max_num*np.arange(a.shape[0])[:,None]
p = np.searchsorted( (a+r).ravel(), (b+r).ravel() ).reshape(m,-1)
res = p - n*(np.arange(m)[:,None]).T

assert(np.all(u_cdf_indices == res.T+1))
# %%