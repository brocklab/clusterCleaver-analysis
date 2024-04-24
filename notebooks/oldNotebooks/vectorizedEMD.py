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
# %%
def vectorizedSearchSorted(a, b):
    """
    A vectorized implementation of the numpy function `searchsorted`.
    This is a columnwise implementation
    """
    m,n = a.shape
    max_num = np.maximum(a.max() - a.min(), b.max() - b.min()) + 1
    r = max_num*np.arange(a.shape[1])[None,:]
    p = np.searchsorted((a+r).ravel(order = 'F'), (b+r).ravel(order = 'F'), 'right').reshape(-1, n, order = 'F')
    res = p - m*(np.arange(n))
    return res
    # z = res
    # print(z)
def vectorizedWasserstein(u_values, v_values):
    """
    Computes the wasserstein distance for two values. This is based heavily
    on the scipy code, however it does not have weights. 

    Note: Wasserstein distances are currently only evaluated over columns, there is no axis value

    Inputs:
    u_values, v_values: Matrices to be computed

    Outputs:
    distance: Computed columnwise distance
    """
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

    u_cdf_indices = vectorizedSearchSorted(u_values, all_values[:-1])
    v_cdf_indices = vectorizedSearchSorted(v_values, all_values[:-1])

    # Calculate the CDFs of u and v using their weights, if specified.
    u_cdf = u_cdf_indices / u_values.shape[0]

    v_cdf = v_cdf_indices / v_values.shape[0]


    wd = np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas), axis = 0)

    return wd

def wasserstein2d(u_values, v_values):
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
    # u_cdf_indices = u_values.searchsorted(all_values[:-1], 'right')
    # v_cdf_indices = v_values.searchsorted(all_values[:-1], 'right')

    u_cdf_indices = np.array([np.searchsorted(u_values[:, i], all_values[:-1, i], 'right') for i in range(all_values.shape[1])]).T
    v_cdf_indices = np.array([np.searchsorted(v_values[:, i], all_values[:-1, i], 'right') for i in range(all_values.shape[1])]).T
    
    # Calculate the CDFs of u and v using their weights, if specified.
    u_cdf = u_cdf_indices / u_values.shape[0]

    v_cdf = v_cdf_indices / v_values.shape[0]


    wd = np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas), axis = 0)

    return wd
# %%
x = np.random.randn(1000)
y = np.random.randn(1000) + 1

print(wasserstein_distance(x, y))
print(vectorizedWasserstein(x[:,None], y[:,None]))
print(wasserstein2d(x[:,None], y[:,None]))
# %%
x = np.random.randn(1000, 2)
y = np.random.randn(2000, 2)

# print(wasserstein_distance(x, y))
print(vectorizedWasserstein(x, y))
print(wasserstein2d(x, y))
# %%
