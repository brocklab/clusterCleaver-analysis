# %%
import numpy as np
import diptest
import matplotlib.pyplot as plt
# %%
# generate some bimodal random draws
N = 1000000
hN = N // 2
x = np.empty(N, dtype=np.float64)
x[:hN] = np.random.normal(1.5, 1.0, hN)
x[hN:] = np.random.normal(-1, 1.0, hN)

# only the dip statistic
dip = diptest.dipstat(x)

# both the dip statistic and p-value
dip, pval = diptest.diptest(x)

plt.hist(x)
plt.title(f'p = {pval:0.2f}')
# %%
