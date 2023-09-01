# %%
import os
# %%
files = os.listdir('./data/h5ads/')
files = [fileName for fileName in files if fileName.startswith("dave") and fileName.endswith('-raw.h5ad') ]
# %%
for fileName in files:
    fileName = fileName.split('-raw.h5ad')[0]
    basePath = f'./data/h5ads/{fileName}'
    postqcPath = f'{basePath}-postqc'
    normalizedPath = f'{postqcPath}-normalized'
    clusteredPath = f'{normalizedPath}-clustered'

    os.system(f'snakemake {postqcPath}.h5ad --cores all')
    os.system(f'snakemake {normalizedPath}.h5ad --cores all')
    os.system(f'snakemake {clusteredPath}.h5ad --cores all')