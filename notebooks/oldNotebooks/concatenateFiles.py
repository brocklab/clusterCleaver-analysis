
# %%
import os
from collections import defaultdict
# %%
files = os.listdir('../../data/rawData/JA24106/data')
# print(files)

samples = {}

for file in files:
    sampleIdentifier = '-'.join(file.split('-')[0:2])
    if not file.endswith('.gz'):
        continue
    print(sampleIdentifier)
    if sampleIdentifier not in samples.keys():
        samples[sampleIdentifier] = [file]
    else:
        samples[sampleIdentifier].append(file)
# %%
commands = []
for sample in samples.keys():
    assert len(samples[sample]) == 2, f'{sample}'
    command = f'echo {sample}; cat {sample}* > {sample}.fastq.gz'
    print(command)
    commands.append(command)

print(';'.join(commands))
# %%
# %% Make Config
for sample in samples.keys():
    print(f'- {sample}')
# %% Make snakemake commands
commands = []
for sample in samples.keys():
    echo = f'echo {sample};'
    command = f'snakemake --until data/trim/{sample}.trimmed.fastq --cores 4'
    commands.append(echo+command)
print(';'.join(commands))
# %% Make samplesheet
for sample in samples.keys():
    samplefq = f'{sample}.trimmed.fastq.gz'
    trimmedPath = f'../../data/rawData/JA24106/data/trim/{samplefq}'
    assert os.path.exists(trimmedPath), f'{trimmedPath} does not exist'
    print(f'{sample},{samplefq},auto')