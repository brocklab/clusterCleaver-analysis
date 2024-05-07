# Computational identification of surface markers for isolating distinct subpopulations from heterogeneous cancer cell populations

This is a set of notebooks/scripts used to process single-cell data and find surface markers for the paper (biorxiv reference). We also released a package for this called `scDistRank` available on pypi at https://pypi.org/project/scDistRank/. 
## Environment Setup

To get started generate and activate the `conda` environment.

> **Note**
> must have one of `conda`, `mamba` or `micromamba`
```sh
make env
```

Then activate the environment with the `conda` package manager you use:

```sh
micromamba activate ./env
```

Environments should only be made using the `locked.yml`
If you need to add an additional package which namespace is directly accessed.
Then it should be added to the `env.yml` and the `locked.yml` re-resolved with the below command. Which will generate a new environment and update the `locked.yml`

```sh
make locked.yml
```
