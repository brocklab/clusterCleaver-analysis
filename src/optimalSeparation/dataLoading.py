from pathlib import Path
import pandas as pd
from tqdm import tqdm
import scanpy as sc
import pickle

from scrna.cluster.main import compute_dimensionality_reductions


def cleanSurfaceGenes(homePath):
    """Loads csv of surface genes and their scores"""
    homePath = Path(homePath)
    surfaceGenes = pd.read_csv(homePath / 'data/optimalGenes/surfaceGenes.csv', skiprows=2)
    surfaceGenes = surfaceGenes.iloc[2:,]
    surfaceCols = list(surfaceGenes.columns)
    surfaceCols[0:5] = ['ensemble', 'gene', 'location', 'coreScore', 'finalScore']
    surfaceGenes.columns = surfaceCols
    surfaceGenes['finalScore'] = surfaceGenes['finalScore'].astype('float')

    return surfaceGenes

def loadClusteredAnndatas(adata, cellLines, adataPath = '../data/h5ads/kinder-leiden2ClustLines.pickle'):
    """
    Recomputes dimensionality reduction by cell line
    
    Inputs:
        - adata: Anndata object with multiple cell lines
        - cellLines: Cell lines to extract

    Outputs:
        - adatas: A dictionary of anndatas with cell line keys
    """
    
    adataPath = Path(adataPath)
    
    if adataPath.exists():
        print('Loading existing dictionary of anndatas')
        adatas = pickle.load(open(adataPath,"rb"))
    else:
        adatas = {}

        for cellLine in tqdm(cellLines):
            adata_sub = adata[adata.obs['Cell_line'] == cellLine]
            compute_dimensionality_reductions(adata_sub)
            sc.tl.leiden(adata_sub, resolution=0.1)
            adatas[cellLine] = adata_sub

    return adatas

def processFullAnndata(adataFull):
    """
    Loads, preprocesses, and clusters concatenated anndatas loaded from adataPath
    Preprocessing consists of removing doublets, calculating highly variable genes, and computing dimensionality reduction
    Clustering comes from selectively reducing the leiden resolution until two clusters remain
    Inputs:
        - adataPath: Path to .h5ad file
    Outpts:
        - adatas: Dictionary where keys are cell lines and values are corresponding anndatas
    
    """
    samples = adataFull.obs['sample'].unique()
    adatas = {}
    for sample in samples:
        adataSub = adataFull[adataFull.obs['sample'].isin([sample])]
        adataSub = adataSub[adataSub.obs['scDblFinder_class'] == 'singlet']
        sc.pp.highly_variable_genes(adataSub, min_mean=0.0125, max_mean=3, min_disp=0.5)
        compute_dimensionality_reductions(adataSub)
        # sc.pl.umap(adataSub, title = sample)
        adatas[sample] = adataSub

    cellLineRes = {}
    for cellLine, adataSub in adatas.items():
        leidenResolution = 2
        nLeiden = 5
        c = 1
        while nLeiden != 2:
            leidenResolution /= 1.3
            sc.tl.leiden(adataSub, resolution= leidenResolution)
            nLeiden = len(adataSub.obs['leiden'].unique())
            c += 1
            if c > 20:
                leidenResolution = 0
                break
        cellLineRes[cellLine] = leidenResolution
        print(adataSub.obs['leiden'].unique())
        print(f'Cell Line: {cellLine} \t Resolution: {leidenResolution} \t nLeiden: {nLeiden}')