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
    