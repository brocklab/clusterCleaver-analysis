# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.sparse
# %%
colors = ['#BB4E44', '#44B1BB', '#76BB44', '#8944BB']
fullPalette = list(colors + sns.color_palette("tab10"))
sns.set_palette(sns.color_palette(fullPalette))
# %%
def plotExpression(adata, genes, colorCol = 'leiden'):
    X = adata[:, genes].X
    if scipy.sparse.issparse(X):
        X = X.toarray()
    
    dfExpr = pd.DataFrame([X[:,0], X[:, 1], adata.obs['leiden']]).T
    dfExpr.columns = [genes[0], genes[1], 'leiden']
    sns.jointplot(data = dfExpr, x = genes[0], y = genes[1], hue = 'leiden')
    
def plotHists(adata, gene, colorCol = 'leiden', logScale = False, saveFig = ''):
    surfaceIdx = np.where(adata.var.index.isin([gene]))[0][0]
    expression = adata.X[:, surfaceIdx]
    if scipy.sparse.issparse(expression):
        expression = expression.toarray()

    dfHist = pd.DataFrame(expression, adata.obs[colorCol]).reset_index()
    dfHist.columns = [colorCol, 'expression']

    dfHist[colorCol] = dfHist[colorCol].astype('category')
    # , log_scale=(False, True)
    plt.figure(figsize=(7,6))
    plt.subplot(211)
    sns.histplot(
                data=dfHist, 
                x='expression', 
                hue=colorCol, 
                element="poly", 
                stat='proportion',
                log_scale = (False, logScale)).set(xlabel='')
            
    plt.subplot(212)
    sns.stripplot(
            data=dfHist, 
            x='expression', 
            hue=colorCol, 
            native_scale=True,
            legend = False,
            # jitter = 0.45
            jitter = True).set(
        xlabel = f'{gene} Expression'
    )
    if len(saveFig) > 0:
        saveDirectory = Path(saveFig).parents[0]
        if saveDirectory.exists():
            plt.savefig(saveFig, dpi = 500)
        else:
            print('Save path directory {saveDirectory} does not exist.')

def plotParetoOptimal(optimalGenes, paretoOptimalGenes, nGenes = 1, metric = 'auc'):
    if nGenes not in [1, 2]:
        print('Number of genes must be one or two for plotting')
        return
    if nGenes == 1:
        plt.scatter(optimalGenes[metric], optimalGenes['finalScore1'])
        plt.scatter(paretoOptimalGenes[metric], paretoOptimalGenes['finalScore1'], c = 'red')
        plt.xlabel('AUC')
        plt.ylabel('Surface Score')
    elif nGenes == 2:
        ax = plt.figure(figsize=(8,8)).add_subplot(projection='3d')
        ax.scatter(optimalGenes['finalScore1'], optimalGenes['finalScore2'], optimalGenes[metric], c = 'blue', alpha = 0.05)
        ax.scatter(paretoOptimalGenes['finalScore1'], paretoOptimalGenes['finalScore2'], paretoOptimalGenes[metric], c = 'red')
        ax.dist = 12
        ax.set_xlabel('Surface Score 1')
        ax.set_ylabel('Surface Score 2')
        ax.set_zlabel('Combined AUC')
        ax.set_title('2-Gene Pareto Optimal')

