# %%
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.sparse
# %%
def plotExpression(adata, genes, colorCol = 'leiden'):
    X = adata[:, genes].X
    plt.scatter(X[:,0], X[:, 1], c = adata.obs[colorCol].astype(int))
    plt.xlabel(genes[0])
    plt.ylabel(genes[1])

def plotHists(adata, gene, colorCol = 'leiden'):
    surfaceIdx = np.where(adata.var.index.isin([gene]))[0][0]
    expression = adata.X[:, surfaceIdx]
    if scipy.sparse.issparse(expression):
        expression = expression.toarray()

    dfHist = pd.DataFrame(expression, adata.obs['leiden']).reset_index()
    dfHist.columns = ['leiden', 'expression']

    # , log_scale=(False, True)
    
    sns.histplot(data=dfHist, x='expression', hue=colorCol, element="poly", stat='proportion').set(
        xlabel = f'{gene} Expression'
    )

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

