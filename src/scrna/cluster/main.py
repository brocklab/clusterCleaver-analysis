import argparse

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import scrna.cfg as cfg
from scrna.utils import get_logger


ROOT = cfg.get_root_dir()
(OUTS := ROOT / "outs" / "cluster").mkdir(exist_ok=True, parents=True)
CFG = cfg.load()

logger = get_logger(__name__, "cluster.log")

cfg.setup_scanpy(OUTS)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset", help="name of dataset to perfrom qc", required=True
    )
    return parser.parse_args()


def load_adata(dataset):
    return sc.read(ROOT / f"data/h5ads/{dataset}-postqc-normalized.h5ad")


def get_highly_variable_genes(adata, dataset):
    logger.info("computing highly variable genes")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, save=f"_{dataset}.svg")


def compute_dimensionality_reductions(adata):
    logger.info("computing pca")
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    # logger.info("computing tsne")
    # sc.tl.tsne(adata, use_rep="X_pca")
    logger.info("computing nearest neighbors")
    sc.pp.neighbors(adata)
    logger.info("computing umap")
    sc.tl.umap(adata)


def plot_dimensionality_reductions(adata, dataset):
    logger.info("plotting dimensionality reductions")
    sc.pl.pca_scatter(adata, color="total_counts", save=f"_counts_{dataset}.svg")
    # sc.pl.tsne(adata, color="total_counts", save=f"_counts_{dataset}.svg")
    sc.pl.umap(
        adata,
        color=["total_counts", "pct_counts_mt"],
        save=f"_counts_pct_mt_{dataset}.png",
    )


def compute_leiden_clusters(adata):
    logger.info("computing leiden clusters")
    sc.tl.leiden(adata)


def main():
    args = get_args()

    adata = load_adata(args.dataset)

    get_highly_variable_genes(adata, args.dataset)

    compute_dimensionality_reductions(adata)
    plot_dimensionality_reductions(adata, args.dataset)
    compute_leiden_clusters(adata)

    outpath = ROOT / f"data/h5ads/{args.dataset}-postqc-normalized-clustered.h5ad"
    logger.info(f"Writing anndata to {outpath}")
    adata.write(outpath)
