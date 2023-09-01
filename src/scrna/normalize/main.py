import argparse

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

import scrna.cfg as cfg
from scrna.utils import get_logger


ROOT = cfg.get_root_dir()
(OUTS := ROOT / "outs" / "normalize").mkdir(exist_ok=True, parents=True)
CFG = cfg.load()

logger = get_logger(__name__, "normalize.log")

cfg.setup_scanpy(OUTS)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset", help="name of dataset to perfrom qc", required=True
    )
    return parser.parse_args()


def normalize_counts(adata):
    """normalize counts with shifted logarithm"""
    logger.info("normalizing counts")

    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)


def load_postqc(dataset):
    return sc.read(ROOT / f"data/h5ads/{dataset}-postqc.h5ad")


def plot_normalized_counts(adata, dataset):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")

    fig.savefig(OUTS / f"normalized_counts_{dataset}.svg")


def main():
    args = get_args()

    adata = load_postqc(args.dataset)

    logger.info("saving raw counts in counts layer")
    adata.layers["counts"] = adata.X.copy()

    logger.info(f"Performing normalization for dataset: {args.dataset}")
    normalize_counts(adata)

    plot_normalized_counts(adata, args.dataset)

    logger.info("updating counts with normalized counts")
    adata.X = adata.layers.pop("log1p_norm")

    logger.info("setting raw annadata")
    adata.raw = adata.copy()

    outpath = ROOT / f"data/h5ads/{args.dataset}-postqc-normalized.h5ad"
    logger.info(f"Writing anndata to {outpath}")

    adata.write(outpath)
