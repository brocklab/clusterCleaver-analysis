import argparse

import scanpy as sc
import seaborn as sns
import numpy as np
from scipy.stats import median_abs_deviation

import scrna.cfg as cfg
from scrna.utils import get_logger


ROOT = cfg.get_root_dir()
(OUTS := ROOT / "outs" / "qc").mkdir(exist_ok=True, parents=True)
CFG = cfg.load()

cfg.setup_scanpy(OUTS)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset", help="name of dataset to perfrom qc", required=True
    )
    return parser.parse_args()


def load_raw(name):
    return sc.read(ROOT / f"data/h5ads/{name}-raw.h5ad")


def add_var_categories(adata):
    adata.var_names_make_unique()
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))


def plot_qc_plots(adata, name, suffix=""):
    hist = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    hist.savefig(OUTS / f"distplot_counts_{name}{suffix}.svg")

    sc.pl.violin(adata, "pct_counts_mt", save=f"_pct_counts_mt_{name}{suffix}.svg")
    sc.pl.scatter(
        adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        save=f"_counts_{name}{suffix}.svg",
    )


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def main():
    args = get_args()

    logger = get_logger(__name__, "qc.log")

    logger.info(f"Performing QC for datasets: {args.dataset}")

    adata = load_raw(args.dataset)

    logger.info(adata)

    add_var_categories(adata)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    plot_qc_plots(adata, args.dataset)

    logger
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    logger.info(f"Outlier cells by counts: \n{adata.obs.outlier.value_counts()}")
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 20
    )
    logger.info(f"Outlier cells by mt: \n{adata.obs.mt_outlier.value_counts()}")

    logger.info(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    logger.info(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    plot_qc_plots(adata, args.dataset, suffix="_post_qc")

    adata_out = ROOT / f"data/h5ads/{args.dataset}-postqc.h5ad"
    logger.info(f"saving anndata to {adata_out}")
    adata.write(adata_out)
