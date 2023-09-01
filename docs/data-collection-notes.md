# General Notes From Collecting Data Sources

## Kinker et al.

Raw umi's and cell assignments were downloaded from the [Broad Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity).

Files:

- umi-counts.txt: UMIcount_data.txt
- metadata.txt: Metadata.txt

After downloading the `umi-counts.txt`.
Rows 2,3 were manually removed with the below `sed` command:

```sh
sed -i -e '/Cell_line.*/d' -e '/Pool_ID.*/d' umi-counts.txt
```

## 231-1KB3

This is an internal dataset consisting of 6 samples.

For now the h5 can be copied directly from a local directory on the POD.

```sh
mkdir -p data/raw/231-1KB3
cp /stor/scratch/Brock/daylin/barcoded231/scrna/fastq2mtx-dm2103/aggregate/outs/count/filtered_feature_bc_matrix.h5 data/raw/231-1KB3/
```

