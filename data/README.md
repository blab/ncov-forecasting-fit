# Data preparation

## Preparing metadata for sequence counts

To begin data preparation, we first download the Nextstrain-curated metadata TSV of the GISAID database.
This metadata file is used to construct the datasets in the sub-folders.

1. Nextstrain-curated metadata TSV of GISAID database was downloaded using [nextstrain-cli](https://docs.nextstrain.org/projects/cli/en/stable/). Uncompressing and renaming this file resulted in `gisaid_metadata.tsv` via:
```
nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
gzip -d metadata.tsv.gz -c > gisaid_metadata.tsv
```

2. The metadata file was pruned to only relevant columns using the `tsv-select` command from [tsv-utils](https://opensource.ebay.com/tsv-utils/):
```
tsv-select -H -f strain,date,country,division,QC_overall_status,Nextstrain_clade gisaid_metadata.tsv > gisaid_metadata_pruned.tsv
```

These operations can be found in the script `./download_prune_metadata.sh`.

## Creating data sets as of observation date

This generates the sequence counts for variants of interest from available sequence data as of a given date as well as an additional retrospective truth data set i.e. including all sequences up to now. 

This process can be found in the notebook `./creating-data-sets.ipynb` and generated data sets can be found in `./time_stamped/`

## Creating data sets under down-scaled sequencing efforts

We generate sequence counts for variants of interest as of various observation dates for the United Kingdom under various versions of down-scaled sequencing efforts including increased submission delay and decreasing total sequences per week.

This process for generating data sets can be found in the notebook `./simulating-downscaling-data-sets.ipynb` and generated data sets can be found in `./down_scaled/`.
