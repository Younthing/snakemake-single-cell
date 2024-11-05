import os
import sys
import warnings

import decoupler as dc
import liana as li
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from liana.method import cellphonedb, rank_aggregate

# Ignore warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

# Get parameters from Snakemake
unique_prefix = snakemake.params.unique_prefix
input_file = snakemake.input[0]
output_file = snakemake.output[0]
figure_dir = snakemake.params.figure_dir
table_dir = snakemake.params.table_dir
log_file = snakemake.log[0]

# Set log output
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

# 特殊参数

cell_type_col = snakemake.params.cell_type
groups_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label
STIM = treatment_label


def analyze_cell_communication(
    adata,
    group_filter,
    anno_col,
    target_organism=9606,
):
    """
    Perform cell communication analysis using CellPhoneDB and LIANA

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    group_filter : str
        Group to filter for analysis
    anno_col : str
        Column name for cell type annotations
    target_organism : int
        NCBI Taxonomy ID for the target organism
    """
    # Filter data
    adata_stim = adata[adata.obs["group"] == group_filter].copy()
    adata_stim.obs[anno_col] = adata_stim.obs[anno_col].cat.remove_unused_categories()
    adata_stim = adata_stim[adata_stim.obs[anno_col].notnull()].copy()

    # Run CellPhoneDB analysis
    cellphonedb(
        adata_stim,
        groupby=anno_col,
        use_raw=False,
        return_all_lrs=True,
        verbose=True,
        # resource=resource, # TODO notebook逻辑
    )

    # Run LIANA rank aggregate analysis
    rank_aggregate(
        adata_stim,
        groupby=anno_col,
        return_all_lrs=True,
        use_raw=False,
        verbose=True,
        # resource=resource,
    )

    return adata_stim


def plot_communication_results(adata_stim, output_prefix):
    """
    Generate and save visualization plots

    Parameters:
    -----------
    adata_stim : AnnData
        Analyzed data with communication results
    output_prefix : str
        Prefix for output files
    """
    # CellPhoneDB dotplot
    source_labels = adata_stim.uns["liana_res"]["source"].unique()
    target_labels = adata_stim.uns["liana_res"]["target"].unique()

    dotplot_1 = li.pl.dotplot(
        adata=adata_stim,
        colour="lr_means",
        size="cellphone_pvals",
        inverse_size=True,
        source_labels=source_labels,
        target_labels=target_labels,
        filterby="cellphone_pvals",
        filter_lambda=lambda x: x <= 0.01,
        orderby="lr_means",
        orderby_ascending=False,
        top_n=20,
        figure_size=(17, 8),
        size_range=(1, 5),
    )
    dotplot_1.save(f"{figure_dir}/14-{output_prefix}_{STIM}-cellphonedb-dotplot.pdf")

    # LIANA dotplot
    dotplot_liana = li.pl.dotplot(
        adata=adata_stim,
        colour="magnitude_rank",
        size="specificity_rank",
        inverse_colour=True,
        inverse_size=True,
        source_labels=source_labels,
        target_labels=target_labels,
        filterby="specificity_rank",
        filter_lambda=lambda x: x <= 0.01,
        orderby="magnitude_rank",
        orderby_ascending=True,
        top_n=20,
        figure_size=(17, 8),
        size_range=(1, 5),
    )
    dotplot_liana.save(f"{figure_dir}/14-{output_prefix}_{STIM}-liana-dotplot.pdf")


def main():
    # Setup
    # Load data
    adata = sc.read(input_file)
    adata.X = adata.layers["log1p_norm"]

    # Run analysis
    adata_stim = analyze_cell_communication(
        adata,
        group_filter=treatment_label,
        anno_col=cell_type_col,
        target_organism=9606,
    )

    # Generate plots
    plot_communication_results(adata_stim, unique_prefix)

    # Save results
    adata_stim.uns["liana_res"].to_csv(
        f"{table_dir}/14-{unique_prefix}_cellphonedb_res.csv", index=False
    )
    adata_stim.uns["liana_res"].to_csv(
        f"{table_dir}/14-{unique_prefix}_liana_res.csv", index=False
    )

    # Save processed AnnData
    adata_stim.write(output_file)


main()
