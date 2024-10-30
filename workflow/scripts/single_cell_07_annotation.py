import shutil
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

# Set up input and output
input_file = snakemake.input[0]
output_file = snakemake.output[0]
figure_dir = snakemake.params.figure_dir
method = snakemake.params.method
anno_type = snakemake.params.anno_type
plot_params = snakemake.params.plot_params
log_file = snakemake.log[0]
obs_annotation = snakemake.params.method
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

sc.settings.figdir = figure_dir


def setup_scanpy_params(n_jobs: int = -1, verbosity: int = 0) -> None:
    """设置scanpy的基本参数"""
    warnings.filterwarnings("ignore")
    warnings.simplefilter("ignore")
    mpl.rcParams["pdf.fonttype"] = 42
    sc.settings.verbosity = verbosity
    sc.settings.plot_suffix = ""
    sc._settings.ScanpyConfig.n_jobs = n_jobs
    sc.settings.set_figure_params(
        dpi=80,
        dpi_save=600,
        facecolor="white",
        frameon=False,
    )


setup_scanpy_params()

# Load and prepare data
adata = sc.read(input_file)
adata_celltypist = adata.raw.to_adata()

sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)
sc.pp.log1p(adata_celltypist)
adata_celltypist.X = adata_celltypist.X.toarray()

# Download and load celltypist model
import os

import celltypist
from celltypist import models

os.environ["http_proxy"] = ""
os.environ["https_proxy"] = ""

models.download_models(force_update=False, model=["Immune_All_High.pkl"])
model_high = models.Model.load(model="Immune_All_High.pkl")

# Annotate cell types
predictions_low = celltypist.annotate(
    adata_celltypist, model=model_high, majority_voting=True
)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "majority_voting"
]
adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
    adata.obs.index, "conf_score"
]
adata.obs["celltypist_cell_label_fine"] = adata.obs[
    "celltypist_cell_label_fine"
].cat.remove_unused_categories()

# Visualize results

default_plot_params = {
    "vmax": "p90",
    "vmin": "p10",
    "cmap": "bwr",
    "ncols": 3,
    "wspace": 0.2,
    "legend_loc": "on data",
}
if plot_params is None:
    plot_params = {}
plot_params = {**default_plot_params, **plot_params}

sc.pl.umap(
    adata,
    color=["celltypist_cell_label_fine", "celltypist_conf_score_fine"],
    **plot_params,
    save=f"{method}-{anno_type}.pdf",
)

shutil.move(
    f"{figure_dir}/umap{method}-{anno_type}.pdf",
    f"{figure_dir}/07-annotation-{anno_type}-{method}.pdf",
)

print(adata.obs["celltypist_cell_label_fine"].value_counts())

# 结果收集
adata.obs[obs_annotation] = adata.obs["celltypist_cell_label_fine"]
# Save the annotated data
adata.write(output_file, compression="gzip")
