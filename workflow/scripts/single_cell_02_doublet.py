# %%
import os
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

# Suppress warnings
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

# Fixed settings for matplotlib and scanpy
mpl.rcParams["pdf.fonttype"] = 42
sc.settings.verbosity = 0
sc._settings.ScanpyConfig.n_jobs = -1
sc.settings.set_figure_params(
    dpi=80,
    dpi_save=600,
    facecolor="white",
    frameon=False,
)

# Parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]
figure_dir = snakemake.params.figure_dir
table_dir = snakemake.params.table_dir
log_file = snakemake.log[0]
random_state = 123

# Redirect stdout and stderr to log file
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

# Ensure directories exist
os.makedirs(figure_dir, exist_ok=True)
os.makedirs(table_dir, exist_ok=True)

# %%
# Read data
adata = sc.read(input_file)

# %%
# Doublet detection
sc.external.pp.scrublet(
    adata,
    batch_key="batch",
    random_state=random_state,
)

# %%
# Plotting
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

sns.scatterplot(
    data=adata.obs,
    x="n_genes_by_counts",
    y="doublet_score",
    hue="predicted_doublet",
    ax=axes[0],
)
sns.scatterplot(
    data=adata.obs,
    x="total_counts",
    y="doublet_score",
    hue="predicted_doublet",
    ax=axes[1],
)

plt.tight_layout()
plt.savefig(f"{figure_dir}/02-校正-双连体.pdf")
plt.savefig(f"{figure_dir}/02-校正-双连体.tiff", dpi=600)
plt.show()

# %%
# Doublet filtering
print(f"过滤前的细胞总数：{adata.shape[0]}")
print(f"预测为 doublet 的细胞数目：{adata.obs['predicted_doublet'].sum()}")

adata = adata[adata.obs["predicted_doublet"] == False, :]

print(f"过滤后的细胞总数：{adata.shape[0]}")

# %%
# Save the filtered data
adata.write(output_file, compression="gzip")  # type: ignore
