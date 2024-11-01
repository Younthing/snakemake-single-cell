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

# Redirect stdout and stderr to log file
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

# Ensure directories exist
os.makedirs(figure_dir, exist_ok=True)
os.makedirs(table_dir, exist_ok=True)

# %% [markdown]
# # 正式分析
#
# - 4.标准化
#

# %% [markdown]
# 4.标准化
#
# - 4.1 读取质控后的数据
# - 4.2 标准化前检查
# - 4.3 基本对数标准化
# - 4.4 保存
#
# ---

# %%
## 4.1 读取所有质控后的数据
adata = sc.read(input_file)

# %%
## 4.2 标准化前检查
sns.jointplot(
    adata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex"
)

# %%
## 4.3 基本对数标准化
"""Normalization using sc.pp.normalize_total
   - 基于delta方法的移位对数
   - 优于其他揭示数据集潜在结构的方法（特别是在进行主成分分析时）
   - 并且有利于稳定方差，以进行后续的降维和差异表达基因的识别
"""
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.tight_layout()
plt.savefig(f"{figure_dir}/03-normalization.pdf")
plt.show()

# Save the filtered data
adata.write(output_file, compression="gzip")  # type: ignore
