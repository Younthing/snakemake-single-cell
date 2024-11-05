# %%
import os
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import omicverse as ov
import pandas as pd
import scanpy as sc

# 忽略警告
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

# 从 Snakemake 获取参数
unique_prefix = snakemake.params.unique_prefix
input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_dir = os.path.dirname(output_file)
figure_dir = snakemake.params.figure_dir
table_dir = snakemake.params.table_dir
log_file = snakemake.log[0]

# 设置日志输出
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

# 设置 Scanpy 参数
sc.settings.figdir = figure_dir


def setup_scanpy_params(n_jobs: int = -1, verbosity: int = 0) -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    sc.settings.verbosity = verbosity
    sc.settings.n_jobs = n_jobs
    sc.settings.set_figure_params(
        dpi=80,
        dpi_save=600,
        facecolor="white",
        frameon=False,
    )


setup_scanpy_params()

# 基础参数
cell_type_col = snakemake.params.cell_type
groups_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label
sub_cell_type = snakemake.params.sub_cell_type


ov.plot_set()

# %%
# 定义您的分组列和细胞类型注释列名称
# unique_prefix=""
# figure_dir = ""
# table_dir = ""
GROUP_COL = groups_col  # 替换为您的实际分组列名，例如 'ctrl' 和 'stim'
SIMPLE_COL = sample_group_column
ANNO_COL = cell_type_col
STIM = treatment_label
CTROL = condition_label


# %%
adata = sc.read(input_file)

# %%
adata.obs[ANNO_COL].unique()

# %%
adata = adata[adata.obs[ANNO_COL] == sub_cell_type].copy()

# %%
adata = ov.pp.preprocess(
    adata,
    mode="shiftlog|pearson",
    n_HVGs=2000,
)

# %%
ov.pp.scale(adata)
ov.pp.pca(adata)

# %%
# ov.pp.umap(adata)

# %% [markdown]
# - 初始化模型

# %%
import numpy as np

## Initialize the cnmf object that will be used to run analyses
cnmf_obj = ov.single.cNMF(
    adata,
    components=np.arange(3, 20),
    n_iter=200,
    seed=123,
    num_highvar_genes=2000,
    output_dir=output_dir + f"/{sub_cell_type}",
    name="dg_cNMF",
)

# %%
## Specify that the jobs are being distributed over a single worker (total_workers=1) and then launch that worker
cnmf_obj.factorize(worker_i=0, total_workers=20)

# %%
cnmf_obj.combine(skip_missing_files=True)

# %% [markdown]
# - Compute the stability and error at each choice of K to see if a clear choice jumps out
cnmf_obj.k_selection_plot(close_fig=False)
# %%
selected_K = 8

# %%
density_threshold = 0.10

# %%

cnmf_obj.consensus(
    k=selected_K,
    density_threshold=density_threshold,
    show_clustering=True,
    close_clustergram_fig=False,
)

import matplotlib.pyplot as plt

# %%
import seaborn as sns
from matplotlib import gridspec

width_ratios = [0.2, 4, 0.5, 10, 1]
height_ratios = [0.2, 4]
fig = plt.figure(figsize=(sum(width_ratios), sum(height_ratios)))
gs = gridspec.GridSpec(
    len(height_ratios),
    len(width_ratios),
    fig,
    0.01,
    0.01,
    0.98,
    0.98,
    height_ratios=height_ratios,
    width_ratios=width_ratios,
    wspace=0,
    hspace=0,
)

D = cnmf_obj.topic_dist[cnmf_obj.spectra_order, :][:, cnmf_obj.spectra_order]
dist_ax = fig.add_subplot(
    gs[1, 1],
    xscale="linear",
    yscale="linear",
    xticks=[],
    yticks=[],
    xlabel="",
    ylabel="",
    frameon=True,
)
dist_im = dist_ax.imshow(
    D, interpolation="none", cmap="viridis", aspect="auto", rasterized=True
)

left_ax = fig.add_subplot(
    gs[1, 0],
    xscale="linear",
    yscale="linear",
    xticks=[],
    yticks=[],
    xlabel="",
    ylabel="",
    frameon=True,
)
left_ax.imshow(
    cnmf_obj.kmeans_cluster_labels.values[cnmf_obj.spectra_order].reshape(-1, 1),
    interpolation="none",
    cmap="Spectral",
    aspect="auto",
    rasterized=True,
)

top_ax = fig.add_subplot(
    gs[0, 1],
    xscale="linear",
    yscale="linear",
    xticks=[],
    yticks=[],
    xlabel="",
    ylabel="",
    frameon=True,
)
top_ax.imshow(
    cnmf_obj.kmeans_cluster_labels.values[cnmf_obj.spectra_order].reshape(1, -1),
    interpolation="none",
    cmap="Spectral",
    aspect="auto",
    rasterized=True,
)

cbar_gs = gridspec.GridSpecFromSubplotSpec(
    3, 3, subplot_spec=gs[1, 2], wspace=0, hspace=0
)
cbar_ax = fig.add_subplot(
    cbar_gs[1, 2],
    xscale="linear",
    yscale="linear",
    xlabel="",
    ylabel="",
    frameon=True,
    title="Euclidean\nDistance",
)
cbar_ax.set_title("Euclidean\nDistance", fontsize=12)
vmin = D.min().min()
vmax = D.max().max()
fig.colorbar(
    dist_im,
    cax=cbar_ax,
    ticks=np.linspace(vmin, vmax, 3),
)
cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), fontsize=12)
plt.tight_layout()  # 确保布局紧凑
plt.savefig(
    f"{output_dir}/16-{unique_prefix}-NMF_cluster.pdf",
    bbox_inches="tight",
)
plt.close()

# %%
density_filter = cnmf_obj.local_density.iloc[:, 0] < density_threshold
fig, hist_ax = plt.subplots(figsize=(4, 4))

# hist_ax = fig.add_subplot(hist_gs[0,0], xscale='linear', yscale='linear',
#   xlabel='', ylabel='', frameon=True, title='Local density histogram')
hist_ax.hist(cnmf_obj.local_density.values, bins=np.linspace(0, 1, 50))
hist_ax.yaxis.tick_right()

xlim = hist_ax.get_xlim()
ylim = hist_ax.get_ylim()
if density_threshold < xlim[1]:
    hist_ax.axvline(density_threshold, linestyle="--", color="k")
    hist_ax.text(
        density_threshold + 0.02, ylim[1] * 0.95, "filtering\nthreshold\n\n", va="top"
    )
hist_ax.set_xlim(xlim)
hist_ax.set_xlabel(
    "Mean distance to k nearest neighbors\n\n%d/%d (%.0f%%) spectra above threshold\nwere removed prior to clustering"
    % (sum(~density_filter), len(density_filter), 100 * (~density_filter).mean())
)
hist_ax.set_title("Local density histogram")

# %%
result_dict = cnmf_obj.load_results(
    K=selected_K,
    density_threshold=density_threshold,
    n_top_genes=None,
)

# %%
result_dict["top_genes"].to_csv(f"{output_dir}/16-{unique_prefix}-NMF-gene_score.csv")
result_dict["top_genes"]

# %%


# %%


# 自带的的索引不是字符串模式 会出错，转换成str
def get_results(adata, result_dict):
    import pandas as pd

    # 将索引转换为字符串以确保匹配
    adata.obs.index = adata.obs.index.astype(str)
    result_dict["usage_norm"].index = result_dict["usage_norm"].index.astype(str)

    if result_dict["usage_norm"].columns[0] in adata.obs.columns:
        # 删除已存在的列
        adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.startswith("cNMF")]

    # 合并 obs
    adata.obs = pd.merge(
        left=adata.obs,
        right=result_dict["usage_norm"],
        how="left",
        left_index=True,
        right_index=True,
    )

    # 合并 var
    adata.var = pd.merge(
        left=adata.var,
        right=result_dict["gep_scores"].loc[adata.var.index],
        how="left",
        left_index=True,
        right_index=True,
    )

    # 计算最大主题
    df = adata.obs[result_dict["usage_norm"].columns].copy()
    max_topic = df.idxmax(axis=1)

    # 将结果添加到 DataFrame 中
    adata.obs["cNMF_cluster"] = max_topic
    print("cNMF_cluster is added to adata.obs")
    print("gene scores are added to adata.var")


# %%

get_results(adata, result_dict)


# %% [markdown]
# - 直接用，细胞类型重叠

# %%
cnmf_obj.get_results_rfc(
    adata, result_dict, use_rep="scaled|original|X_pca", cNMF_threshold=0.5
)

# %%


# %%
plot_genes = []
for i in result_dict["top_genes"].columns:
    plot_genes += result_dict["top_genes"][i][:5].values.reshape(-1).tolist()

# %%
adata.obs["cNMF_cluster"] = adata.obs["cNMF_cluster_rfc"].apply(lambda x: "cNMF_" + x)

# %%

fig, ax = plt.subplots()
ov.pl.embedding(
    adata,
    basis="X_umap",
    color=["cNMF_cluster"],
    frameon="small",
    # title="Celltypes",
    # legend_loc='on data',
    legend_fontsize=14,
    legend_fontoutline=2,
    # size=10,
    # legend_loc=True,
    add_outline=False,
    # add_outline=True,
    outline_color="black",
    outline_width=1,
    ax=ax,
    show=False,
)
plt.savefig(
    f"{figure_dir}/16-{unique_prefix}-cNMF-umap.pdf",
    bbox_inches="tight",
)
plt.close()
# %%
fig, ax = plt.subplots(figsize=(14, 8))
sc.pl.dotplot(
    adata,
    plot_genes,
    "cNMF_cluster",
    # dendrogram=True,
    log=True,
    standard_scale="var",
    show=False,
    ax=ax,
)
fig.savefig(f"{figure_dir}/16-{unique_prefix}-cNMF-dotplot.pdf")

# %%
adata.write(output_file)
