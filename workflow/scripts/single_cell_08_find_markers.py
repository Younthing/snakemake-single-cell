import shutil
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

# Set up input and output
unique_prefix = snakemake.params.unique_prefix
input_file = snakemake.input[0]
output_file = snakemake.output[0]
figure_dir = snakemake.params.figure_dir
table_dir = snakemake.params.table_dir

plot_params = snakemake.params.plot_params
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

sc.settings.figdir = figure_dir

group_by = snakemake.params.group_by
test_method = snakemake.params.test_method


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
        frameon=False,  # TODO 这个好看些 frameon='small',
    )


setup_scanpy_params()

# Load and prepare data
adata = sc.read(input_file)

# %%
adata.var_names_make_unique()


sc.tl.rank_genes_groups(
    adata,
    groupby=group_by,
    method=test_method,  # 't-test', 'wilcoxon'
    use_raw=False,
    layer="log1p_norm",
)

# %%
sc.pl.rank_genes_groups(adata, ncols=3)

# %%
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=5,
    use_raw=False,
    # vmin=0,
    # vmax=20,
)

# %%
sc.pl.rank_genes_groups_matrixplot(
    adata,
    n_genes=5,
    standard_scale="var",
    use_raw=False,
)

# %%
## 过滤
# min_in_group_fraction: 0.25 min_fold_change: 2, max_out_group_fraction: 0.5
# 为了保留 adata.uns['rank_genes_groups'] 的原始结构，过滤基因设置为NaN。
sc.tl.filter_rank_genes_groups(
    adata,
    min_fold_change=1,  # 最小折叠变化阈值
    min_in_group_fraction=0.25,  # 在组内的最小基因表达比例
    max_out_group_fraction=0.5,  # 在组外的最大基因表达比例
    key="rank_genes_groups",  # 基因组数据的键
    key_added="rank_genes_groups_filtered",  # 过滤后的基因组数据的键
    use_raw=False,
    compare_abs=False,  # 是否比较绝对值
)

# # %%
# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     groupby=group_by,
#     standard_scale="var",
#     n_genes=5,
#     key="rank_genes_groups_filtered",
#     use_raw=False,
# )

# %%
sc.pl.rank_genes_groups_matrixplot(
    adata,
    n_genes=5,
    standard_scale="var",
    key="rank_genes_groups_filtered",
    use_raw=False,
    cmap="Blues",
    save=f"{unique_prefix}.pdf",
)

shutil.move(
    f"{figure_dir}/matrixplot_{unique_prefix}.pdf",
    f"{figure_dir}/08-markers-matrixplot_{unique_prefix}.pdf",
)

# %%
##  保存全部maker列表

deg_table = sc.get.rank_genes_groups_df(
    adata,
    group=adata.obs[group_by].unique(),
    key="rank_genes_groups_filtered",
)
deg_table = deg_table.dropna(subset=["names"])
deg_table.to_csv(
    f"{table_dir}/08-all_filtered_rank_genes_{unique_prefix}.csv", index=False
)

# # %%
# ## 分别保存

for group in adata.obs[group_by].unique():
    rank_genes_df = sc.get.rank_genes_groups_df(
        adata, group=group, key="rank_genes_groups_filtered"
    )
    # 前面的过滤会生成NaN
    rank_genes_df = rank_genes_df.dropna(subset=["names"])

    # 将数据框保存为CSV文件
    rank_genes_df.to_csv(
        f"{table_dir}/08-{unique_prefix}_{group}_rank_genes_.csv", index=False
    )

# # 不转换的话可能无法进行序列化保存
adata.uns["rank_genes_groups_filtered"] = {
    key: str(value) for key, value in adata.uns["rank_genes_groups_filtered"].items()
}
adata.uns["rank_genes_groups"] = {
    key: str(value) for key, value in adata.uns["rank_genes_groups"].items()
}

adata.write_h5ad(output_file)
