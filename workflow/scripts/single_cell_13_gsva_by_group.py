import os
import sys
import warnings

import decoupler as dc
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc

# 忽略警告
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

# 从 Snakemake 获取参数
unique_prefix = snakemake.params.unique_prefix
input_file = snakemake.input[0]
output_file = snakemake.output[0]
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

# %%
# base param
cell_type_col = snakemake.params.cell_type
groups_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label

# %%
ANNO_COL = cell_type_col
GROUPs_COL = groups_col
CTROL = condition_label
STIM = treatment_label
adata = sc.read(input_file)


heatmap_fig_path = f"{figure_dir}/13-{unique_prefix}-GSVA-heatmap.pdf"
test_csv = f"{table_dir}/12-{unique_prefix}-GSVA-test_cell_group.csv"
test_csv_sig = f"{table_dir}/12-{unique_prefix}-GSVA-test-significant.csv"
sub_fig_dir = figure_dir

adata = sc.read(input_file)

# %%
adata.obs.group.unique()

# %%
""" 
    - 需要原始数据raw,格式也需要处理
"""

import numpy as np
from scipy.sparse import csc_matrix

## 这里只支持csr稀疏矩阵，不支持csc稀疏矩阵
raw = adata.raw.to_adata()
raw.X = adata.raw.X.astype(np.float32).tocsr().toarray()
adata.raw = raw.copy()
adata

# %%
# Retrieving via python
msigdb = dc.get_resource("MSigDB")
msigdb["collection"].unique().tolist()
# Get reactome pathways
hallmark = msigdb.query("collection == 'hallmark'")
# Filter duplicates
hallmark = hallmark[~hallmark.duplicated(("geneset", "genesymbol"))]

# %%
msigdb.head()

# %%
# 这玩意不支持稀疏矩阵
dc.run_gsva(
    mat=adata,
    net=hallmark,
    source="geneset",
    target="genesymbol",
    use_raw=True,
    seed=123,
    kcdf="poisson",  # 重现R的原始行为
    min_n=10,
    mx_diff=True,
    abs_rnk=True,
    # weight="weight",
    verbose=False,
)

# %%
# adata.obsm["gsva_estimate"].to_csv(gsva_csv,index=False)
adata.obsm["gsva_estimate"].head()

# %%
# 这里不使用copy，更改直接反应到adata，不小心重复运行不稳键
estimates = adata.obsm["gsva_estimate"].copy()
# 为列名添加前缀 "path"
estimates.columns = [col for col in estimates.columns]
estimates.columns.to_list()

# %%
## 转换成基因集评分矩阵
acts = dc.get_acts(adata, obsm_key="gsva_estimate")
acts

# %%
acts.obs

# %%
# # 删除包含 "HALLMARK" 的列
# acts.obs = acts.obs.drop(columns=[col for col in acts.obs.columns if "HALLMARK" in col])

# 创建一个新的列用于分组
acts.obs["celltype_group"] = (
    acts.obs[ANNO_COL].astype(str) + "_" + acts.obs[GROUPs_COL].astype(str)
).astype("category")

# 创建一个图形和一个坐标轴
fig, ax = plt.subplots(figsize=(14, 8))

# 绘制矩阵图
sc.pl.matrixplot(
    acts,
    var_names=acts.var_names,
    groupby="celltype_group",  # 使用新的分组列
    # dendrogram=True,  # 如果需要显示树状图
    standard_scale="var",  # 标准化
    colorbar_title="zscores",  # 颜色条标题
    cmap="RdBu_r",
    show=False,
    # swap_axes=True,  # 翻转 x 和 y 轴
    ax=ax,
)

plt.tight_layout()

plt.savefig(heatmap_fig_path, bbox_inches="tight")
plt.show()

# %%
adata.obs[estimates.columns.to_list()] = estimates

# %%
import pandas as pd

# %%
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# 用于存储结果
all_results = []

# 针对每种细胞类型循环
for cell_type in adata.obs[ANNO_COL].unique():
    try:
        # 筛选当前细胞类型的数据
        cell_data = adata.obs[adata.obs[ANNO_COL] == cell_type]

        # 针对每个途径进行检验
        for pathway in estimates.columns:
            group1_scores = cell_data.query(f"{GROUPs_COL} == '{STIM}'")[pathway]
            group2_scores = cell_data.query(f"{GROUPs_COL} == '{CTROL}'")[pathway]
            t_stat, p_val = ttest_ind(group1_scores, group2_scores, nan_policy="omit")

            # 将结果保存到列表
            all_results.append(
                {
                    "Cell_Type": cell_type,
                    "Pathway": pathway,
                    "t_stat": t_stat,
                    "p_val": p_val,
                }
            )
    except Exception as e:
        print(f"Error processing {cell_type}: {e}")

# 转换为 DataFrame
results_df = pd.DataFrame(all_results)
# FDR 校正
results_df["adj_p_val"] = multipletests(results_df["p_val"], method="fdr_bh")[1]

results_df.to_csv(test_csv, index=False)
# 筛选显著结果
significant_results = results_df[results_df["adj_p_val"] < 0.05]
significant_results.to_csv(test_csv_sig, index=False)

# %%
from scipy.stats import ttest_ind


def calculate_t_test(group_a, group_b):
    """计算t检验并返回显著性标志"""
    t_stat, p_val = ttest_ind(group_a, group_b, nan_policy="omit")
    if p_val < 0.001:
        return "***"
    elif p_val < 0.01:
        return "**"
    elif p_val < 0.05:
        return "*"
    else:
        return "ns"


# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


def plot_gsva_violin_box(
    adata,
    estimates,
    cell_type_col,
    group_col,
    stim_group,
    control_group,
    output_dir=None,
    prefix="",
    show_significant_only=False,
    max_pathways=None,
    figsize=(14, 6),
    alpha=0.5,
    p_threshold=0.05,
):
    """
    Plot violin and box plots for GSVA results with significance filtering options.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    estimates : pandas.DataFrame
        GSVA estimates
    cell_type_col : str
        Column name for cell type annotations
    group_col : str
        Column name for group assignments
    stim_group : str
        Name of treatment/stimulation group
    control_group : str
        Name of control group
    output_dir : str, optional
        Directory to save plots
    prefix : str, optional
        Prefix for output filenames
    show_significant_only : bool, optional
        Only show pathways with significant differences (default: False)
    max_pathways : int, optional
        Maximum number of pathways to display
    figsize : tuple, optional
        Figure size (width, height)
    alpha : float, optional
        Transparency for violin plots
    p_threshold : float, optional
        P-value threshold for significance
    """

    def calculate_t_test(group1, group2):
        """Calculate t-test and return significance markers"""
        _, p_value = stats.ttest_ind(group1, group2)
        if p_value >= 0.05:
            return "ns"
        elif p_value >= 0.01:
            return "*"
        elif p_value >= 0.001:
            return "**"
        else:
            return "***"

    # 定义配色方案
    palette = sns.color_palette("husl", len(adata.obs[group_col].unique()))

    # 按细胞类型循环
    for cell_type in adata.obs[cell_type_col].unique():
        try:
            # 筛选当前细胞类型的数据
            cell_data = adata.obs[adata.obs[cell_type_col] == cell_type]

            # 计算每个pathway的显著性和效应量
            pathways = estimates.columns.to_list()
            pathway_stats = []

            for pathway in pathways:
                group1_scores = cell_data.query(f"{group_col} == '{stim_group}'")[
                    pathway
                ]
                group2_scores = cell_data.query(f"{group_col} == '{control_group}'")[
                    pathway
                ]

                # 计算t检验和效应量
                t_stat, p_value = stats.ttest_ind(group1_scores, group2_scores)
                effect_size = abs(np.mean(group1_scores) - np.mean(group2_scores))

                pathway_stats.append(
                    {
                        "pathway": pathway,
                        "p_value": p_value,
                        "effect_size": effect_size,
                        "significant": p_value < p_threshold,
                    }
                )

            # 转换为DataFrame并排序
            stats_df = pd.DataFrame(pathway_stats)
            stats_df = stats_df.sort_values(
                ["significant", "effect_size"], ascending=[False, False]
            )

            # 筛选通路
            if show_significant_only:
                selected_pathways = stats_df[stats_df["significant"]][
                    "pathway"
                ].tolist()
            else:
                selected_pathways = stats_df["pathway"].tolist()

            # 限制通路数量
            if max_pathways is not None and max_pathways > 0:
                selected_pathways = selected_pathways[:max_pathways]

            if not selected_pathways:
                print(f"No significant pathways found for {cell_type}")
                continue

            # 准备绘图数据
            cell_data_melted = cell_data.melt(
                id_vars=[group_col],
                value_vars=selected_pathways,
                var_name="Pathway",
                value_name="gsva_estimate",
            )

            # 创建图形
            fig, ax = plt.subplots(figsize=figsize)

            # 添加小提琴图
            sns.violinplot(
                data=cell_data_melted,
                x="Pathway",
                y="gsva_estimate",
                hue=group_col,
                split=True,
                inner=None,
                palette=palette,
                alpha=alpha,
                ax=ax,
            )

            # 添加箱线图
            sns.boxplot(
                data=cell_data_melted,
                x="Pathway",
                y="gsva_estimate",
                hue=group_col,
                showfliers=False,
                width=0.3,
                palette=palette,
                ax=ax,
            )

            # 避免图例重复
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles[: len(handles) // 2],
                labels[: len(labels) // 2],
                title="Group",
                loc="upper left",
                bbox_to_anchor=(1.01, 1),
            )

            # 添加显著性标记
            y_max = cell_data_melted["gsva_estimate"].max()

            for i, pathway in enumerate(selected_pathways):
                group1_scores = cell_data.query(f"{group_col} == '{stim_group}'")[
                    pathway
                ]
                group2_scores = cell_data.query(f"{group_col} == '{control_group}'")[
                    pathway
                ]

                significance = calculate_t_test(group1_scores, group2_scores)

                x = i
                y = y_max - 0.01 * y_max

                ax.text(
                    x,
                    y,
                    significance,
                    ha="center",
                    color="black",
                    fontsize=12,
                    fontweight="bold",
                )

            # 设置标题
            plt.title(f"{cell_type}")

            # 修复x轴标签对齐
            ax.set_xticks(range(len(selected_pathways)))
            ax.set_xticklabels(selected_pathways, rotation=45, ha="right")

            # 调整布局以确保标签可见
            plt.tight_layout()

            # 微调底部边距，确保旋转的标签完全可见
            plt.subplots_adjust(bottom=0.2)

            # 保存图形
            if output_dir:
                plt.savefig(
                    f"{output_dir}/13-{prefix}-{cell_type}-GSVA-boxplot.pdf",
                    bbox_inches="tight",
                    dpi=600,
                )

            plt.show()
            plt.close()

        except Exception as e:
            print(f"Error processing {cell_type}: {e}")


# %%

# %%
plot_gsva_violin_box(
    adata=adata,
    estimates=estimates,
    cell_type_col=ANNO_COL,
    group_col=GROUPs_COL,
    stim_group=STIM,
    control_group=CTROL,
    output_dir=sub_fig_dir,
    prefix=unique_prefix,
    show_significant_only=True,
    max_pathways=10,
)

# %%


with open(output_file, mode="wt") as f:
    f.write("Analysis completed successfully.")
