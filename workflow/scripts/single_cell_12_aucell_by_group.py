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
# ANNO_COL = "celltypist"
# SAMPEL_COL = sample_group_column
# GROUPs_COL = groups_col

# CTROL = condition_label
# STIM = treatment_label

# # 定义路径
# IMAGE_DIR = figure_dir
# TABLE_DIR = table_dir

# %%
# ANNO_COL = "celltypist"
# GROUPs_COL =
# CTROL = "CCL"
# STIM = "LCL"
# heatmap_fig_path = "xx.pdf"
# test_csv = "xx.csv"
# test_csv_sig = "xxx-significant.csv"
# sub_fig_dir = "group"

# %%
ANNO_COL = cell_type_col
GROUPs_COL = groups_col
CTROL = condition_label
STIM = treatment_label
adata = sc.read(input_file)


heatmap_fig_path = f"{figure_dir}/12-{unique_prefix}-AUCell-heatmap.pdf"
test_csv = f"{table_dir}/12-{unique_prefix}-AUCell-test_cell_group.csv"
test_csv_sig = f"{table_dir}/12-{unique_prefix}-AUCell-test-significant.csv"
sub_fig_dir = figure_dir

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
""" 
    PROGENy是一个综合资源，包含精选的路径及其目标基因集合，
    以及每个相互作用的权重。对于此示例，我们将使用人类（也可以使用其他生物体），
    并且我们将使用按 p 值排名的前 500 个响应基因。以下是每个途径的简要说明：
---
    雄激素（Androgen）：参与男性生殖器官的生长和发育。
    表皮生长因子受体（EGFR）：调节哺乳动物细胞的生长、存活、迁移、凋亡、增殖和分化。
    雌激素（Estrogen）：促进女性生殖器官的生长和发育。
    缺氧（Hypoxia）：在氧气水平低时促进血管生成和代谢重编程。
    JAK-STAT：参与免疫、细胞分裂、细胞死亡和肿瘤形成。
    MAPK（丝裂原活化蛋白激酶）：整合外部信号，促进细胞生长和增殖。
    NFκB（核因子κB）：调节免疫应答、细胞因子产生和细胞存活。
    p53：调节细胞周期、凋亡、DNA修复和肿瘤抑制。
    PI3K（磷脂酰肌醇3激酶）：促进生长和增殖。
    TGFβ（转化生长因子β）：参与大多数组织的发育、稳态维持和修复。
    TNFα（肿瘤坏死因子α）：介导造血、免疫监视、肿瘤退缩和防止感染。
    Trail：诱导凋亡。
    VEGF（血管内皮生长因子）：介导血管生成、血管通透性和细胞迁移。
    WNT：在发育过程中和组织修复时调节器官形态生成。
"""
progeny = dc.get_progeny(organism="human", top=500)
# progeny = dc.get_progeny(organism="mouse", top=500)

progeny.head()

# %%
dc.run_aucell(
    mat=adata,
    net=progeny,
    source="source",
    target="target",
    use_raw=True,
    seed=123,
    # weight="weight",
    verbose=True,
)
adata.obsm["aucell_estimate"].head()

# %%
# 这里不使用copy，更改直接反应到adata，不小心重复运行不稳键
estimates = adata.obsm["aucell_estimate"].copy()
# 为列名添加前缀 "path"
estimates.columns = ["Pathway-" + col for col in estimates.columns]
estimates.columns

# %%
adata.obs[estimates.columns.to_list()] = estimates

## 查看
[key for key in adata.obs.keys() if "Pathway" in key]
# adata.obs.filter(like="Path")

# %%
# """
#     - 聚类图可视化
# """

# sc.pl.umap(
#     adata[adata.obs[GROUPs_COL]=="LCL"],
#     color=estimates.columns.to_list(),
#     frameon=False,
#     ncols=4,
#     wspace=0.2,

#     cmap="Reds",
#     vmax="p99",
#     vmin="p01",
# )

# %%
acts = dc.get_acts(adata, obsm_key="aucell_estimate")
acts

# %%
import matplotlib.pyplot as plt

# %%
acts.obs["celltype_group"] = (
    acts.obs[ANNO_COL].astype(str) + "_" + acts.obs[GROUPs_COL].astype(str)
).astype("category")
# 创建一个图形和一个坐标轴fig, ax = plt.subplots()

fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.matrixplot(
    acts,
    var_names=acts.var_names,
    groupby="celltype_group",
    num_categories=-1,
    # dendrogram=True,  # 如果需要显示树状图
    standard_scale="var",  # 如果需要标准化
    colorbar_title="zscore",  # 如果需要显示颜色条标题
    cmap="RdBu_r",
    show=False,
    ax=ax,
)
# 获取当前x轴标签

plt.tight_layout()
plt.savefig(heatmap_fig_path, bbox_inches="tight")
plt.show()



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




import matplotlib.pyplot as plt

# %%
import seaborn as sns

# 定义配色方案
palette = sns.color_palette("husl", len(adata.obs[GROUPs_COL].unique()))  # 使用husl色板

# 按细胞类型循环
for cell_type in adata.obs[ANNO_COL].unique():
    try:
        # 筛选当前细胞类型的数据
        cell_data = adata.obs[adata.obs[ANNO_COL] == cell_type]

        # 将数据转换为长格式
        cell_data_melted = cell_data.melt(
            id_vars=[GROUPs_COL],
            value_vars=estimates.columns.to_list(),
            var_name="Pathway",
            value_name="AUCELL Score",
        )

        # 创建图形
        plt.figure(figsize=(14, 6))

        # 添加小提琴图
        sns.violinplot(
            data=cell_data_melted,
            x="Pathway",
            y="AUCELL Score",
            hue=GROUPs_COL,
            split=True,  # 分组显示
            inner=None,  # 不显示小提琴内部的箱线
            palette=palette,
            alpha=0.5,  # 透明度设置
        )

        # 添加箱线图，去掉离群值
        ax = sns.boxplot(
            data=cell_data_melted,
            x="Pathway",
            y="AUCELL Score",
            hue=GROUPs_COL,
            showfliers=False,  # 不显示离群值
            width=0.3,  # 调整箱线宽度
            palette=palette,
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
        pathways = estimates.columns.to_list()
        y_max = cell_data_melted[
            "AUCELL Score"
        ].max()  # 取得Y轴最大值，以便放置显著性标记

        for i, pathway in enumerate(pathways):
            # 获取两组数据
            group1_scores = cell_data.query(f"{GROUPs_COL} == '{STIM}'")[pathway]
            group2_scores = cell_data.query(f"{GROUPs_COL} == '{CTROL}'")[pathway]

            # 计算显著性标记
            significance = calculate_t_test(group1_scores, group2_scores)

            # 显著性标记的 x 坐标为途径位置，y 坐标为图表最大值上方偏移量
            x = i
            y = y_max - 0.01 * y_max  # 调整标记位置

            # 添加显著性标记
            ax.text(
                x,
                y,
                significance,
                ha="center",
                color="black",
                fontsize=12,
                fontweight="bold",
            )
        
           # 修复x轴标签对齐
            ax.set_xticks(range(len(pathways)))
            ax.set_xticklabels(pathways, rotation=45, ha="right")

        # 设置标题和旋转 x 轴标签
        plt.title(f"{cell_type}")
        plt.xticks(rotation=45)
        

        plt.tight_layout()  # 调整布局避免重叠
        plt.savefig(
            f"{sub_fig_dir}/12-{unique_prefix}-{cell_type}-AUC-boxplot.pdf",
            bbox_inches="tight",
            dpi=600,
        )
        plt.show()
    except Exception as e:
        print(f"Error processing {cell_type}: {e}")

with open(output_file, mode="wt") as f:
    f.write("Analysis completed successfully.")
