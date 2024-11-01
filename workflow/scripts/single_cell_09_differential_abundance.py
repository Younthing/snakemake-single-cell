# %%
import os
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import pertpy as pt
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
label_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label

# 其他输出文件名
nhood_boxplot_fig = os.path.join(figure_dir, f"09-{unique_prefix}_nhood_boxplot.pdf")

nhood_beeswarmplot_fig = os.path.join(
    figure_dir, f"09-{unique_prefix}_nhood_beeswarmplot.pdf"
)


# %%
# 加载您的数据到 AnnData 对象，这里假设为 'adata'
# 请替换以下代码为您的实际数据加载步骤
adata = sc.read_h5ad(input_file)


# %%
# %%
# 定义分组列和细胞类型注释列名称
GROUP_COL = label_col  # 使用 Snakemake 参数
ANNO_COL = cell_type_col  # 使用 Snakemake 参数
SIMPLE_COL = sample_group_column
GROUP_CONTROL = condition_label  # 使用 Snakemake 参数
GROUP_STIM = treatment_label  # 使用 Snakemake 参数

# 其他中间变量
ALL_CELL_PDF = nhood_boxplot_fig
DA_PDF = nhood_beeswarmplot_fig
# 定义您的分组列和细胞类型注释列名称
# GROUP_COL = "group"  # 替换为您的实际分组列名，例如 'ctrl' 和 'stim'
# ANNO_COL = "celltypist"  # 替换为您的实际细胞类型列名
# SIMPLE_COL = "patients_organ"
# GROUP_CONTROL = "CNL"
# GROUP_STIM = "CCL"
# ALL_CELL_PDF = "results/all_cell_pdf.pdf"
# DA_PDF = "results/da_pdf.pdf"


# %%

# 确保您的分组和细胞类型列存在于 adata.obs 中
if GROUP_COL not in adata.obs.columns or ANNO_COL not in adata.obs.columns:
    raise ValueError(f"确保 {GROUP_COL} 和 {ANNO_COL} 存在于 adata.obs 中。")

# 打印分组和细胞类型信息
print(adata.obs[[GROUP_COL, ANNO_COL]].head())


# %%
# 初始化 Milo 分析对象
milo = pt.tl.Milo()
mdata = milo.load(adata)

# 构建 KNN 图
sc.pp.neighbors(mdata["rna"], use_rep="X_pca", n_neighbors=150)  # 根据需要调整参数

# 构建邻域
milo.make_nhoods(mdata["rna"], prop=0.1)

# %%
mdata = milo.count_nhoods(mdata, sample_col=SIMPLE_COL)

# %%
adata.obs["patients_organ"].unique()

# %%
# (by default, the last category is taken as the condition of interest)
# 这里的样本指的是单个独立分组的样本
mdata["rna"].obs[GROUP_COL] = (
    mdata["rna"].obs[GROUP_COL].cat.reorder_categories([GROUP_CONTROL, GROUP_STIM])
)
milo.da_nhoods(mdata, design=f"~{GROUP_COL}")


# %%
milo.build_nhood_graph(mdata)

# %%
plt.rcParams["figure.figsize"] = [7, 7]
milo.plot_nhood_graph(
    mdata,
    alpha=1,  ## SpatialFDR level (1%)
    min_size=1,  ## Size of smallest dot
    save=f"{unique_prefix}.pdf",
)
import shutil

shutil.move(
    f"{figure_dir}/X_milo_graph{unique_prefix}.pdf",
    f"{figure_dir}/09-{unique_prefix}-milo_gra.pdf",
)
# %%

# 注释邻域（按照主要的细胞类型）
milo.annotate_nhoods(mdata, anno_col=ANNO_COL)


# %%
# Ensure 'nhood_annotation' is a categorical type
mdata["milo"].var["nhood_annotation"] = (
    mdata["milo"].var["nhood_annotation"].astype("category")
)

# Add the new category 'Mixed'
mdata["milo"].var["nhood_annotation"] = (
    mdata["milo"].var["nhood_annotation"].cat.add_categories("Mixed")
)

# Now make the assignment
mdata["milo"].var.loc[
    mdata["milo"].var["nhood_annotation_frac"] < 0.6, "nhood_annotation"
] = "Mixed"

# 去除未使用的类别
mdata["milo"].var["nhood_annotation"] = (
    mdata["milo"].var["nhood_annotation"].cat.remove_unused_categories()
)

import math

# %%
import matplotlib.pyplot as plt

# 假设 cell_types 是一个包含所有细胞类型的列表
cell_types = mdata["milo"].var["nhood_annotation"].unique()

# 动态计算行数和列数
num_types = len(cell_types)
cols = 4
rows = math.ceil(num_types / cols)

# 设置画板大小
plt.figure(figsize=(cols * 4, rows * 4))

for i, cell_type in enumerate(cell_types):
    # 按细胞类型筛选邻域
    cell_nhoods = mdata["milo"].var_names[
        (mdata["milo"].var["nhood_annotation"] == cell_type)
    ]

    plt.subplot(rows, cols, i + 1)
    # 删除点并只绘制箱线图
    milo.plot_nhood_counts_by_cond(
        mdata,
        test_var=GROUP_COL,
        subset_nhoods=cell_nhoods,
        log_counts=False,
    )
    plt.title(f"{cell_type}")

plt.tight_layout()

# 保存图像
plt.savefig(ALL_CELL_PDF, dpi=300)
plt.close()

# %%

plt.figure(figsize=(7, 7))
# 绘制图像
ax = milo.plot_da_beeswarm(mdata, alpha=0.5, palette="Set2")

# 设置标题和标签（如果需要）
ax.set_title("DA Beeswarm Plot")
ax.set_xlabel("logFC")
ax.set_ylabel("nhood_annotation")

# 保存图像
plt.savefig(DA_PDF, dpi=300)

plt.close()
with open(output_file, "wt") as f:
    f.write("Milo analysis completed successfully.")
