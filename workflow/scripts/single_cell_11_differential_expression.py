import os
import sys
import warnings

import decoupler as dc
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
groups_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label


# %%
adata = sc.read_h5ad(input_file)
adata = adata.raw.to_adata()

# %%
# 将原始计数存储在一个层中，因为大多数DGE模型需要原始计数
adata.layers["counts"] = adata.X.copy()

# %%
ANNO_COL = cell_type_col
SAMPEL_COL = sample_group_column
GROUPs_COL = groups_col

CTROL = condition_label
STIM = treatment_label

# 定义路径
IMAGE_DIR = figure_dir
TABLE_DIR = table_dir

# 确保目录存在
os.makedirs(IMAGE_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# %%
# 创建伪批量样本
ps = pt.tl.PseudobulkSpace()
pdata = ps.compute(
    adata,
    target_col=SAMPEL_COL,
    groups_col=ANNO_COL,
    layer_key="counts",
    mode="sum",
    min_cells=10,
    min_counts=300,
)

# %%
# 对伪批量样本进行预处理和PCA分析
pdata.layers["counts"] = pdata.X.copy()
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# %%
# 返回原始计数到X
dc.swap_layer(pdata, "counts", X_layer_key=None, inplace=True)

# %%
sc.pl.pca_variance_ratio(pdata)

# %%
sc.pl.pca(pdata, color=[SAMPEL_COL, ANNO_COL, GROUPs_COL], ncols=1, size=300)

# %%
# 方差分析
dc.get_metadata_associations(
    pdata,
    obs_keys=[
        SAMPEL_COL,
        ANNO_COL,
        GROUPs_COL,
        "psbulk_n_cells",
        "psbulk_counts",
    ],
    obsm_key="X_pca",
    uns_key="pca_anova",
    inplace=True,
)
dc.plot_associations(
    pdata,
    uns_key="pca_anova",
    obsm_key="X_pca",
    stat_col="p_adj",
    obs_annotation_cols=[
        SAMPEL_COL,
        ANNO_COL,
        GROUPs_COL,
    ],
    titles=["Principle component scores", "Adjusted p-values from ANOVA"],
    save=f"{IMAGE_DIR}/11-{unique_prefix}-ANOVA-heatmap.pdf",
    return_fig=True,
)

# %%
# 使用edgeR进行差异基因表达分析
edgr = pt.tl.EdgeR(
    pdata,
    design=f"~{GROUPs_COL}",
)
edgr.fit()

# %%
# 计算对比以确定不同处理之间的差异基因
res_df = edgr.test_contrasts(
    edgr.contrast(column=GROUPs_COL, baseline=CTROL, group_to_compare=STIM)
)
res_df.to_csv(
    f"{TABLE_DIR}/11-{unique_prefix}-differential_expression_EdgeR.csv", index=False
)

# %%
edgr.plot_volcano(
    res_df,
    log2fc_thresh=0,
    save=False,
)
plt.savefig(
    f"{IMAGE_DIR}/11-{unique_prefix}-differential_expression_volcano_EdgeR.pdf"
)  # Save manually


# %%
def analyze_and_plot(pdata, anno_col, groups_col, control, stim):
    """
    Analyze differential expression and plot volcano plots for each cell type.

    Args:
        pdata: The input data, expected to be an AnnData object.
        anno_col: The column in pdata.obs containing cell type annotations.
        groups_col: The column in pdata.obs used for grouping in the design formula.
        control: The baseline group for the contrast.
        stim: The group to compare against the baseline.
    """
    for cell in pdata.obs[anno_col].unique():
        try:
            # Filter data for the specific cell type
            cell_data = pdata[pdata.obs[anno_col] == cell]

            # Perform differential expression analysis using EdgeR
            edgr = pt.tl.EdgeR(
                cell_data,
                design=f"~{groups_col}",
            )
            edgr.fit()

            # Test contrasts
            res_df = edgr.test_contrasts(
                edgr.contrast(
                    column=groups_col, baseline=control, group_to_compare=stim
                )
            )

            # Save differential expression results
            table_filename = os.path.join(
                TABLE_DIR, f"11-{cell}_differential_expression_EdgeR.csv"
            )
            res_df.to_csv(table_filename, index=False)

            # Plot volcano plot and annotate with cell type
            plot_filename = os.path.join(
                IMAGE_DIR, f"11-{unique_prefix}-{cell}_volcano_plot.pdf"
            )
            edgr.plot_volcano(
                res_df,
                log2fc_thresh=0,
                save=False,
            )
            plt.savefig(plot_filename)  # Save manually
        except Exception as e:
            print(f"Error analyzing cell type {cell}: {e}")


# 示例调用
analyze_and_plot(pdata, ANNO_COL, GROUPs_COL, CTROL, STIM)

with open(snakemake.output[0], "w") as f:
    f.write("Analysis completed successfully.")
