import os
import sys
import warnings

import matplotlib as mpl
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

# Augur 相关参数
# 条件列名称

cell_type_col = snakemake.params.cell_type
label_col = snakemake.params.condition_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label

# 输出文件名
lollipop_fig = os.path.join(figure_dir, f"10-{unique_prefix}_lollipop.pdf")
important_features_fig = os.path.join(
    figure_dir, f"10-{unique_prefix}_important_features.pdf"
)
augur_all_csv = os.path.join(table_dir, f"10-{unique_prefix}_augur_all.csv")
augur_summary_csv = os.path.join(table_dir, f"10-{unique_prefix}_augur_summary.csv")
feature_importances_csv = os.path.join(
    table_dir, f"10-{unique_prefix}_feature_importances.csv"
)
features_importances_csv = os.path.join(
    table_dir, f"10-{unique_prefix}_topn_features_importances.csv"
)

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

# 读入数据
adata = sc.read(input_file)
adata_raw = adata.raw.to_adata()
adata_raw.obs = adata.obs
adata = adata_raw.copy()  # 不copy 就会出现错误


# 创建 Augur 实例
ag_rfc = pt.tl.Augur("random_forest_classifier")

# 加载数据到 Augur 所需的格式

loaded_data = ag_rfc.load(
    adata,
    label_col=label_col,
    cell_type_col=cell_type_col,
    condition_label=condition_label,
    treatment_label=treatment_label,
)
print(loaded_data)

# 运行 Augur 预测
v_adata, v_results = ag_rfc.predict(
    loaded_data,
    subsample_size=50,
    n_threads=40,
    select_variance_features=False,
    span=0.75,
    key_added="augurpy_results",
    random_state=123,
)

# 绘制 lollipop 图
lollipop = ag_rfc.plot_lollipop(
    v_results,
    key="augurpy_results",
    return_fig=True,
)
lollipop.set_size_inches(8, 5)
lollipop.savefig(lollipop_fig, dpi=300, bbox_inches="tight")

# 保存结果到 CSV 文件
v_results["full_results"].to_csv(augur_all_csv, index=True)
v_results["summary_metrics"].to_csv(augur_summary_csv, index=True)

# 计算 UMAP 并绘制 Augur 得分
sc.pp.neighbors(v_adata, random_state=123)
sc.tl.umap(v_adata, random_state=123)

sc.pl.umap(
    adata=v_adata,
    color=["augur_score", "cell_type"],
    wspace=0.2,
    save=f"{unique_prefix}_augur_score.pdf",
    ncols=2,
    cmap="Reds",
)
import shutil

shutil.move(
    f"{figure_dir}/umap{unique_prefix}_augur_score.pdf",
    f"{figure_dir}/10-{unique_prefix}_augur_score.pdf",
)

# 绘制重要特征图
important_features = ag_rfc.plot_important_features(
    v_results, top_n=10, return_fig=True
)
important_features.set_size_inches(8, 5)
important_features.savefig(important_features_fig, dpi=300, bbox_inches="tight")

# 保存重要特征到 CSV 文件
v_results["feature_importances"].to_csv(feature_importances_csv, index=True)

# 计算平均特征重要性并保存
n_features = (
    v_results["feature_importances"]
    .groupby("genes", as_index=False)
    .feature_importances.mean()
    .sort_values(by="feature_importances", ascending=False)
)
n_features.to_csv(features_importances_csv, index=False)

adata.write_h5ad(output_file)
