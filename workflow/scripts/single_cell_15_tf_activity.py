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

os.environ["http_proxy"] = ""
os.environ["https_proxy"] = ""
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

# 基础参数
cell_type_col = snakemake.params.cell_type
groups_col = snakemake.params.condition_column
sample_group_column = snakemake.params.sample_group_column
condition_label = snakemake.params.ctrol_label
treatment_label = snakemake.params.treat_label

ANNO_COL = cell_type_col
GROUPs_COL = groups_col
CTROL = condition_label
STIM = treatment_label

# 加载数据
adata = sc.read(input_file)

# 检查数据是否成功加载
if adata is None or adata.n_obs == 0:
    raise ValueError(f"Failed to load data from {input_file}")

# 检查必要的列是否存在
required_columns = [GROUPs_COL, ANNO_COL]
for col in required_columns:
    if col not in adata.obs.columns:
        raise ValueError(f"Required column {col} not found in adata.obs")

# 过滤数据
adata = adata[adata.obs[GROUPs_COL] == STIM].copy()

# 检查过滤后是否还有数据
if adata.n_obs == 0:
    raise ValueError(f"No cells remaining after filtering for {STIM}")

# 检索 CollecTRI 基因调控网络
net = dc.get_collectri(organism="human", split_complexes=False)


def run_ulm_inference(adata, net):
    """运行 ULM 推断并存储结果"""
    try:
        dc.run_ulm(
            mat=adata,
            net=net,
            source="source",
            target="target",
            weight="weight",
            verbose=True,
        )
        adata.obsm["collectri_ulm_estimate"] = adata.obsm["ulm_estimate"].copy()
        adata.obsm["collectri_ulm_pvals"] = adata.obsm["ulm_pvals"].copy()
    except Exception as e:
        raise RuntimeError(f"ULM inference failed: {str(e)}")


def get_tf_activities(adata, obsm_key="ulm_estimate"):
    """提取转录因子活性"""
    if obsm_key not in adata.obsm:
        raise KeyError(f"Required key {obsm_key} not found in adata.obsm")
    return dc.get_acts(adata, obsm_key=obsm_key)


def visualize_tf_activities(acts, genes=None, groupby=ANNO_COL):
    """可视化转录因子活性"""
    if genes is None:
        genes = ["PAX5"]

    # 检查指定的基因是否存在
    missing_genes = [g for g in genes if g not in acts.var_names]
    if missing_genes:
        warnings.warn(f"Genes not found in data: {missing_genes}")
        genes = [g for g in genes if g not in missing_genes]

    if not genes:
        raise ValueError("No valid genes to visualize")

    sc.pl.umap(acts, color=genes + [groupby], cmap="RdBu_r", vcenter=0)
    sc.pl.violin(acts, keys=genes, groupby=groupby, rotation=45)


def rank_and_extract_top_markers(acts, groupby=ANNO_COL, n_markers=3):
    """确定每个细胞类型的顶级转录因子"""
    if groupby not in acts.obs.columns:
        raise ValueError(f"Groupby column {groupby} not found in data")

    df = dc.rank_sources_groups(
        acts, groupby=groupby, reference="rest", method="t-test_overestim_var"
    )

    source_markers = (
        df.groupby("group")
        .head(n_markers)
        .groupby("group")["names"]
        .apply(lambda x: list(x))
        .to_dict()
    )
    return source_markers


def plot_tf_network(net, n_sources=None, n_targets=15):
    """绘制网络图"""
    if n_sources is None:
        n_sources = ["PAX5", "EBF1", "RFXAP"]

    # 检查指定的转录因子是否存在于网络中
    missing_tfs = [tf for tf in n_sources if tf not in net["source"].unique()]
    if missing_tfs:
        warnings.warn(f"Transcription factors not found in network: {missing_tfs}")
        n_sources = [tf for tf in n_sources if tf not in missing_tfs]

    if not n_sources:
        raise ValueError("No valid transcription factors to plot")

    dc.plot_network(
        net=net,
        n_sources=n_sources,
        n_targets=n_targets,
        node_size=100,
        s_cmap="white",
        t_cmap="white",
        c_pos_w="darkgreen",
        c_neg_w="darkred",
        figsize=(5, 5),
        save=".pdf",
    )


try:
    # 执行主要分析流程
    run_ulm_inference(adata, net)
    acts = get_tf_activities(adata)
    visualize_tf_activities(acts)
    source_markers = rank_and_extract_top_markers(acts)
    plot_tf_network(net, source_markers)

    # 保存完成状态
    with open(output_file, mode="wt") as f:
        f.write("Analysis completed successfully.")
except Exception as e:
    # 记录错误并确保失败状态被保存
    error_message = f"Analysis failed: {str(e)}"
    with open(output_file, mode="wt") as f:
        f.write(error_message)
    raise RuntimeError(error_message)
finally:
    # 确保日志文件被正确关闭
    sys.stderr.close()
    sys.stdout.close()
