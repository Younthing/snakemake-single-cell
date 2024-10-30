# %%
import os
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

# %%
# 设置参数
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
mpl.rcParams["pdf.fonttype"] = 42  # 保留字体
sc.settings.verbosity = 0  # 4 输出细节
sc.settings.plot_suffix = ""
sc._settings.ScanpyConfig.n_jobs = 96  # 设置并行处理的线程数
sc.settings.set_figure_params(
    dpi=80,
    dpi_save=600,
    facecolor="white",
    frameon=False,  # 移除图框
)

input_file = snakemake.input[0]
output_file = snakemake.output[0]
figure_dir = snakemake.params.figure_dir
batch_key = snakemake.params.batch_key
covariate_keys = snakemake.params.covariate_keys
n_top_pcs = snakemake.params.n_top_pcs
n_neighbors = snakemake.params.n_neighbors
log_file = snakemake.log[0]
sc.settings.figdir = figure_dir

# 重定向标准输出和错误输出到日志文件
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

# 确保目录存在
os.makedirs(figure_dir, exist_ok=True)

# %%
adata_batch = sc.read(input_file)

# %%
## 5.4 Harmony去批次
"""Harmony
    - 稳健性能
    - 无需猜测哪些簇对应于特定细胞类型
    - 这里崩溃可能是openblas，线程数太多不支持
"""
## 5.4.1 预先降维
# 使用反卷积标准化的数据
"""注意文件中的.X的变化,现在是对数化后"""
adata_batch.X = adata_batch.layers["log1p_norm"]

# regressing out 可能会过度校正，见5.5
## sc.pp.regress_out(adata_batch, ["log1p_total_counts_ribo","log1p_total_counts_hb"]) # 可选，修改adata.X，核糖体细胞周期有关
## sc.pp.scale(adata_batch, max_value=10) # 可选，默认修改adata.X，可定义layer
sc.pp.pca(adata_batch, n_comps=50, random_state=123)
sc.pl.pca_variance_ratio(
    adata_batch,
    n_pcs=50,
    # log=True,
    save=".pdf",
)
import shutil

shutil.move(f"{figure_dir}/pca_variance_ratio.pdf", f"{figure_dir}/05-pca-方差占比.pdf")

# %%
## 5.4.2 Harmony
import scanpy.external as sce

# os.environ["OPENBLAS_NUM_THREADS"] = "1"
# os.environ["NUM_THREADS"] = "1"
sce.pp.harmony_integrate(
    adata_batch, key=[batch_key] + covariate_keys, max_iter_kmeans=100
)

# %%
## 5.4.3 整合去批次结果
adata_batch.obsm["X_pca"] = adata_batch.obsm["X_pca_harmony"]
sc.pp.neighbors(adata_batch, n_neighbors=n_neighbors, n_pcs=n_top_pcs, random_state=123)
sc.tl.leiden(adata_batch, random_state=123)
sc.tl.umap(adata_batch, random_state=123)  # tsne可以指定use_rep，umap不可以
# sc.tl.tsne(adata_batch, use_rep="X_pca", random_state=123)

# %%
## 5.4.4 降维聚类可视化批次效应
sc.pl.umap(
    adata_batch,
    color=["leiden"] + [batch_key] + covariate_keys,
    save=".pdf",
)
shutil.move(f"{figure_dir}/umap.pdf", f"{figure_dir}/05-harmony-umap.pdf")

sc.pl.pca(
    adata_batch,
    color=["leiden"] + [batch_key] + covariate_keys,
    save=".pdf",
)
shutil.move(f"{figure_dir}/pca.pdf", f"{figure_dir}/05-harmony-pca.pdf")


# sc.pl.tsne(
#     adata_batch,
#     color=["leiden"] + [batch_key] + covariate_keys,
#     save=".pdf",
# )

# shutil.move(f"{figure_dir}/tsne.pdf", f"{figure_dir}/05-harmony-tsne.pdf")

adata = adata_batch.copy()
adata.write(output_file, compression="gzip")  # type: ignore
