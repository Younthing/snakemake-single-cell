# %%
import os
import shutil
import sys
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns

# %%
# 设置参数
warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
mpl.rcParams["pdf.fonttype"] = 42  # 保留字体
sc.settings.verbosity = 0  # 4 输出细节
sc.settings.plot_suffix = ""
sc._settings.ScanpyConfig.n_jobs = -1  # 设置并行处理的线程数
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
## 5.4 scvi去批次
adata_batch.X = adata_batch.layers["counts"]

# %%
## 5.4.1 建模去批次
scvi.model.SCVI.setup_anndata(
    adata_batch,
    layer="counts",
    batch_key=batch_key,
    categorical_covariate_keys=covariate_keys,
)

# 建立模型
vae = scvi.model.SCVI(
    adata_batch,
    n_layers=2,
    n_latent=10,
    gene_likelihood="zinb",  #
    dispersion="gene-batch",  #
    dropout_rate=0.3,
)

# 训练
vae.train(
    early_stopping=True,  # TODO 多gpu冲突，需要仔细设置
)

# %%
## 5.4.2 整合去批次结果
adata_batch.obsm["X_scVI"] = vae.get_latent_representation()

adata_batch.obsm["X_pca"] = adata_batch.obsm["X_scVI"]
sc.pp.neighbors(
    adata_batch, use_rep="X_scVI", n_neighbors=n_neighbors, random_state=123
)  # TODO，后面都用scVI，也可以直接把PCA替换成scVI
sc.tl.leiden(adata_batch, random_state=123)
sc.tl.umap(adata_batch, random_state=123)
# sc.tl.tsne(adata_batch, use_rep="X_scVI", random_state=123)

# %%
## 5.4.3 降维聚类可视化批次效应
sc.pl.umap(
    adata_batch,
    color=["leiden", batch_key] + covariate_keys,
    save=".pdf",
)
shutil.move(f"{figure_dir}/umap.pdf", f"{figure_dir}/05-scvi-umap.pdf")

sc.pl.pca(
    adata_batch,
    color=["leiden", batch_key] + covariate_keys,
    save=".pdf",
)
shutil.move(f"{figure_dir}/pca.pdf", f"{figure_dir}/05-scvi-pca.pdf")

# sc.pl.tsne(adata_batch, color=["leiden", batch_key] + covariate_keys, save=".pdf")
# shutil.move(f"{figure_dir}/tsne.pdf", f"{figure_dir}/05-scvi-tsne.pdf")


adata_batch.write(output_file, compression="gzip")  # type: ignore
