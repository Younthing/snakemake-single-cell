import os
import shutil
import sys
import warnings
from typing import Any, Dict, List, Optional, Union

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


def setup_scanpy_params(
    n_jobs: int = 96, verbosity: int = 0, figure_dir: Optional[str] = None
) -> None:
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
        frameon=False,
    )
    if figure_dir:
        sc.settings.figdir = figure_dir


def perform_leiden_clustering(
    adata: ad.AnnData,
    resolutions: List[float],
    n_neighbors: int,
    method,
    random_state: int,
    figure_dir: str,
    marker_genes: Optional[List[str]] = None,
    species: str = "human",
    plot_params: Optional[Dict[str, Any]] = None,
) -> ad.AnnData:
    """
    执行leiden聚类并可视化结果

    参数:
    - adata: AnnData对象
    - resolutions: leiden聚类分辨率列表
    - n_neighbors: 邻域图的邻居数
    - random_state: 随机种子
    - figure_dir: 图片保存目录
    - marker_genes: 用于可视化的标记基因列表
    - species: 物种类型（human/mouse/rat）
    - plot_params: 绘图参数字典，包含vmax, vmin, cmap, ncols等设置
    """
    # 设置默认绘图参数
    default_plot_params = {
        "vmax": "p90",
        "vmin": "p10",
        "cmap": "bwr",
        "ncols": 3,
        "wspace": 0.2,
        "legend_loc": "on data",
    }

    if plot_params is None:
        plot_params = {}
    # 更新默认参数
    plot_params = {**default_plot_params, **plot_params}

    # 确保输出目录存在
    os.makedirs(figure_dir, exist_ok=True)

    # 计算邻域图
    sc.pp.neighbors(
        adata,
        use_rep="X_pca",
        n_neighbors=n_neighbors,
        random_state=random_state,
    )

    # 执行不同分辨率的leiden聚类
    leiden_keys = []
    for res in resolutions:
        # 将分辨率转换为字符串形式的key
        res_str = str(res).replace(".", "_")
        key = f"leiden_{res_str}"
        leiden_keys.append(key)
        sc.tl.leiden(adata, key_added=key, random_state=random_state, resolution=res)

    # 可视化leiden聚类结果
    sc.pl.umap(
        adata,
        color=leiden_keys,
        **plot_params,
        save=f"{method}.pdf",
    )

    # 移动并重命名保存的文件
    shutil.move(
        f"{figure_dir}/umap{method}.pdf", f"{figure_dir}/06-leiden-{method}-umap.pdf"
    )

    # 如果未指定marker基因，使用默认值
    if not marker_genes:
        marker_genes = {
            "human": ["ARG1", "CD33", "HLA-DRA", "CD14", "FCGR3A", "CD68", "CD163"],
            "mouse": ["Ly6c1", "Ly6c2"],
            "rat": ["Cd14", "Fcgr3a", "Cd68", "Cd163"],
        }.get(
            species, ["gene"]
        )  # 默认值不存在返回gene

    # 可视化leiden聚类结果和marker基因表达
    visualization_keys = leiden_keys + marker_genes  # 使用前两个分辨率的聚类结果

    sc.pl.umap(
        adata,
        color=visualization_keys,
        use_raw=True,
        **plot_params,
        save=f"{method}.pdf",
    )

    # 移动并重命名保存的文件
    shutil.move(
        f"{figure_dir}/umap{method}.pdf",
        f"{figure_dir}/06-leiden-{method}-umap-gene.pdf",
    )

    return adata


# 主程序
if __name__ == "__main__":
    # 从snakemake获取所有参数
    input_file: str = snakemake.input[0]
    output_file: str = snakemake.output[0]
    log_file: str = snakemake.log[0]

    # 从params获取参数，设置类型和默认值
    params = snakemake.params
    figure_dir: str = params.figure_dir
    n_neighbors: int = params.n_neighbors
    random_state: int = 123
    n_jobs: int = params.get("n_jobs", -1)
    resolutions: List[float] = params.get("resolutions", [0.25, 0.5, 1, 2])
    species: str = params.get("species", "human")
    method = params.method

    # 获取绘图参数
    plot_params: Dict[str, Any] = params.get("plot_params", {})

    # 设置日志
    sys.stderr = open(log_file, "w")
    sys.stdout = open(log_file, "w")

    # 设置scanpy参数
    setup_scanpy_params(n_jobs=n_jobs, figure_dir=figure_dir)

    # 读取数据
    adata = sc.read(input_file)

    # 执行聚类和可视化
    adata = perform_leiden_clustering(
        adata=adata,
        method=method,
        resolutions=resolutions,
        n_neighbors=n_neighbors,
        random_state=random_state,
        figure_dir=figure_dir,
        species=species,
        plot_params=plot_params,
    )

    # 保存结果
    adata.write(output_file, compression="gzip")
