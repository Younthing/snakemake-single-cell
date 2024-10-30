"""单细胞转录组数据分析工作流

本工作流实现了单细胞转录组数据的标准分析流程，包含以下步骤：
1. 质量控制(QC):
   - 过滤低质量细胞和基因
   - 计算并过滤线粒体基因比例
2. 双细胞去除:
   - 使用 DoubletFinder 识别和去除双细胞
3. 数据标准化:
   - 对基因表达数据进行标准化处理
4. 识别高变基因(HVG):
   - 识别在细胞间表达变化较大的基因
   - 用于后续降维分析

使用方法：
1. 配置 `config/config.yaml` 文件
2. 运行：`snakemake --cores all`
"""

from pathlib import Path
from typing import Union, List


configfile: "config/config.yaml"


# 全局通用函数
def construct_path(base_dir: Union[str, Path], *args: str) -> str:
    """构建文件路径

    Args:
        base_dir: 基础目录路径
        *args: 子路径组件

    Returns:
        str: 构建好的完整路径字符串
    """
    return str(Path(base_dir).joinpath(*args))


def get_output_path(*args: str) -> str:
    """构建输出路径

    Args:
        *args: 路径组件

    Returns:
        str: 输出路径字符串
    """
    return construct_path(config["output_dir"], *args)


def get_input_path(*args: str) -> str:
    """构建输入路径

    Args:
        *args: 路径组件

    Returns:
        str: 输入路径字符串
    """
    return construct_path(config["input_dir"], *args)


def get_log_path(*args: str) -> str:
    """构建日志路径

    Args:
        *args: 路径组件

    Returns:
        str: 日志路径字符串
    """
    return construct_path(LOG_DIR, *args)


def get_benchmark_path(*args: str) -> str:
    """构建基准测试路径

    Args:
        *args: 路径组件

    Returns:
        str: 基准测试路径字符串
    """
    return construct_path(BENCHMARK_DIR, *args)


def create_directories(dirs: List[Union[str, Path]]) -> None:
    """创建多个目录

    Args:
        dirs: 需要创建的目录列表
    """
    for directory in dirs:
        Path(directory).mkdir(parents=True, exist_ok=True)


# 全局常量
OUTPUT_DIR = config["output_dir"]
FIGURE_DIR = get_output_path(config["figure_dir"])
TABLE_DIR = get_output_path(config["table_dir"])
LOG_DIR = "logs"
BENCHMARK_DIR = "benchmarks"

# 创建必要的目录
create_directories([OUTPUT_DIR, FIGURE_DIR, TABLE_DIR, LOG_DIR, BENCHMARK_DIR])


# QC规则
rule quality_control:
    input:
        adata=get_input_path(config["quality_control"]["input_file"]),
    output:
        adata=get_output_path("anndata_qc.h5ad"),
    params:
        mad_general=config["quality_control"]["mad_general"],
        mad_mt=config["quality_control"]["mad_mt"],
        mt_threshold=config["quality_control"]["mt_threshold"],
        figure_dir=FIGURE_DIR,
        table_dir=TABLE_DIR,
    log:
        get_log_path("quality_control.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("quality_control.txt")
    script:
        "../scripts/single_cell_01_qc.py"


# 双细胞去除规则
rule doublet_removal:
    input:
        adata=rules.quality_control.output.adata,
    output:
        adata=get_output_path("anndata_doublet_removed.h5ad"),
    params:
        random_state=config["doublet_removal"]["random_state"],
        batch_key=config["doublet_removal"]["batch_key"],
        figure_dir=FIGURE_DIR,
        table_dir=TABLE_DIR,
    log:
        get_log_path("doublet_removal.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("doublet_removal.txt")
    script:
        "../scripts/single_cell_02_doublet.py"


# 标准化规则
rule normalization:
    input:
        adata=rules.doublet_removal.output.adata,
    output:
        adata=get_output_path("anndata_normalized.h5ad"),
    params:
        figure_dir=FIGURE_DIR,
        table_dir=TABLE_DIR,
    log:
        get_log_path("normalization.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("normalization.txt")
    script:
        "../scripts/single_cell_03_normalization.py"


# 寻找高变基因规则
rule find_hvg:
    input:
        adata=rules.normalization.output.adata,
    output:
        adata=get_output_path("anndata_hvg.h5ad"),
    params:
        n_top_genes=config["find_hvg"]["n_top_genes"],
        batch_key=config["find_hvg"]["batch_key"],
        figure_dir=FIGURE_DIR,
    log:
        get_log_path("find_hvg.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("find_hvg.txt")
    notebook:
        "../notebooks/04-高变基因.ipynb"


# 批次效应去除规则
rule batch_removal:
    input:
        adata=rules.find_hvg.output.adata,
    output:
        adata=get_output_path("anndata_batch_{method}.h5ad"),
    params:
        figure_dir=FIGURE_DIR,
        batch_key=config["batch_removal"]["batch_key"],
        covariate_keys=config["batch_removal"]["covariate_keys"],
        n_top_pcs=config["batch_removal"]["n_top_pcs"],
        n_neighbors=config["batch_removal"]["n_neighbors"],
        script_path=lambda wildcards: f"../scripts/single_cell_05_batch_removal_{wildcards.method}.py",
    log:
        get_log_path("batch_removal_{method}.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("batch_removal_{method}.txt")
    script:
        "../scripts/single_cell_05_batch_removal_{wildcards.method}.py"  # 不需要使用format,也不能使用lamada


# 批次效应去除方法列表
batch_methods = config["batch_removal"]["methods"]


rule cluster:
    input:
        adata="results/anndata_batch_{method}.h5ad",
    output:
        adata="results/anndata_cluster_{method}.h5ad",
    params:
        figure_dir=FIGURE_DIR,
        table_dir=TABLE_DIR,
        method=lambda wildcards: wildcards.method,
        # 聚类参数
        n_neighbors=config["cluster"]["leiden"]["n_neighbors"],
        resolutions=config["cluster"]["leiden"]["resolutions"],
        # 物种和marker基因
        species=config["cluster"]["species"],
        # 可视化参数
        plot_params=config["cluster"]["plot_params"],
    log:
        "logs/cluster_{method}.log",
    conda:
        config["env"]["conda"]
    benchmark:
        "benchmarks/cluster_{method}.txt"
    script:
        "../scripts/single_cell_06_cluster.py"


# 合并所有批次效应去除方法的聚类结果
rule all_clusters:
    input:
        expand(
            "results/anndata_cluster_{method}.h5ad",
            method=config["batch_removal"]["methods"],
        ),


rule cell_annotation:
    input:
        expand(get_output_path("anndata_batch_{method}.h5ad"), method=batch_methods),
    output:
        expand(get_output_path("anndata_annot_{method}.h5ad"), method=batch_methods),
    params:
        figure_dir=FIGURE_DIR,
