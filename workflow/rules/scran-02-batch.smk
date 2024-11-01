from pathlib import Path
from typing import Union, List


configfile: "config/config.yaml"


# 全局通用函数
def construct_path(base_dir: Union[str, Path], *args: str) -> str:
    """构建文件路径"""
    return str(Path(base_dir).joinpath(*args))


def get_output_path(sample: str, *args: str) -> str:
    """构建样本的输出路径"""
    return construct_path(config["output_dir"], sample, *args)


def get_log_path(sample: str, *args: str) -> str:
    """构建日志路径"""
    return construct_path(LOG_DIR, sample, *args)


def get_benchmark_path(sample: str, *args: str) -> str:
    """构建基准测试路径"""
    return construct_path(BENCHMARK_DIR, sample, *args)


def create_directories(dirs: List[Union[str, Path]]) -> None:
    """创建多个目录"""
    for directory in dirs:
        Path(directory).mkdir(parents=True, exist_ok=True)


def get_sample_label(sample, label_type):
    """根据样本获取特定标签"""
    try:
        return config["samples"][sample][label_type]
    except KeyError:
        raise ValueError(
            f"Label '{label_type}' for sample '{sample}' is not specified in config."
        )


# 全局常量
OUTPUT_DIR = config["output_dir"]
LOG_DIR = "logs"
BENCHMARK_DIR = "benchmarks"

# 从配置中读取输入文件列表并映射样本
INPUT_FILE_LIST = config["input_file_list"]
SAMPLE_MAP = {Path(input_file).stem: input_file for input_file in INPUT_FILE_LIST}
SAMPLES = list(SAMPLE_MAP.keys())

# 为每个样本创建独立的输出、日志和基准目录
for sample in SAMPLES:
    sample_output_dir = get_output_path(sample)
    sample_log_dir = get_log_path(sample)
    sample_benchmark_dir = get_benchmark_path(sample)
    create_directories([sample_output_dir, sample_log_dir, sample_benchmark_dir])


# 修改后的规则
rule quality_control:
    input:
        adata=lambda wildcards: SAMPLE_MAP.get(wildcards.sample),
    output:
        adata=get_output_path("{sample}", "anndata_qc.h5ad"),
    params:
        mad_general=config["quality_control"]["mad_general"],
        mad_mt=config["quality_control"]["mad_mt"],
        mt_threshold=config["quality_control"]["mt_threshold"],
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
    log:
        get_log_path("{sample}", "quality_control.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "quality_control.txt")
    script:
        "../scripts/single_cell_01_qc.py"


rule doublet_removal:
    input:
        adata=rules.quality_control.output.adata,
    output:
        adata=get_output_path("{sample}", "anndata_doublet_removed.h5ad"),
    params:
        random_state=config["doublet_removal"]["random_state"],
        batch_key=config["doublet_removal"]["batch_key"],
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
    log:
        get_log_path("{sample}", "doublet_removal.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "doublet_removal.txt")
    script:
        "../scripts/single_cell_02_doublet.py"


# 标准化规则
rule normalization:
    input:
        adata=rules.doublet_removal.output.adata,
    output:
        adata=get_output_path("{sample}", "anndata_normalized.h5ad"),
    params:
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
    log:
        get_log_path("{sample}", "normalization.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "normalization.txt")
    script:
        "../scripts/single_cell_03_normalization.py"


# 寻找高变基因规则
rule find_hvg:
    input:
        adata=rules.normalization.output.adata,
    output:
        adata=get_output_path("{sample}", "anndata_hvg.h5ad"),
    params:
        n_top_genes=config["find_hvg"]["n_top_genes"],
        batch_key=config["find_hvg"]["batch_key"],
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
    log:
        get_log_path("{sample}", "find_hvg.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "find_hvg.txt")
    notebook:
        "../notebooks/04-高变基因.ipynb"


# 批次效应去除规则
rule batch_removal:
    input:
        adata=rules.find_hvg.output.adata,
    output:
        adata=get_output_path("{sample}", "anndata_batch_{method}.h5ad"),
    params:
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        batch_key=config["batch_removal"]["batch_key"],
        covariate_keys=config["batch_removal"]["covariate_keys"],
        n_top_pcs=config["batch_removal"]["n_top_pcs"],
        n_neighbors=config["batch_removal"]["n_neighbors"],
        script_path=lambda wildcards: f"../scripts/single_cell_05_batch_removal_{wildcards.method}.py",
    log:
        get_log_path("{sample}", "batch_removal_{method}.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "batch_removal_{method}.txt")
    script:
        "../scripts/single_cell_05_batch_removal_{wildcards.method}.py"


# 聚类规则
rule cluster:
    input:
        adata=rules.batch_removal.output.adata,
    output:
        adata=get_output_path("{sample}", "anndata_cluster_{method}.h5ad"),
    params:
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        method=lambda wildcards: wildcards.method,
        n_neighbors=config["cluster"]["leiden"]["n_neighbors"],
        resolutions=config["cluster"]["leiden"]["resolutions"],
        species=config["cluster"]["species"],
        plot_params=config["cluster"]["plot_params"],
    log:
        get_log_path("{sample}", "cluster_{method}.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "cluster_{method}.txt")
    script:
        "../scripts/single_cell_06_cluster.py"


# 细胞注释规则
rule cell_annotation:
    input:
        adata=rules.cluster.output.adata,
    output:
        adata=get_output_path(
            "{sample}", "anndata_annotation_{method}_{anno_type}.h5ad"
        ),
    params:
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        anno_type=lambda wildcards: wildcards.anno_type,
        method=lambda wildcards: wildcards.anno_type,
        plot_params=config["cell_annotation"]["plot_params"],
    log:
        get_log_path("{sample}", "cell_annotation_{method}_{anno_type}.log"),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path("{sample}", "cell_annotation_{method}_{anno_type}.txt")
    script:
        "../scripts/single_cell_07_annotation.py"


# 找marker基因规则
rule find_markers:
    input:
        adata=rules.cell_annotation.output.adata,
    output:
        adata=get_output_path(
            "{sample}",
            "find_markers",
            "anndata_find_markers_{group_by}_{method}_{anno_type}.h5ad",
        ),
    params:
        unique_prefix=lambda wildcards: "_".join(
            str(wildcards[name]) for name in wildcards.keys()
        ),
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        plot_params=config["cell_annotation"]["plot_params"],
        test_method=config["find_markers"]["test_methods"],
        group_by=lambda wildcards: wildcards.group_by,
    log:
        get_log_path(
            "{sample}",
            "find_markers_{method}_{anno_type}_{group_by}.log",
        ),
    conda:
        config["env"]["conda"]
    benchmark:
        get_benchmark_path(
            "{sample}",
            "find_markers_{method}_{anno_type}_{group_by}.txt",
        )
    script:
        "../scripts/single_cell_08_find_markers.py"


# 差异丰度分析规则
rule differential_abundance:
    input:
        adata=rules.cell_annotation.output.adata,
    output:
        adata=get_output_path(
            "{sample}",
            "differential_abundance",
            "anndata_differential_abundance_{method}_{anno_type}.h5ad",
        ),
    params:
        unique_prefix=lambda wildcards: "_".join(
            str(wildcards[name]) for name in wildcards.keys()
        ),
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        plot_params=config["cell_annotation"]["plot_params"],
        cell_type=lambda wildcards: wildcards.anno_type,
        condition_column=config["condition_column"],
        sample_group_column=config["sample_group_column"],
        ctrol_label=lambda wildcards: config["samples"][wildcards.sample][
            "control_label"
        ],
        treat_label=lambda wildcards: config["samples"][wildcards.sample][
            "treatment_label"
        ],
    log:
        get_log_path(
            "{sample}",
            "differential_abundance_{method}_{anno_type}.log",
        ),
    conda:
        config["env"]["pertpy"]
    benchmark:
        get_benchmark_path(
            "{sample}",
            "differential_abundance_{method}_{anno_type}.txt",
        )
    script:
        "../scripts/single_cell_09_differential_abundance.py"


# Augur优先级分析规则
rule augur_prioritization:
    input:
        adata=rules.cell_annotation.output.adata,
    output:
        adata=get_output_path(
            "{sample}", "augur", "anndata_augur_{method}_{anno_type}.h5ad"
        ),
    params:
        unique_prefix=lambda wildcards: "_".join(
            str(wildcards[name]) for name in wildcards.keys()
        ),
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        plot_params=config["cell_annotation"]["plot_params"],
        cell_type=lambda wildcards: wildcards.anno_type,
        condition_column=config["condition_column"],
        ctrol_label=lambda wildcards: config["samples"][wildcards.sample][
            "control_label"
        ],
        treat_label=lambda wildcards: config["samples"][wildcards.sample][
            "treatment_label"
        ],
    log:
        get_log_path("{sample}", "augur_{method}_{anno_type}.log"),
    conda:
        config["env"]["pertpy"]
    benchmark:
        get_benchmark_path("{sample}", "augur_{method}_{anno_type}.txt")
    script:
        "../scripts/single_cell_10_augur.py"


# 伪批量差异分析规则
rule differential_expression:
    input:
        adata=rules.cell_annotation.output.adata,
    output:
        adata=get_output_path(
            "{sample}",
            "differential_expression",
            "anndata_differential_expression_{method}_{anno_type}.h5ad",
        ),
    params:
        unique_prefix=lambda wildcards: "_".join(
            str(wildcards[name]) for name in wildcards.keys()
        ),
        figure_dir=get_output_path("{sample}", config["figure_dir"]),
        table_dir=get_output_path("{sample}", config["table_dir"]),
        plot_params=config["cell_annotation"]["plot_params"],
        cell_type=lambda wildcards: wildcards.anno_type,
        condition_column=config["condition_column"],
        sample_group_column=config["sample_group_column"],
        ctrol_label=lambda wildcards: config["samples"][wildcards.sample][
            "control_label"
        ],
        treat_label=lambda wildcards: config["samples"][wildcards.sample][
            "treatment_label"
        ],
    log:
        get_log_path(
            "{sample}",
            "differential_expression_{method}_{anno_type}.log",
        ),
    conda:
        config["env"]["pertpy"]
    benchmark:
        get_benchmark_path(
            "{sample}",
            "differential_expression_{method}_{anno_type}.txt",
        )
    script:
        "../scripts/single_cell_11_differential_expression.py"


# 使用expand生成所有样本的所有组合# 使用 expand 生成所有样本的所有组合
BATCH_METHODS = config["batch_removal"]["methods"]
ANNO_TYPES = config["cell_annotation"]["anno_type"]


rule all:
    input:
        # doublet removal 结果
        expand("results/{sample}/anndata_doublet_removed.h5ad", sample=SAMPLES),
        # normalization 结果
        expand("results/{sample}/anndata_normalized.h5ad", sample=SAMPLES),
        # 高变基因结果
        expand("results/{sample}/anndata_hvg.h5ad", sample=SAMPLES),
        # 批次效应去除的不同方法
        expand(
            "results/{sample}/anndata_batch_{method}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
        ),
        # 聚类的不同方法
        expand(
            "results/{sample}/anndata_cluster_{method}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
        ),
        # 细胞注释的不同类型和方法
        expand(
            "results/{sample}/anndata_annotation_{method}_{anno_type}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
            anno_type=ANNO_TYPES,
        ),
        # find_markers 的不同 group_by 和方法
        expand(
            "results/{sample}/find_markers/anndata_find_markers_{group_by}_{method}_{anno_type}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
            anno_type=ANNO_TYPES,
            group_by=config["cell_annotation"]["anno_type"] + ["leiden_0_25"],
        ),
        # 差异丰度分析的不同方法和注释类型
        expand(
            "results/{sample}/differential_abundance/anndata_differential_abundance_{method}_{anno_type}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
            anno_type=ANNO_TYPES,
        ),
        # Augur 优先级分析的不同方法和注释类型
        expand(
            "results/{sample}/augur/anndata_augur_{method}_{anno_type}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
            anno_type=ANNO_TYPES,
        ),
        expand(
            "results/{sample}/differential_expression/anndata_differential_expression_{method}_{anno_type}.h5ad",
            sample=SAMPLES,
            method=BATCH_METHODS,
            anno_type=ANNO_TYPES,
        ),
