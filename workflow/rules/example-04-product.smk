import os
import itertools


configfile: "config/config_product.yaml"


# 获取样本列表和参数列表
samples = config["samples"]
n_top_genes_list = config["parameters"]["n_top_genes"]
mt_threshold_list = config["parameters"]["mt_threshold"]

# 生成所有参数组合
param_combinations = list(itertools.product(n_top_genes_list, mt_threshold_list))

# 生成所有样本和参数组合的笛卡尔积
all_combinations = list(itertools.product(samples, n_top_genes_list, mt_threshold_list))

# 定义输入和输出目录
input_dir = config["input_dir"]
output_dir = config["output_dir"]


# 生成所有目标输出文件
rule example_all_product:
    input:
        expand(
            os.path.join(
                output_dir,
                "{sample}",
                "n_top_genes_{n_top_genes}",
                "mt_threshold_{mt_threshold}",
                "{sample}_results.txt",
            ),
            sample=[combo[0] for combo in all_combinations],
            n_top_genes=[combo[1] for combo in all_combinations],
            mt_threshold=[combo[2] for combo in all_combinations],
        ),


rule analyze_sample:
    input:
        os.path.join(input_dir, "{sample}.fastq"),
    output:
        os.path.join(
            output_dir,
            "{sample}",
            "n_top_genes_{n_top_genes}",
            "mt_threshold_{mt_threshold}",
            "{sample}_results.txt",
        ),
    params:
        n_top_genes=lambda wc: int(wc.n_top_genes),
        mt_threshold=lambda wc: float(wc.mt_threshold),
    log:
        os.path.join(
            "logs",
            "{sample}",
            "n_top_genes_{n_top_genes}",
            "mt_threshold_{mt_threshold}",
            "{sample}.log",
        ),
    shell:
        """
        echo "Processing {input} with n_top_genes={params.n_top_genes} and mt_threshold={params.mt_threshold}" > {output}
        """
