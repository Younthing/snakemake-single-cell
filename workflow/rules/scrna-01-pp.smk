configfile: "config/config.yaml"


input_dir = config["input_dir"]
output_dir = config["output_dir"]
figure_dir = f"{output_dir}/{config['figure_dir']}"
table_dir = f"{output_dir}/{config['table_dir']}"


rule quality_control:
    input:
        f"{input_dir}/{config['qc']['input_file']}",
    output:
        f"{output_dir}/{config['qc']['output_file']}",
    params:
        mad_general=config["qc"]["mad_general"],
        mad_mt=config["qc"]["mad_mt"],
        mt_threshold=config["qc"]["mt_threshold"],
        figure_dir=figure_dir,
        table_dir=table_dir,
    conda:
        config["env"]["conda"]
    log:
        f"logs/quality_control.log",
    script:  # 在最后定义
        "../scripts/single_cell_01_qc.py"


rule doublet_removal:
    input:
        f"{output_dir}/{config['doublet_removal']['input_file']}",
    output:
        f"{output_dir}/{config['doublet_removal']['output_file']}",
    params:
        random_state=config["doublet_removal"]["random_state"],
        batch_key=config["doublet_removal"]["batch_key"],
        figure_dir=figure_dir,
        table_dir=table_dir,
    conda:
        config["env"]["conda"]
    log:
        f"logs/doublet_removal.log",
    script:  # 在最后定义脚本
        "../scripts/single_cell_02_doublet.py"


rule normalization:
    input:
        f"{output_dir}/{config['normalization']['input_file']}",
    output:
        f"{output_dir}/{config['normalization']['output_file']}",
    params:
        figure_dir=figure_dir,
        table_dir=table_dir,
    conda:
        config["env"]["conda"]
    log:
        f"logs/normalization.log",
    script:  # 在最后定义脚本
        "../scripts/single_cell_03_normalization.py"
