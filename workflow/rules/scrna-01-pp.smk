configfile: "config/config.yaml"

input_dir = config["input_dir"]
output_dir = config["output_dir"]
adata_raw = config["adata_raw"]


adata_raw = f"{input_dir}/{adata_raw}.h5ad"
figure_dir = f"{output_dir}/figures"
table_dir = f"{output_dir}/tables"

rule run_single_cell_pp:
    input:
        config["adata_raw"],
    output_dir:
        figures        
