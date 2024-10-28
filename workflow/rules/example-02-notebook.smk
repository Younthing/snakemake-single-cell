rule run_notebook_papermill:
    input:
        notebook="workflow/notebooks/example_notebook_py.ipynb",
    output:
        notebook_out="results/example_papermill.ipynb",
    conda:
        "lang"
    params:
        param1=42,
        param2="'Hello, Snakemake and Jupyter!'",
        kernel="python3",  # 指定所需的 Jupyter 内核
    shell:
        """
        papermill {input.notebook} {output.notebook_out} -p param1 {params.param1} -p param2 {params.param2} -k {params.kernel}
        """


rule run_notebook_native:
    conda:
        "lang"
    params:
        param1=42,
        param2="Hello, Snakemake and Jupyter!",
    notebook:
        # 注意：对笔记本nbformat版本有要求，不然会提示id错误，vscode默认创建的笔记本版本都不对
        "../notebooks/example_notebook_py_naive.ipynb"  # echo '{"cells":[],"metadata":{},"nbformat":4,"nbformat_minor":5}' > new_notebook.ipynb
