rule process_data:
    output:
        "results/python_output.txt",
    params:
        param="Hello from Python",
    script:
        "../scripts/script1.py"


rule analyze_data:
    output:
        "results/r_output.txt",
    params:
        param="Hello from R",
    script:
        "../scripts/script2.R"


rule process_data_wrapper:
    output:
        "results/example/wrapper.txt",
    params:
        param="Hello from wrapper",
    log:
        "logs/tool1.log",
    wrapper:
        "file:wrapper/bio/tool1"


# 或者使用远程包装器:
# wrapper:
#     "0.70.0/bio/tool1"
