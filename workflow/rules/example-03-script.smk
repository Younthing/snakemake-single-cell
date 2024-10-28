rule run_script_py:
    output:
        "results/python_output.txt",
    params:
        param="Hello from Python",
    script:
        "../scripts/script1.py"


rule run_script_R:
    output:
        "results/example/R.txt",  # 输出文件
    params:
        param1="Hello",
        param2="World",
    threads: 1
    script:
        "../scripts/script2.R"


rule run_script_wrapper:
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
