rule run_python_script:
    output:
        "results/hello_python.txt",
    script:
        "../../scripts/script1.py"
