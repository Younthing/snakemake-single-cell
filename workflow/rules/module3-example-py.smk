rule run_python_script:
    output:
        "results/hello_python.txt",
    script:
        "workflow/scripts/script1.py"
