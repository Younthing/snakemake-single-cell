rule run_python_script:
    output:
        "results/hello_python.txt",
    run:
        shell("python ../../scripts/script1.py {output}")
