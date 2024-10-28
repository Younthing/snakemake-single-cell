rule process_hello_world:
    input:
        "results/hello_world.txt",
    output:
        "results/outputs.txt",
    shell:
        """
        cat {input} | tr 'a-z' 'A-Z' > {output}
        """
