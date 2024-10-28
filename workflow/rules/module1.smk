rule hello_world:
    output:
        "results/hello_world.txt",
    shell:
        """
        echo "Hello World" > {output}
        date >> {output}
        """
