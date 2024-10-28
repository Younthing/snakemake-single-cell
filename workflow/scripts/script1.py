def main(param, output_file):
    with open(output_file, "w") as f:
        f.write(f"Python parameter: {param}\n")


if __name__ == "__main__":
    main(snakemake.params.param, snakemake.output[0])
