from snakemake.shell import shell

# 获取参数并打印调试信息
output = snakemake.output[0]
param = snakemake.params.get("param", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)  # 修改日志设置

# 打印调试信息到日志
shell(
    """
    echo "Debug: Output path is {output}" >> {snakemake.log[0]}
    echo "Debug: Parameter value is {param}" >> {snakemake.log[0]}
    echo "{param}" > {output} 2>> {snakemake.log[0]}
    """
)
