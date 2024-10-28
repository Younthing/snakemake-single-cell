# 打印输出文件路径
cat("\nOutput files:\n")
for (i in seq_along(snakemake@output)) {
  cat(sprintf("Output %d: %s\n", i, snakemake@output[[i]]))
}

# 打印参数
cat("\nParameters:\n")
cat(sprintf("param1: %s\n", snakemake@params[["param1"]]))
cat(sprintf("param2: %s\n", snakemake@params[["param2"]]))

# 打印线程数
cat(sprintf("\nNumber of threads: %d\n", snakemake@threads))

# 如果需要创建输出文件
writeLines(
  c(
    paste("param1:", snakemake@params[["param1"]]),
    paste("param2:", snakemake@params[["param2"]])
  ),
  snakemake@output[[1]]
)
