args <- commandArgs(trailingOnly = TRUE)
param <- args[1]
output_file <- args[2]

writeLines(paste("R parameter:", param), con = output_file)
