#!/bin/env Rscript
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(filelock))
suppressPackageStartupMessages(library(RCircos))

file_lock <- function(path, FUN, ..., exclusive = TRUE, timeout = 5000) {
  FUN <- match.fun(FUN)
  lock_file <- paste0(path, ".lock")
  lock <- lock(lock_file, exclusive = exclusive, timeout = timeout)
  unlock <- unlock
  if (is.null(lock)) {
    stop(paste0("The file lock cannot be obtained: ", lock_file))
  } else {
    res <- tryCatch(
      forceAndCall(1, FUN, path, ...),
      error = function(e) stop(e),
      finally = unlock(lock)
    )
  }
  invisible(res)
}


args <- commandArgs(trailingOnly = TRUE)
mat <- args[1]
locate <- args[2]
plot <- args[3]
tax_id <- args[4]
confirm_file <- args[5]


chr <- read.csv(locate)
gene <- read.csv(mat)[[1]]

chr_gene <- chr[which(chr$Gene %in% gene), ]
pdf(file = plot, width = 8, height = 8)
if (tax_id == "9606") {
  data(UCSC.HG38.Human.CytoBandIdeogram)
  cyto_info <- UCSC.HG38.Human.CytoBandIdeogram
  tax_name <- "人类"
} else if (tax_id == "10090") {
  data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
  cyto_info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
  tax_name <- "小鼠"
} else {
  stop("Unsupported tax_id")
}


RCircos.Set.Core.Components(cyto_info)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Gene.Connector.Plot(chr_gene, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(chr_gene, name.col = 4, track.num = 2, side = "in")
dev.off()
# 结束

file_lock(confirm_file, function(confirm_file) {
  confirm <- if (file.exists(confirm_file)) read_yaml(confirm_file) else list()
  confirm["locate_method"] <- paste0("并使用R包RCircos （Version ", packageVersion("RCircos"), "）绘制染色体定位图展示")
  confirm["locate_result"] <- paste0(
    "最后，通过R包RCircos分析{rdegs}个{pheno_cn}相关差异表达基因（{pheno_abbr}RDEGs）在", tax_name, "染色体上的位置，绘制染色体定位图（{fig.}D）。",
    "染色体定位图展示了{pheno_cn}相关差异表达基因在", tax_name, "染色体上的分布情况。"
  )
  write_yaml(confirm, confirm_file)
})
