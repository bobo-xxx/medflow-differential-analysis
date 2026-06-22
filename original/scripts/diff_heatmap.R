#!/bin/env Rscript
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(filelock))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))

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
map <- args[2]
rdegs <- args[3]
out_mat <- args[4]
plot <- args[5]
top <- as.integer(args[6])
color_heat <- strsplit(args[7], ",")[[1]]
confirm_file <- args[8]

mat <- data.table::fread(mat, data.table = FALSE)
rownames(mat) <- mat[[1]]
mat <- mat[-1]
map <- data.table::fread(map, data.table = FALSE)
dat_rdegs <- data.table::fread(rdegs, data.table = FALSE)


# 当表型相关差异基因较多时选择logFC_top20展示热图
if (nrow(dat_rdegs) > top * 2) {
  rdegs_sorted <- dat_rdegs[order(dat_rdegs[["logFC"]]), ]
  top_mat <- rdegs_sorted[c(1:top, (nrow(rdegs_sorted) - top + 1):nrow(rdegs_sorted)), ] # nolint
  deg_heatmap <- top_mat[[1]]
} else {
  deg_heatmap <- dat_rdegs[[1]]
}
dat_heatmap <- mat[deg_heatmap, map[[1]]]
map[[2]] <- factor(map[[2]], levels = unique(map[[2]]))
treat_name <- levels(map[[2]])[2]
con_name <- levels(map[[2]])[1]
color <- unique(map[[3]])
names(color) <- unique(map[[2]])
pdf(NULL) # pheatmap 隐式调用 pdf 问题
p_heatmap <- pheatmap(dat_heatmap,
  # 对行归一化（每个基因在样本中的表达量）
  scale = "row",
  # 列注释
  annotation_col = data.frame(
    Group = map[[2]],
    row.names = map[[1]]
  ),
  # 列注释颜色
  annotation_colors = list(Group = color),
  # 热图颜色
  color = colorRampPalette(color_heat)(50),
  # 颜色范围（注意这里的50要和上面括号中的50保持一致）
  breaks = c(seq(-3, 3, length = 50)),
  # 行聚类和列聚类
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = "",
  border_color = NA,
  main = "Heatmap"
)
dev.off()
plain <- 1
names(plain) <- "plain"
main_index <- which(sapply(p_heatmap$gtable$grobs, function(x) x$label) == "Heatmap")
group_index1 <- which(sapply(p_heatmap$gtable$grobs, function(x) x$label) == "Group")
group_index2 <- which(sapply(p_heatmap$gtable$grobs, function(x) x$children[[1]]$label) == "Group")
p_heatmap$gtable$grobs[[main_index]]$gp$font <- plain
p_heatmap$gtable$grobs[[group_index1]]$gp$font <- plain
p_heatmap$gtable$grobs[[group_index2]]$children[[1]]$gp$font <- plain
ggsave(
  file = plot,
  p_heatmap, width = 8, height = 8
)

write.csv(dat_heatmap, out_mat)

file_lock(confirm_file, function(confirm_file) {
  confirm <- if (file.exists(confirm_file)) read_yaml(confirm_file) else list()
  confirm["top"] <- as.character(top)
  write_yaml(confirm, confirm_file)
})
