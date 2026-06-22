#!/bin/env Rscript
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(ggplot2))

color_map <- function(color, group) {
  color <- as.vector(color)
  group <- as.vector(group)
  g <- unique(group)
  n <- length(g)
  color <- head(rep(color, ceiling(n / length(color))), n)
  names(color) <- g
  color
}

args <- commandArgs(trailingOnly = TRUE)
mat <- args[1]
rgs <- args[2]
venn <- args[3]
pheno_abbr <- args[4]
color_panel <- strsplit(args[5], ",")[[1]]


dif <- data.table::fread(mat, data.table = FALSE)
rgs <- data.table::fread(rgs, data.table = FALSE)[[1]]

dif <- dif[dif[["group"]] == "Up" | dif[["group"]] == "Down", ]
p_venn <- ggvenn(
  setNames(
    list(dif[[1]], rgs),
    c("DEGs", pheno_abbr)
  ),
  # 下面要与列表中的命名一致
  c("DEGs", pheno_abbr),
  # 不展示比例
  show_percentage = FALSE,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = color_map(color_panel, c("DEGs", pheno_abbr))
)
ggsave(file = venn, p_venn, width = 8, height = 5)
