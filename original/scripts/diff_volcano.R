#!/bin/env Rscript
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
mat <- args[1]
out_mat <- args[2]
vocano <- args[3]
p_name <- ifelse(args[4] == "p", "Pvalue", "Padj")
p_value <- as.numeric(args[5])
logfc_cutoff <- as.numeric(args[6])
top <- args[7]
gene <- args[8]
fc_name <- "logFC"

dif <- data.table::fread(mat, data.table = FALSE)
rownames(dif) <- dif[[1]]
dif_up <- dif %>%
  filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
dif_down <- dif %>%
  filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)
dif <- dif %>% mutate(group = case_when(
  dif[[1]] %in% dif_up[[1]] ~ "Up",
  dif[[1]] %in% dif_down[[1]] ~ "Down",
  TRUE ~ "Not"
))

dif$logP <- -log10(dif[[p_name]])
dif$change <- factor(dif$group,
  levels = c("Down", "Not", "Up")
)
dif <- dif[order(abs(dif$logFC), decreasing = TRUE), ]


if (gene == "NULL" || gene == "None") {
  if (top == "NULL" || top == "None") {
    top <- NULL
  } else {
    top <- as.integer(top)
    gene <- head(dif[[1]], top)
  }
} else {
  gene <- strsplit(gene, ",")[[1]]
}


lab_y <- paste0("-Log10(", p_name, ")")
x_lim <- max(abs(dif$logFC)) * 1.1
p <- ggplot(dif, aes(logFC, logP, color = change)) + # nolint
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = "LogFC", y = lab_y, color = "Significance") +
  geom_hline(yintercept = -log10(p_value), lty = 2) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 2) +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +
  scale_color_manual(values = c(
    "Down" = "#4DBBD5",
    "Not" = "grey", "Up" = "#E64B35"
  )) +
  ggtitle("Volcano") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12)
  ) +
  theme(panel.grid = element_blank())

if (!is.null(gene)) {
  dif2 <- dif
  dif2$label <- mapply(function(x) {
    ifelse(x %in% gene, x, "")
  }, dif2[[1]])

  p <- p + ggrepel::geom_text_repel(
    data = dif2, aes(label = label),
    segment.alpha = 0.4, # 连线透明度
    box.padding = 0.5,
    force = 1,
    max.overlaps = Inf,
    min.segment.length = 0.25,
    show.legend = FALSE
  )
}

ggsave(
  file = vocano,
  p, width = 8, height = 5
)
write.csv(dif, out_mat, row.names = FALSE)
