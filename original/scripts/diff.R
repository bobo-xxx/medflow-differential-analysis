#!/bin/env Rscript
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(filelock))
suppressPackageStartupMessages(library(dplyr))

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

create_file_dir <- function(file) {
  path <- dirname(file)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

## 数据读取和预处理
deal_data <- function(df, group, oirgs) {
  df <- data.table::fread(df, data.table = FALSE)
  rownames(df) <- df[[1]]
  df <- df[, -1]
  group <- data.table::fread(group,
    data.table = FALSE, header = TRUE
  )
  rownames(group) <- group[[1]]
  group[2] <- factor(group[[2]], levels = unique(group[[2]]))
  treat_name <- levels(group[[2]])[2]
  con_name <- levels(group[[2]])[1]
  df <- df[, group[[1]]]
  list(
    mat = df,
    map = group,
    treat_name = treat_name,
    con_name = con_name
  )
}
## deseq2
diff_deseq2 <- function(df, group) {
  suppressPackageStartupMessages(library(DESeq2))
  col_data <- data.frame(
    row.names = colnames(df),
    group_list = group[[2]]
  )
  dds <- DESeqDataSetFromMatrix(
    countData = round(df),
    colData = col_data,
    design = ~group_list
  )
  dds2 <- DESeq(dds)
  res <- results(dds2)
  dif <- res %>%
    as.data.frame() %>%
    dplyr::rename(
      logFC = "log2FoldChange",
      Pvalue = "pvalue", Padj = "padj"
    ) %>%
    na.omit() %>%
    arrange(Pvalue) # nolint
  dif
}
## limma
diff_limma <- function(df, map) {
  pdf(NULL) # limma 隐式调用 plotSA 问题
  suppressPackageStartupMessages(library(limma))
  # treat_name and con_name 这两个变量实际上只保存到了 fit2 中间变量中
  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]
  # 实验设计矩阵
  design <- model.matrix(~ 0 + map[[2]])
  rownames(design) <- map[[1]]
  colnames(design) <- levels(map[[2]])

  # 线性建模
  fit <- lmFit(df, design)
  cont_matrix <- makeContrasts(
    contrasts = paste0(treat_name, "-", con_name),
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont_matrix)

  # 经验贝叶斯调整
  fit2 <- eBayes(fit2)
  plotSA(fit2)

  # 筛选差异基因
  dif <- topTable(fit2, coef = 1, n = Inf)
  dif <- dif %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "P.Value", Padj = "adj.P.Val") %>%
    na.omit() %>%
    arrange(Pvalue) # nolint
  dev.off()
  dif
}
## edgeR
diff_edger <- function(df, map, norm = "TMM", model = "glmFit") {
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(statmod))

  # 数据预处理
  # （1）构建 DGEList 对象
  dgelist <- DGEList(counts = df, group = map[[2]])
  # （2）过滤 low count 数据，例如 CPM 标准化（推荐）
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  # （3）标准化，以 TMM 标准化为例
  dgelist_norm <- calcNormFactors(dgelist, method = norm)
  # 差异表达基因分析
  design <- model.matrix(~ map[[2]])
  # （1）估算基因表达值的离散度
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  # （2）模型拟合，edgeR 提供了多种拟合算法
  func <- get(model)
  fit <- func(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  dif <- lrt %>%
    as.data.frame() %>%
    dplyr::rename(Pvalue = "PValue", Padj = "FDR") %>%
    na.omit() %>%
    arrange(Pvalue) # nolint
  dif
}
## diff stat
### test_row
test_row <- function(row, group, func = t.test) {
  df <- data.frame(x = row, Group = group)
  res <- tryCatch(
    {
      result <- func(x ~ Group, data = df)
      c(result$p.value, result$statistic)
    },
    error = function(e) {
      print(e)
      c(1, 0)
    }
  )
  res
}
### logfc_row
logfc_row <- function(row, treat, control) {
  mean_treatment <- mean(row[treat])
  mean_control <- mean(row[control])
  fc <- mean_treatment / mean_control
  logfc <- log2(fc)
  logfc
}
### diff_stat
diff_stat <- function(df, map, stat) {
  treat_name <- levels(map[[2]])[2]
  con_name <- levels(map[[2]])[1]
  treat <- map[map[[2]] == treat_name, ][[1]]
  control <- map[map[[2]] == con_name, ][[1]]
  res <- data.frame(logFC = apply(df, 1, logfc_row,
    treat = treat, control = control
  ))
  if (stat == "t") {
    df <- apply(df, 1, test_row,
      group = map[[2]], func = t.test
    )
  } else {
    df <- apply(df, 1, test_row,
      group = map[[2]], func = wilcox.test
    )
  }
  df <- t(df)
  colnames(df) <- c("Pvalue", "stat")
  res <- cbind(res, df)
  res["Padj"] <- p.adjust(res[["Pvalue"]], method = "BH")
  dif <- res %>%
    na.omit() %>%
    arrange(Pvalue) # nolint
  dif
}

diff_analysis <- function(mat, map, diff_path, confirm, diff, norm, model) {
  dat <- deal_data(mat, map)
  if (diff == "limma") {
    dif <- diff_limma(dat$mat, dat$map)
  } else if (diff == "deseq2") {
    dif <- diff_deseq2(dat$mat, dat$map)
  } else if (diff == "edgeR") {
    confirm["norm"] <- norm
    confirm["model"] <- model
    dif <- diff_edger(dat$mat, dat$map, norm, model)
  } else {
    dif <- diff_stat(dat$mat, dat$map, diff)
  }
  write.csv(dif, diff_path)

  g <- as.character(unique(dat$map[[2]]))
  if (any(grepl("HighRisk", g, ignore.case = TRUE))) {
    data_set_name <- "{data_set_name_raw}疾病组"
    treat_cn <- "高风险组"
    treat <- "HighRisk"
    control_cn <- "低风险组"
    control <- "LowRisk"
  } else if (any(grepl("Cluster", g, ignore.case = TRUE))) {
    data_set_name <- "{data_set_name_raw}疾病组"
    treat_cn <- "亚组2"
    treat <- "Cluster2"
    control_cn <- "亚组1"
    control <- "Cluster1"
  } else {
    data_set_name <- NULL
    treat_cn <- NULL
    treat <- NULL
    control_cn <- NULL
    control <- NULL
  }
  confirm["data_set_name"] <- data_set_name
  confirm["treat_cn"] <- treat_cn
  confirm["treat"] <- treat
  confirm["control_cn"] <- control_cn
  confirm["control"] <- control
  confirm
}
if (!exists("NAME")) {
  args <- commandArgs(trailingOnly = TRUE)
  mat <- args[1]
  map <- args[2]
  diff_path <- args[3]
  confirm_file <- args[4]
  diff <- args[5]
  norm <- args[6]
  model <- args[7]
  for (file in c(mat, map, confirm_file, diff_path)) {
    create_file_dir(file)
  }
  diif_method <- list(
    limma = paste0("R包limma{ref.limma}（Version ", packageVersion("limma"), "）"),
    deseq2 = paste0("R包DESeq2{ref.deseq2}（Version ", packageVersion("DESeq2"), "）"),
    edgeR = paste0("R包edgeR{ref.edgeR}（Version ", packageVersion("edgeR"), "）"),
    t = "t检验",
    wilcox = "Wilcoxon符号秩检验"
  )
  file_lock(confirm_file, function(confirm_file) {
    confirm <- if (file.exists(confirm_file)) read_yaml(confirm_file) else list()
    confirm <- diff_analysis(mat, map, diff_path, confirm, diff, norm, model)
    confirm["diff_method"] <- diif_method[diff]
    write_yaml(confirm, confirm_file)
  })
}
