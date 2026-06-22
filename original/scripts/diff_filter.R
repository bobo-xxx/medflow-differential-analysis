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


test_cutoff <- function(dif, fc_name, p_name, logfc_test, p_value = 0.05) {
  test <- c()
  for (value in logfc_test) {
    test <- append(
      test,
      sum((abs(dif[fc_name]) > value) & (dif[p_name] < p_value))
    )
  }
  names(test) <- logfc_test
  return(test)
}
args <- commandArgs(trailingOnly = TRUE)
mat <- args[1]
in_gene <- args[2]
out_mat <- args[3]
out_rdegs <- args[4]
confirm_file <- args[5]
p_set <- args[6]
p_value <- as.numeric(args[7])
cutoff <- as.numeric(args[8])
logfc_cutoff <- as.numeric(args[9])
logfc_test <- seq(4, 0, -0.5)
fc_name <- "logFC"

p_name <- ifelse(p_set == "p", "Pvalue", "Padj")
dif <- data.table::fread(mat, data.table = FALSE)
rownames(dif) <- dif[[1]]
test <- test_cutoff(dif, fc_name, p_name, logfc_test, p_value)
num <- sum((abs(dif[fc_name]) > logfc_cutoff) & (dif[p_name] < p_value))
if (num < cutoff) {
  rlang::abort(
    message = paste0(
      "使用 ", p_name, " 过滤 p=", p_value, " 后，基因数小于 ", cutoff
    ),
    class = "value_error"
  )
}

dif_up <- dif %>%
  filter(!!sym(fc_name) > logfc_cutoff & !!sym(p_name) < p_value)
dif_down <- dif %>%
  filter(!!sym(fc_name) < -logfc_cutoff & !!sym(p_name) < p_value)
dif2 <- dif %>% mutate(group = case_when(
  dif[[1]] %in% dif_up[[1]] ~ "Up",
  dif[[1]] %in% dif_down[[1]] ~ "Down",
  TRUE ~ "Not"
))
degs <- dif2[dif2[["group"]] == "Up" | dif2[["group"]] == "Down", ]

if (in_gene != "None") {
  rgs <- data.table::fread(in_gene, data.table = FALSE)[[1]]
  genes <- intersect(degs[[1]], rgs)
  rdegs <- degs[genes, ]
  if (nrow(rdegs) < cutoff) {
    rlang::abort(
      message = paste0(
        "与表型交集的上下调基因数不足 ", cutoff
      ),
      class = "value_error"
    )
  }
} else {
  rdegs <- degs
}

write.csv(rdegs, out_rdegs, row.names = FALSE)
write.csv(degs, file = out_mat, row.names = FALSE)

file_lock(confirm_file, function(confirm_file) {
  confirm <- if (file.exists(confirm_file)) read_yaml(confirm_file) else list()
  confirm["rdegs"] <- nrow(rdegs)
  confirm["test"] <- list(test)
  confirm["diff_p_name"] <- ifelse(p_set == "p", "p value", "p adjust")
  confirm["diff_p_cut"] <- p_value
  confirm["cutoff"] <- cutoff
  confirm["logfc_cut"] <- logfc_cutoff
  confirm["up"] <- nrow(dif_up)
  confirm["down"] <- nrow(dif_down)
  confirm["total"] <- nrow(dif_up) + nrow(dif_down)
  confirm <- confirm[unique(names(confirm))]
  write_yaml(confirm, confirm_file)
})
