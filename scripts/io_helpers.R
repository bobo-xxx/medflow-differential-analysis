# io_helpers.R — File I/O and color utilities for differential-analysis
#
# Provides advisory file locking, directory creation, and color mapping.

suppressPackageStartupMessages(library(filelock))

#' Acquire an exclusive file lock, run FUN, then release
#'
#' @param path Path to the file to lock
#' @param FUN Function to call with path as first argument
#' @param ... Additional arguments passed to FUN
#' @param exclusive Logical, whether lock is exclusive (default TRUE)
#' @param timeout Timeout in milliseconds (default 5000)
#' @return Result of FUN(path, ...), invisibly
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

#' Create parent directory for a file if it does not exist
#'
#' @param file File path whose parent directory should exist
create_file_dir <- function(file) {
  path <- dirname(file)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#' Map colors to named groups
#'
#' Cycles colors if there are more groups than colors.
#'
#' @param color Character vector of color values
#' @param group Character vector of group names
#' @return Named character vector mapping group -> color
color_map <- function(color, group) {
  color <- as.vector(color)
  group <- as.vector(group)
  g <- unique(group)
  n <- length(g)
  color <- head(rep(color, ceiling(n / length(color))), n)
  names(color) <- g
  color
}
