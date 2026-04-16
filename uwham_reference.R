#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("UWHAM", quietly = TRUE)) {
    stop("The UWHAM package is required. Install it from CRAN before running this script.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript scripts/uwham_reference.R <input-dir> [base-state-index]")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
base_state <- if (length(args) == 2) {
  as.integer(args[[2]])
} else {
  NA_integer_
}

read_logq <- function(path) {
  as.matrix(read.csv(path, header = FALSE, check.names = FALSE))
}

read_size <- function(path) {
  size_df <- read.csv(path, header = TRUE, check.names = FALSE)
  if (ncol(size_df) != 1L) {
    stop("size.csv must contain exactly one data column")
  }
  as.integer(size_df[[1L]])
}

write_matrix_csv <- function(path, values) {
  write.table(
    values,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

write_vector_csv <- function(path, values, header) {
  write.table(
    data.frame(setNames(list(values), header)),
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

logQ <- read_logq(file.path(input_dir, "logQ.csv"))
size <- read_size(file.path(input_dir, "size.csv"))

if (ncol(logQ) != length(size)) {
  stop(sprintf(
    "column mismatch: logQ has %d states but size has %d entries",
    ncol(logQ),
    length(size)
  ))
}

if (is.na(base_state)) {
  sampled_states <- which(size > 0L)
  if (length(sampled_states) == 0L) {
    stop("size must contain at least one sampled state")
  }
  base_state <- sampled_states[[length(sampled_states)]]
}

if (base_state < 1L || base_state > length(size)) {
  stop("base-state-index is out of range")
}

fit <- UWHAM::uwham(logQ = logQ, size = size, base = base_state, fisher = TRUE)

ze <- as.numeric(fit$ze)
if (base_state > length(ze)) {
  stop(sprintf("base-state-index %d is out of range for %d states", base_state, length(ze)))
}

ze <- ze - ze[[base_state]]
delta_f <- outer(ze, ze, function(a, b) b - a)

write_vector_csv(file.path(input_dir, "ze.csv"), ze, "ze")
write_matrix_csv(file.path(input_dir, "delta_f.csv"), delta_f)

message(sprintf("Wrote ze.csv and delta_f.csv to %s", input_dir))
