#!/usr/bin/env Rscript

# compare_outputs.R -- verify a fresh end-to-end pipeline run reproduces the
# committed golden fixtures in example/output/.
#
# Usage:  Rscript scripts/compare_outputs.R <fresh_dir> <golden_dir>
#
# Compares each pipeline output table in <fresh_dir> against the same-named file
# in <golden_dir>: numeric columns must agree to 6 significant figures (the
# tolerance used by tests/testthat/test-IVEA.R), other columns exactly. Exits 0
# if every file matches, 1 (with a per-file report) otherwise.
#
# The comparison is deliberately insensitive to two kinds of *meaningless* ordering,
# both of which drift with the pandas/bioframe version and carry no information:
#   1. Row order within a table (these intermediate files are not all sorted; R's
#      load_data() keys rows by name, never by position).
#   2. Token order inside comma-joined cells -- e.g. the set of gene names
#      overlapping an element, or the promoter's regulatory-element coordinate
#      list. These are unordered sets concatenated with ",".
# To defeat both without hiding real differences, every comma-joined cell is
# canonicalised (split on ",", sort, re-join) and rows are then sorted by their
# character columns (which uniquely identify each row here) before the value
# comparison. A genuine change -- a different gene set, coordinate, count, or score
# -- still surfaces as a mismatch.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: compare_outputs.R <fresh_dir> <golden_dir>", call. = FALSE)
}
fresh_dir  <- args[[1]]
golden_dir <- args[[2]]

# file -> whether it has a header row (TRUE for the labelled tables, FALSE for BEDPE)
manifest <- c(
  "enhancer_elements.txt"            = TRUE,
  "gene_promoter_elements.txt"       = TRUE,
  "enhancer-gene_contacts.chr22.txt" = TRUE,
  "gene_expressions.txt"             = TRUE,
  "predictions_info.chr22.txt"       = TRUE,
  "predictions_score.chr22.bedpe"    = FALSE
)

SIGNIF <- 6L  # significant figures for numeric comparison (matches test-IVEA.R)

read_tbl <- function(path, header) {
  read.table(path, sep = "\t", header = header, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

# Order-normalise comma-joined token lists (unordered sets written with ","):
# a cell with no comma, NA, or a numeric value is returned unchanged.
canon_cell <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || !grepl(",", s, fixed = TRUE)) return(s)
    paste(sort(strsplit(s, ",", fixed = TRUE)[[1]]), collapse = ",")
  }, character(1), USE.NAMES = FALSE)
}

# Canonicalise every character column, then sort rows by those columns so that
# row order (and within-cell token order) does not affect the comparison. The
# character columns -- element coordinates, gene names/ids -- uniquely identify
# each row in these tables, giving a deterministic alignment on both sides.
canonicalise <- function(df) {
  char <- !vapply(df, is.numeric, logical(1))
  for (j in which(char)) df[[j]] <- canon_cell(df[[j]])
  if (any(char)) {
    key <- do.call(paste, c(df[char], sep = "\x1f"))
    df <- df[order(key), , drop = FALSE]
  }
  df
}

# Mismatching cells in one column (NA-safe). Numeric pairs compare at SIGNIF
# figures; anything else compares as strings so int/double/char read differences
# do not masquerade as value differences.
compare_column <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    a <- signif(a, SIGNIF)
    b <- signif(b, SIGNIF)
  } else {
    a <- as.character(a)
    b <- as.character(b)
  }
  sum(!((a == b) | (is.na(a) & is.na(b))), na.rm = TRUE)
}

compare_file <- function(name, header) {
  fresh_path  <- file.path(fresh_dir, name)
  golden_path <- file.path(golden_dir, name)
  if (!file.exists(fresh_path))  return(sprintf("%s: MISSING in fresh output", name))
  if (!file.exists(golden_path)) return(sprintf("%s: MISSING golden fixture", name))

  a <- read_tbl(fresh_path, header)
  b <- read_tbl(golden_path, header)

  if (!identical(dim(a), dim(b))) {
    return(sprintf("%s: dimension mismatch -- fresh %s vs golden %s",
                   name, paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x")))
  }
  if (header && !identical(names(a), names(b))) {
    return(sprintf("%s: column-name mismatch", name))
  }

  a <- canonicalise(a)
  b <- canonicalise(b)

  bad <- vapply(seq_len(ncol(a)), function(j) compare_column(a[[j]], b[[j]]), integer(1))
  if (any(bad > 0L)) {
    cols   <- if (header) names(a) else paste0("V", seq_len(ncol(a)))
    detail <- paste(sprintf("%s=%d", cols[bad > 0L], bad[bad > 0L]), collapse = ", ")
    return(sprintf("%s: %d column(s) differ (mismatched cells: %s)",
                   name, sum(bad > 0L), detail))
  }
  NULL  # match
}

problems <- character(0)
for (name in names(manifest)) {
  msg <- compare_file(name, manifest[[name]])
  if (is.null(msg)) {
    cat(sprintf("OK    %s\n", name))
  } else {
    cat(sprintf("FAIL  %s\n", msg))
    problems <- c(problems, msg)
  }
}

if (length(problems) > 0L) {
  cat(sprintf("\n%d file(s) did not reproduce the golden fixtures.\n", length(problems)))
  quit(status = 1L)
}
cat("\nAll pipeline outputs reproduced the golden fixtures.\n")
