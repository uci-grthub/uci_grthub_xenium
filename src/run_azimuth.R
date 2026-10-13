#!/usr/bin/env Rscript
# Map one Xenium sample onto an Azimuth reference (rule xenium_azimuth).
#
# Input is the sample's raw counts exported by src/h5ad_to_mtx.py; the
# reference folder (ref.Rds, idx.annoy) comes from
# src/build_azimuth_reference.R. Writes one CSV row per cell: cell_id,
# predicted.<label> and predicted.<label>.score for each reference label,
# mapping.score, and the projected reference UMAP (refUMAP_1, refUMAP_2).
#
# Usage:
#   Rscript src/run_azimuth.R --mtx-dir <dir> --reference <dir> --out <csv>

suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
})

parse_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(defaults)) stop("unknown argument: ", args[i])
    opts[[key]] <- args[i + 1]
    i <- i + 2
  }
  missing <- names(opts)[vapply(opts, is.null, logical(1))]
  if (length(missing)) stop("missing arguments: --", paste(missing, collapse = ", --"))
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE),
                   list(`mtx-dir` = NULL, reference = NULL, out = NULL))
dir <- opts[["mtx-dir"]]

counts <- ReadMtx(mtx = file.path(dir, "matrix.mtx.gz"),
                  cells = file.path(dir, "barcodes.tsv.gz"),
                  features = file.path(dir, "features.tsv.gz"),
                  feature.column = 1)
query <- CreateSeuratObject(counts = counts)
message(sprintf("query: %d cells x %d genes", ncol(query), nrow(query)))

query <- RunAzimuth(query, reference = opts$reference, assay = "RNA")

meta <- query[[]]
keep <- grep("^predicted\\.|^mapping\\.score$", colnames(meta), value = TRUE)
out <- data.frame(cell_id = rownames(meta), meta[, keep, drop = FALSE], check.names = FALSE)
umap <- Embeddings(query, "ref.umap")
out$refUMAP_1 <- umap[out$cell_id, 1]
out$refUMAP_2 <- umap[out$cell_id, 2]

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
write.csv(out, opts$out, row.names = FALSE)
message(sprintf("wrote %d cells to %s", nrow(out), opts$out))
