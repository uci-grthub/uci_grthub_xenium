#!/usr/bin/env Rscript
# Build an Azimuth reference from a labeled scRNA-seq count matrix exported by
# src/h5ad_to_mtx.py (rule azimuth_reference_build).
#
#   SCTransform (v2) -> PCA -> UMAP (return.model) -> AzimuthReference
#
# Writes into --outdir: ref.Rds and idx.annoy (the folder RunAzimuth takes as
# `reference`), and reference_umap.png colored by each label column.
#
# Usage:
#   Rscript src/build_azimuth_reference.R --mtx-dir <dir> --outdir <dir> \
#       --label-columns cell_type [--dims 50]

suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(ggplot2)
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

read_mtx_dir <- function(dir) {
  counts <- ReadMtx(mtx = file.path(dir, "matrix.mtx.gz"),
                    cells = file.path(dir, "barcodes.tsv.gz"),
                    features = file.path(dir, "features.tsv.gz"),
                    feature.column = 1)
  meta <- read.csv(file.path(dir, "metadata.csv"), row.names = 1, check.names = FALSE)
  CreateSeuratObject(counts = counts, meta.data = meta[colnames(counts), , drop = FALSE])
}

opts <- parse_args(commandArgs(trailingOnly = TRUE),
                   list(`mtx-dir` = NULL, outdir = NULL, `label-columns` = NULL, dims = "50"))
label_columns <- strsplit(opts[["label-columns"]], ",")[[1]]
dims <- seq_len(as.integer(opts$dims))

ref <- read_mtx_dir(opts[["mtx-dir"]])
for (col in label_columns) {
  if (!col %in% colnames(ref[[]])) stop("label column not in metadata: ", col)
  ref[[col]] <- factor(ref[[col, drop = TRUE]])
}
message(sprintf("reference: %d cells x %d genes", ncol(ref), nrow(ref)))

ref <- SCTransform(ref, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
ref <- RunPCA(ref, npcs = max(dims), verbose = FALSE)
ref <- RunUMAP(ref, dims = dims, return.model = TRUE, verbose = FALSE)

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
plots <- lapply(label_columns, function(col)
  DimPlot(ref, group.by = col, label = TRUE, repel = TRUE, raster = FALSE) +
    NoLegend() + ggtitle(sprintf("Reference: %s (%d cells)", col, ncol(ref))))
ggsave(file.path(opts$outdir, "reference_umap.png"),
       patchwork::wrap_plots(plots, ncol = 1), width = 9, height = 8 * length(plots), dpi = 150)

azimuth_ref <- AzimuthReference(ref, refUMAP = "umap", refDR = "pca", refAssay = "SCT",
                                dims = dims, metadata = label_columns,
                                reference.version = "1.0.0", verbose = TRUE)
# SaveAzimuthReference pastes file names onto folder, so it needs the slash
SaveAzimuthReference(azimuth_ref, paste0(normalizePath(opts$outdir), "/"))
message("wrote Azimuth reference to ", opts$outdir)
