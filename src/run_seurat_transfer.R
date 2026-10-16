#!/usr/bin/env Rscript
# Seurat label transfer (FindTransferAnchors + TransferData) of reference
# cell types onto a merged Xenium object.
#
# Query: raw counts of the merged object exported by src/h5ad_to_mtx.py
# (--obs-columns sample_id). Reference: a directory from
# src/export_rds_counts.R (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz,
# obs.csv.gz). The reference is restricted to genes present in the query, its
# variable genes are the anchor features, and its PCA is the space the query
# is projected into. Anchors are found per sample (sample_id), so each sample
# is transferred on its own against the same reference.
#
# --normalization SCT (default): the reference is SCTransformed (v2,
# glmGamPoi) and the query residuals are computed from the reference's SCT
# model (recompute.residuals), as Azimuth does; columns are prefixed
# seurat_sct_. LogNormalize: both are LogNormalized; columns are prefixed
# seurat_.
#
# Reference types whose markers are not on the panel (--drop-labels, comma
# separated, e.g. RBC) are removed first. Cells with fewer than --min-counts
# transcripts are not transferred and are labeled "low_counts".
#
# Writes one CSV row per cell: cell_id, seurat_<label> (top type),
# seurat_<label>_score (its prediction score), total_counts, and
# prob_<type> (prediction score per type), the layout of src/run_tangram.py.
#
# Usage:
#   Rscript src/run_seurat_transfer.R --query-mtx <dir> --reference <dir> \
#       --out <csv> [--label cell_type_015] [--drop-labels RBC] \
#       [--min-counts 50] [--dims 30] [--nfeatures 2000] \
#       [--normalization SCT|LogNormalize]

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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
                   list(`query-mtx` = NULL, reference = NULL, out = NULL,
                        label = "cell_type_015", `drop-labels` = "",
                        `min-counts` = "50", dims = "30", nfeatures = "2000",
                        normalization = "SCT"))
normalization <- match.arg(opts$normalization, c("SCT", "LogNormalize"))
sct <- normalization == "SCT"
label <- opts$label
drop_labels <- strsplit(opts[["drop-labels"]], ",")[[1]]
min_counts <- as.numeric(opts[["min-counts"]])
dims <- seq_len(as.integer(opts$dims))
key <- paste0(if (sct) "seurat_sct_" else "seurat_", label)

read_mtx_dir <- function(dir) {
  ReadMtx(mtx = file.path(dir, "matrix.mtx.gz"),
          cells = file.path(dir, "barcodes.tsv.gz"),
          features = file.path(dir, "features.tsv.gz"),
          feature.column = 1)
}

# query
qdir <- opts[["query-mtx"]]
qcounts <- read_mtx_dir(qdir)
qmeta <- read.csv(file.path(qdir, "metadata.csv"), row.names = 1, check.names = FALSE)
total_counts <- Matrix::colSums(qcounts)
keep <- total_counts >= min_counts
message(sprintf("query: %d cells x %d genes; %d (%.0f%%) with >= %g transcripts are transferred",
                ncol(qcounts), nrow(qcounts), sum(keep), 100 * mean(keep), min_counts))

# reference, restricted to genes in the query
rdir <- opts$reference
rcounts <- read_mtx_dir(rdir)
robs <- data.table::fread(file.path(rdir, "obs.csv.gz"), data.table = FALSE)
rownames(robs) <- robs$cell
robs <- robs[colnames(rcounts), , drop = FALSE]
rkeep <- !robs[[label]] %in% drop_labels
genes <- intersect(rownames(rcounts), rownames(qcounts))
ref <- CreateSeuratObject(counts = rcounts[genes, rkeep],
                          meta.data = robs[rkeep, label, drop = FALSE])
message(sprintf("reference: %d cells x %d shared genes (dropped: %s)", ncol(ref), nrow(ref),
                if (length(drop_labels)) paste(drop_labels, collapse = ", ") else "none"))
print(sort(table(ref[[label]][, 1]), decreasing = TRUE))
types <- sort(unique(ref[[label]][, 1]))

nfeatures <- min(as.integer(opts$nfeatures), nrow(ref))
if (sct) {
  ref <- SCTransform(ref, vst.flavor = "v2", method = "glmGamPoi",
                     variable.features.n = nfeatures, verbose = FALSE)
} else {
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, nfeatures = nfeatures, verbose = FALSE)
  ref <- ScaleData(ref, verbose = FALSE)
}
ref <- RunPCA(ref, npcs = max(dims), verbose = FALSE)
features <- VariableFeatures(ref)
message(sprintf("normalization: %s, %d anchor features", normalization, length(features)))

transfer_sample <- function(cells) {
  query <- CreateSeuratObject(counts = qcounts[genes, cells, drop = FALSE])
  if (sct) {
    # query residuals from the reference's SCT model
    anchors <- FindTransferAnchors(reference = ref, query = query, features = features,
                                   normalization.method = "SCT", reference.assay = "SCT",
                                   query.assay = "RNA", recompute.residuals = TRUE,
                                   reference.reduction = "pca", dims = dims,
                                   k.filter = NA, verbose = FALSE)
  } else {
    query <- NormalizeData(query, verbose = FALSE)
    anchors <- FindTransferAnchors(reference = ref, query = query, features = features,
                                   normalization.method = "LogNormalize",
                                   reference.reduction = "pca", dims = dims,
                                   k.filter = NA, verbose = FALSE)
  }
  n_anchors <- nrow(anchors@anchors)
  # TransferData needs more anchors than k.weight; small samples have few
  pred <- TransferData(anchorset = anchors, refdata = ref[[label]][, 1], dims = dims,
                       k.weight = min(50, n_anchors - 1), verbose = FALSE)
  message(sprintf("  %d cells, %d anchors", length(cells), n_anchors))
  scores <- pred[, paste0("prediction.score.", make.names(types)), drop = FALSE]
  colnames(scores) <- paste0("prob_", types)
  out <- data.frame(pred$predicted.id, pred$prediction.score.max, scores,
                    check.names = FALSE, row.names = rownames(pred))
  colnames(out)[1:2] <- c(key, paste0(key, "_score"))
  out
}

samples <- qmeta[colnames(qcounts), "sample_id"]
parts <- list()
for (s in sort(unique(samples[keep]))) {
  message(s)
  parts[[s]] <- transfer_sample(colnames(qcounts)[keep & samples == s])
}
pred <- do.call(rbind, unname(parts))

# every query cell, low-count ones labeled low_counts with no scores. Indexed
# column by column: indexing the data.frame by row name would make unique row
# names for hundreds of thousands of NA rows, which is quadratic.
idx <- match(colnames(qcounts), rownames(pred))
out <- as.data.frame(lapply(pred, function(col) col[idx]), check.names = FALSE)
out[[key]][!keep] <- "low_counts"
rownames(out) <- colnames(qcounts)
out <- data.frame(cell_id = rownames(out), out[, 1:2], total_counts = total_counts,
                  out[, -(1:2), drop = FALSE], check.names = FALSE)

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, opts$out)
message(sprintf("wrote %d cells to %s", nrow(out), opts$out))
print(sort(table(out[[key]]), decreasing = TRUE))
