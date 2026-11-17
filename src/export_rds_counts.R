# Export raw counts (RNA assay) and labels of a Seurat .rds as Matrix Market,
# one reference per level of a grouping column (e.g. tissue: Skin, Tumor),
# for src/run_tangram.py.
#
# Per group, at most <max_per_label> cells of each label are kept (random,
# seeded), labels with fewer than <min_per_label> cells are dropped, and genes
# are restricted to those in <genes.tsv> (first column with gene symbols, e.g.
# a Xenium cell_feature_matrix/features.tsv.gz, whose second column is used).
#
# <overrides.csv> (optional) relabels cells before subsampling: columns cell
# and <label_column>, e.g. tumor subtypes from a subclustered object whose
# cells are a subset of this one.
#
# Writes <outdir>/<group>/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz,obs.csv.gz}
#
# Usage: Rscript src/export_rds_counts.R <file.rds> <outdir> <label_column>
#            <group_column> <genes.tsv(.gz)> [max_per_label] [min_per_label]
#            [overrides.csv(.gz)]
suppressPackageStartupMessages({
    library(SeuratObject)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5)
rds <- args[1]
outdir <- args[2]
label_col <- args[3]
group_col <- args[4]
genes_file <- args[5]
max_per_label <- if (length(args) >= 6) as.integer(args[6]) else 2000L
min_per_label <- if (length(args) >= 7) as.integer(args[7]) else 20L
overrides_file <- if (length(args) >= 8) args[8] else NA

panel <- data.table::fread(genes_file, header = FALSE)
# Xenium features.tsv: id, symbol, type; keep gene rows, use the symbol
if (ncol(panel) >= 3) panel <- panel[panel[[3]] == "Gene Expression"]
panel <- if (ncol(panel) >= 2) panel[[2]] else panel[[1]]

obj <- readRDS(rds)
meta <- obj@meta.data
rna <- obj[["RNA"]]
# Seurat v3/v4 Assay keeps counts in a slot; v5 Assay5 in a layer
counts <- if (inherits(rna, "Assay5")) LayerData(rna, layer = "counts") else rna@counts
rm(obj, rna)
invisible(gc())

if (!is.na(overrides_file)) {
    ov <- data.table::fread(overrides_file)
    ov <- ov[ov$cell %in% rownames(meta)]
    cat(nrow(ov), "cells relabeled from", overrides_file, "\n")
    meta[[label_col]] <- as.character(meta[[label_col]])
    meta[ov$cell, label_col] <- ov[[label_col]]
}

genes <- intersect(rownames(counts), panel)
cat(length(genes), "of", length(unique(panel)), "panel genes in the reference\n")
counts <- counts[genes, , drop = FALSE]

set.seed(0)
for (group in sort(unique(as.character(meta[[group_col]])))) {
    cells <- rownames(meta)[meta[[group_col]] == group]
    labels <- as.character(meta[cells, label_col])
    keep <- unlist(lapply(split(cells, labels), function(x) {
        if (length(x) < min_per_label) return(character(0))
        if (length(x) > max_per_label) sample(x, max_per_label) else x
    }), use.names = FALSE)
    keep <- cells[cells %in% keep]

    gdir <- file.path(outdir, group)
    dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
    mtx <- file.path(gdir, "matrix.mtx")
    writeMM(counts[, keep, drop = FALSE], mtx)
    R.utils::gzip(mtx, overwrite = TRUE)
    data.table::fwrite(data.table::data.table(genes), file.path(gdir, "features.tsv.gz"),
                       col.names = FALSE)
    data.table::fwrite(data.table::data.table(keep), file.path(gdir, "barcodes.tsv.gz"),
                       col.names = FALSE)
    obs <- cbind(cell = keep, meta[keep, , drop = FALSE])
    data.table::fwrite(obs, file.path(gdir, "obs.csv.gz"))

    cat("\n", group, ": ", length(keep), " cells\n", sep = "")
    print(sort(table(obs[[label_col]]), decreasing = TRUE))
}
