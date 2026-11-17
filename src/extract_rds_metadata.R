# Dump the cell metadata of a Seurat .rds (e.g. cell type labels) without
# exporting its expression data.
#
# Writes into <outdir>, named after the .rds:
#   <name>.meta.csv.gz  meta.data with the cell barcode as column "cell"
#   <name>.summary.txt  dimensions, assays, reductions, active idents, and the
#                       levels and counts of every non-numeric meta.data column
#
# Usage: Rscript src/extract_rds_metadata.R <file.rds> <outdir>
suppressPackageStartupMessages(library(SeuratObject))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
rds <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
name <- sub("\\.rds$", "", basename(rds))

obj <- readRDS(rds)
meta <- obj@meta.data
out <- file(file.path(outdir, paste0(name, ".summary.txt")), "w")
w <- function(...) writeLines(paste0(...), out)
w("file: ", rds)
w("class: ", class(obj)[1], "  version: ", as.character(obj@version))
w("cells: ", ncol(obj), "  features (default assay): ", nrow(obj))
w("assays: ", paste(names(obj@assays), collapse = ", "),
  "  default: ", DefaultAssay(obj))
w("reductions: ", paste(names(obj@reductions), collapse = ", "))
w("active idents: ", paste(utils::head(levels(Idents(obj)), 60), collapse = ", "))
w("meta.data columns: ", paste(colnames(meta), collapse = ", "))
for (col in colnames(meta)) {
    x <- meta[[col]]
    if (is.numeric(x) && length(unique(x)) > 60) next
    tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
    w("\n== ", col, " (", length(tab), " levels)")
    w(paste0("  ", names(tab), ": ", as.integer(tab), collapse = "\n"))
}
close(out)

meta <- cbind(cell = rownames(meta), meta)
data.table::fwrite(meta, file.path(outdir, paste0(name, ".meta.csv.gz")))
cat(name, ": ", nrow(meta), " cells, ", ncol(meta) - 1, " meta.data columns\n", sep = "")
