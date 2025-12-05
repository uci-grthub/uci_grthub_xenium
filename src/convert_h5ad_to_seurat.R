#!/usr/bin/env Rscript
library(Seurat)
library(reticulate)
library(sceasy)

# Simple command-line argument parsing to accept --counts, --scaled_data and --output
args <- commandArgs(trailingOnly = TRUE)

print_usage <- function() {
	cat("convert_h5ad_to_seurat.R - convert an AnnData (.h5ad) file to a Seurat .rds file\n")
	cat("\nUsage:\n")
	cat("  Rscript src/convert_h5ad_to_seurat.R --data <data.h5ad> [--scaled_data <scaled.h5ad>] [--output <file.rds>]\n\n")
	cat("Options:\n")
	cat("  --data  Path to input data .h5ad file (default: output/scanpy/combined_harmony_integrated.h5ad)\n")
	cat("  --scaled_data  Path to scaled matrix (.mtx) with companion barcode/gene CSVs (optional)\n")
	cat("  --output  Path to output .rds file (default: same as counts but with .rds extension)\n")
	invisible(NULL)
}

# defaults
data_file <- NULL
scaled_data_file <- NULL
output_file <- NULL

if (length(args) == 0) {
	# No args provided: use existing defaults
	data_file <- "output/scanpy/combined_harmony_integrated.h5ad"
	output_file <- "output/scanpy/combined_harmony_integrated.rds"
} else {
	i <- 1
	while (i <= length(args)) {
		a <- args[i]
		if (a == "-h" || a == "--help") {
			print_usage()
			quit(status = 0)
		} else if (startsWith(a, "--data=")) {
			data_file <- sub("^--data=", "", a)
		} else if (a == "--data") {
			if ((i+1) <= length(args)) {
				data_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--data requires a value")
			}
		} else if (startsWith(a, "--scaled_data=")) {
			scaled_data_file <- sub("^--scaled_data=", "", a)
		} else if (a == "--scaled_data") {
			if ((i+1) <= length(args)) {
				scaled_data_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--scaled_data requires a value")
			}
		} else if (startsWith(a, "--output=")) {
			output_file <- sub("^--output=", "", a)
		} else if (a == "--output") {
			if ((i+1) <= length(args)) {
				output_file <- args[i+1]
				i <- i + 1
			} else {
				stop("--output requires a value")
			}
		} else {
			stop(paste("Unknown argument:", a))
		}
		i <- i + 1
	}

	# set defaults if missing
	if (is.null(data_file)) {
		data_file <- "output/scanpy/bp_harmony_integrated.h5ad"
	}
	if (is.null(output_file)) {
		if (grepl("\\.h5ad$", data_file, ignore.case = TRUE)) {
			output_file <- sub("\\.h5ad$", ".rds", data_file, ignore.case = TRUE)
		} else {
			output_file <- paste0(data_file, ".rds")
		}
	}
}

# Ensure conda env is used for reticulate (if needed in the environment)
try({
	use_condaenv("scvi-tools", conda = "/opt/apps/mamba/24.3.0/bin/mamba")
}, silent = TRUE)


message("Converting: ", data_file, " -> ", output_file)

# perform conversion
message("Converting counts AnnData: ", data_file, " -> ", output_file)
sceasy::convertFormat(data_file, from = "anndata", to = "seurat", outFile = output_file, main_layer = "data")

if (!is.null(scaled_data_file) && nzchar(scaled_data_file) && file.exists(scaled_data_file)) {
	message("Converting scaled AnnData (if provided): ", scaled_data_file)
	sceasy::convertFormat(scaled_data_file, from = "anndata", to = "seurat", outFile = file.path(dirname(output_file), paste0(tools::file_path_sans_ext(basename(output_file)), "_scaled_data.rds")), main_layer = "scale.data")
} else {
	message("No scaled_data provided or file not found; skipping scaled data conversion.")
}

seu <- readRDS(output_file)
Idents(seu) <- seu$batch

# seu <- NormalizeData(seu)
# seu <- FindVariableFeatures(seu, selection.method = "vst")
# seu <- ScaleData(seu, features = VariableFeatures(seu), block.size = 200)

saveRDS(seu, file = output_file)

message("Done. Saved Seurat object to: ", output_file)
