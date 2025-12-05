#!/usr/bin/env python
# coding: utf-8
# # Integrate (batch correct) single cell datasets with scanpy and scvi-tools

# %% [markdown]
# ## load required packages

# %%
import scanpy as sc
import scanpy.external as sce
from pathlib import Path
import anndata as ad
import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import scvi
import grthub_tools as gt


# %% [markdown]
# ## set scanpy default display settings 

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')

# Parse command line arguments

parser = argparse.ArgumentParser(description="Integrate single cell datasets with scanpy and scvi-tools", add_help=False)
parser.add_argument('--help', '-h', action='help', default=argparse.SUPPRESS,
                    help='Show this help message and exit')
parser.add_argument('--filtered_matrix_dirs', nargs='+', help='A list cellranger output directories to process')
# parser.add_argument('--input_dir', type=str, required=True, help='Directory containing sample CellRanger output folders', default='/dfs9/ucightf-lab/projects/GreiS/250930_GreiS/output/cellranger')
parser.add_argument('--output_prefix', type=str, default="output/scanpy", help='Prefix for output files (default: input_dir name)')
parser.add_argument('--min_genes', type=int, default=300, help='Minimum genes per cell for filtering')
parser.add_argument('--min_cells', type=int, default=5, help='Minimum cells per gene for filtering')
parser.add_argument('--n_top_genes', type=int, default=2000, help='Number of highly variable genes to select')
parser.add_argument('--batch_key', type=str, default='batch', help='Batch key for integration')
parser.add_argument('--metadata', type=str, default='metadata.csv', help='csv file with metadata for samples')
args = parser.parse_args()

# directory_path = Path(args.input_dir)

# if args.output_prefix is None:
#     output_prefix = str(directory_path.name)
# else:
#     output_prefix = args.output_prefix

output_prefix = args.output_prefix

# Look for directories containing outs/filtered_feature_bc_matrix/ with 10x files
data_paths = []

for filtered_matrix_dir in args.filtered_matrix_dirs:
    matrix_file = os.path.join(filtered_matrix_dir, 'matrix.mtx.gz')
    features_file = os.path.join(filtered_matrix_dir, 'features.tsv.gz')
    barcodes_file = os.path.join(filtered_matrix_dir, 'barcodes.tsv.gz')
    
    if (os.path.exists(matrix_file) and 
        os.path.exists(features_file) and 
        os.path.exists(barcodes_file)):
        data_paths.append(filtered_matrix_dir)

print(data_paths)

# ## Reading in data
# 
# After reading in the data we'll perform basic filtering on our expression matrix to remove low-quality cells and uninformative genes. The parameter "min_genes" will keep cells that have at least 300 genes, and similarly, "min_cells" will keep genes that are expressed in at least 5 cells.

combined_adata_path = f"{output_prefix}.h5ad"
if os.path.exists(combined_adata_path):
    print(f"{combined_adata_path} already exists. Skipping data loading and concatenation.")
    combined_adata = ad.read_h5ad(combined_adata_path)
else:
    adatas = []
    for i, data_path in enumerate(data_paths):
        # Extract sample ID from the path (e.g., from /path/to/01001CM_011923/outs/filtered_feature_bc_matrix)
        sample_id = os.path.basename(os.path.dirname(os.path.dirname(data_path)))
        print(f"Processing sample {sample_id} from {data_path}")
        
        # Read 10x data
        adata = gt.read_10x_data(data_path, min_genes=args.min_genes, min_cells=args.min_cells)
        
        # Add batch information
        adata.obs[args.batch_key] = sample_id
        
        # Preprocess the data
        adata = gt.preprocess_adata(adata, n_top_genes=args.n_top_genes)
        
        adatas.append(adata)
    
    combined_adata = ad.concat(adatas, join='outer', uns_merge="first")

    metadata_df = pd.read_csv(args.metadata, index_col=0)

    adata.obs = pd.merge(adata.obs, metadata_df, on="batch", how="left", validate="many_to_one")

    combined_adata.write(combined_adata_path)
    del adatas

print(output_prefix)
combined_adata = gt.preprocess_adata(combined_adata, n_top_genes=args.n_top_genes, batch = True)
gt.integrate_w_harmony(combined_adata, output_prefix)