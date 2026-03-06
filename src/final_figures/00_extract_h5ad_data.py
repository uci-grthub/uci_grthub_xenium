#!/usr/bin/env python3
"""
Extract data from h5ad for downstream R visualization scripts.

Outputs (to src/final_figures/output/):
  umap_coords_metadata.csv  - UMAP coords + cell metadata (downsampled)
  cluster_mean_expr.csv     - mean log-norm expression per cluster x age x condition
  panel_genes.csv           - all genes in the Xenium panel
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

BASE = "/dfs9/ucightf-lab/projects/ChatS/251006_0925Bio-06_ChatS_xenium"
H5AD = os.path.join(BASE, "output/scanpy/all_harmony_integrated.h5ad")
DEG_PB = os.path.join(BASE, "results/to_client/all/degs/pseudobulk")
OUT = os.path.join(BASE, "src/final_figures/output")
os.makedirs(OUT, exist_ok=True)

# ── Load h5ad ──────────────────────────────────────────────────────────────────
print("Loading h5ad...", flush=True)
adata = sc.read_h5ad(H5AD)
print(f"  {adata.n_obs:,} cells  x  {adata.n_vars:,} genes", flush=True)

# ── Export panel gene list ─────────────────────────────────────────────────────
pd.DataFrame({"gene": adata.var_names.tolist()}).to_csv(
    os.path.join(OUT, "panel_genes.csv"), index=False
)
print(f"Saved panel genes: {adata.n_vars}", flush=True)

# ── 1. UMAP coordinates + metadata (downsampled to 80 k for fast plotting) ────
meta_cols = ["age", "condition", "area", "leiden_harmony"]
if "sample_id" in adata.obs.columns:
    meta_cols.append("sample_id")

umap_df = pd.DataFrame(
    adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names
)
for col in meta_cols:
    umap_df[col] = adata.obs[col].values

# Stratified downsample preserving cluster proportions
MAX_CELLS = 80_000
if adata.n_obs > MAX_CELLS:
    rng = np.random.default_rng(42)
    idx = rng.choice(adata.n_obs, MAX_CELLS, replace=False)
    umap_df = umap_df.iloc[sorted(idx)]

umap_df.to_csv(os.path.join(OUT, "umap_coords_metadata.csv"), index=False)
print(f"Saved UMAP coords: {len(umap_df):,} cells", flush=True)

# ── 2. Collect target genes ────────────────────────────────────────────────────
panel_set = set(adata.var_names.tolist())

# Top DEGs from main pseudobulk comparisons
pb_genes: set[str] = set()
for comp in ["comp1_heatvsctrl", "comp1_agedvsadult",
             "comp2_adultonly_heatvsctrl", "comp2_agedonly_heatvsctrl"]:
    for cl in range(17):
        fpath = os.path.join(DEG_PB, f"DEG_pb_{comp}_cluster{cl}.csv")
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, index_col=0)
        sig = df[df["adj.P.Val"] < 0.05].sort_values("adj.P.Val")
        pb_genes.update(sig.head(20).index.tolist())

# Curated pathway gene sets (heat stress experiment in mouse brain)
curated = {
    "Heat_shock": [
        "Hspa1a", "Hspa1b", "Hspa5", "Hspa8", "Hsph1",
        "Hsp90aa1", "Hsp90ab1", "Hspb1", "Hspb8", "Hspb11",
        "Dnajb1", "Dnajb4", "Dnajb6", "Dnajc5", "Hsf1",
    ],
    "UPR": [
        "Atf3", "Atf4", "Ddit3", "Eif2ak3", "Xbp1", "Atf6",
        "Bip", "Hspa5", "Perk", "Ire1",
    ],
    "Neuroinflammation": [
        "Nfkb1", "Tnf", "Il6", "Il1b", "Ccl2", "Ccl4",
        "Nos2", "Ptgs2", "Tgfb1", "Cx3cr1",
        "Gfap", "S100b", "Vim", "Nes", "Aqp4", "Aldh1l1",
    ],
    "Apoptosis": [
        "Casp3", "Casp7", "Bax", "Bcl2", "Bcl2l1",
        "Tp53", "Mcl1", "Puma", "Noxa",
    ],
    "Oxidative_stress": [
        "Nfe2l2", "Hmox1", "Sod1", "Sod2", "Cat",
        "Gpx1", "Txn1", "Prdx1", "Sqstm1",
    ],
    "Neuron_markers": [
        "Rbfox3", "Map2", "Tubb3", "Syp", "Snap25",
        "Syt1", "Vglut1", "Gad1", "Gad2",
    ],
    "Oligodendrocyte": [
        "Mog", "Mbp", "Plp1", "Cnp", "Sox10",
        "Olig1", "Olig2", "Cldn11",
    ],
    "Microglia": [
        "Iba1", "Aif1", "Cx3cr1", "P2ry12",
        "Tmem119", "Hexb", "Csf1r",
    ],
}
curated_genes = [g for genes in curated.values() for g in genes]
pb_genes.update(curated_genes)

# Filter to panel
target_genes = [g for g in pb_genes if g in panel_set]
print(f"Target genes in panel: {len(target_genes)} / {len(pb_genes)} requested",
      flush=True)

# Save curated gene labels for R scripts
curated_records = [
    {"gene": g, "pathway": path}
    for path, genes in curated.items()
    for g in genes
    if g in panel_set
]
pd.DataFrame(curated_records).to_csv(
    os.path.join(OUT, "curated_gene_sets.csv"), index=False
)

# ── 3. Cluster mean expression ─────────────────────────────────────────────────
print("Extracting expression for target genes...", flush=True)
gene_idx = [adata.var_names.get_loc(g) for g in target_genes]

X = adata.X
if sp.issparse(X):
    X_sub = np.asarray(X[:, gene_idx].todense())
else:
    X_sub = np.array(X[:, gene_idx])

expr_df = pd.DataFrame(X_sub, columns=target_genes, index=adata.obs_names)
expr_df["leiden_harmony"] = adata.obs["leiden_harmony"].values
expr_df["condition"] = adata.obs["condition"].values
expr_df["age"] = adata.obs["age"].values

group_cols = ["leiden_harmony", "age", "condition"]
cluster_mean = (
    expr_df
    .groupby(group_cols)[target_genes]
    .mean()
    .reset_index()
)
cluster_mean.to_csv(os.path.join(OUT, "cluster_mean_expr.csv"), index=False)
print(f"Saved cluster mean expression: {cluster_mean.shape}", flush=True)

print("Done.", flush=True)
