#!/usr/bin/env python3
"""
Tangram cell type annotation of a merged Xenium AnnData (from
merge_xenium_anndata.py) with a labeled scRNA-seq reference (from
export_rds_counts.R), and UMAPs colored by the transferred cell type.

Both sides are library-size normalized and log1p transformed from raw counts
(reference matrix.mtx.gz, Xenium layers["counts"]). Training genes are the top
--n-markers marker genes per reference cell type (Wilcoxon) shared with the
panel. Tangram maps the per-type mean profiles onto the Xenium cells
(mode="clusters", uniform density: one segmented cell per spot), then
project_cell_annotations gives each cell a score per type; scores are scaled
to sum to 1 per cell, the top type is the cell's label and its share the
label's score.

Reference types whose markers are not on the panel (--drop-labels, e.g. RBC:
no hemoglobin genes on the Xenium 5K panel) are removed first; otherwise
they attract cells with little signal. Cells with fewer than --min-counts
transcripts are not mapped and are labeled "low_counts".

Writes --predictions, a CSV with cell_id, tangram_<label>,
tangram_<label>_score, total_counts and the per-type probabilities
(prob_<type>), and the plots of src/plot_label_transfer.py into --outdir
(umap_tangram.png, ...). To plot the same predictions on another merged
object of the same cells (e.g. merged_nocorrection/), run
src/plot_label_transfer.py --method tangram.

Usage:
  python src/run_tangram.py --h5ad merged/<group>.h5ad --reference <ref dir> \
      --label cell_type_015 --predictions <csv> --outdir merged/plots/<group>
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import tangram as tg

from plot_label_transfer import LOW_COUNTS, plot_all


def read_reference(refdir, label):
    """Reference counts (cells x genes) with obs from export_rds_counts.R."""
    counts = scipy.io.mmread(os.path.join(refdir, "matrix.mtx.gz")).T.tocsr()
    genes = pd.read_csv(os.path.join(refdir, "features.tsv.gz"), header=None)[0]
    obs = pd.read_csv(os.path.join(refdir, "obs.csv.gz"), index_col="cell", low_memory=False)
    ref = ad.AnnData(X=counts.astype(np.float32), obs=obs)
    ref.var_names = genes.astype(str).to_numpy()
    ref.obs[label] = ref.obs[label].astype("category")
    return ref


def read_spatial(path):
    """Merged Xenium object with X = raw counts (from layers["counts"])."""
    adata = ad.read_h5ad(path)
    counts = adata.layers["counts"]
    return ad.AnnData(X=sp.csr_matrix(counts, dtype=np.float32), obs=adata.obs,
                      var=adata.var[[]], obsm={"X_umap": adata.obsm["X_umap"]})


def lognorm(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def marker_genes(ref, label, n_markers):
    sc.tl.rank_genes_groups(ref, groupby=label, method="wilcoxon", use_raw=False)
    names = pd.DataFrame(ref.uns["rank_genes_groups"]["names"]).head(n_markers)
    return sorted(set(names.to_numpy().ravel()))


def run_tangram(ref, sp_adata, label, genes, epochs, device, seed):
    tg.pp_adatas(ref, sp_adata, genes=genes)
    print(f"{len(sp_adata.uns['training_genes'])} training genes")
    ad_map = tg.map_cells_to_space(ref, sp_adata, mode="clusters", cluster_label=label,
                                   density_prior="uniform", num_epochs=epochs,
                                   device=device, random_state=seed)
    tg.project_cell_annotations(ad_map, sp_adata, annotation=label)
    scores = sp_adata.obsm["tangram_ct_pred"]
    return scores.div(scores.sum(axis=1), axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", required=True, help="merged .h5ad")
    parser.add_argument("--reference", required=True,
                        help="directory from export_rds_counts.R")
    parser.add_argument("--label", default="cell_type_015", help="reference obs label column")
    parser.add_argument("--drop-labels", nargs="*", default=[],
                        help="reference cell types to remove before mapping")
    parser.add_argument("--min-counts", type=int, default=50,
                        help="cells with fewer transcripts are not mapped (low_counts)")
    parser.add_argument("--n-markers", type=int, default=100,
                        help="marker genes per reference cell type used for training")
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--predictions", required=True, help="output CSV")
    parser.add_argument("--outdir", required=True, help="plot directory")
    args = parser.parse_args()

    name = os.path.basename(args.h5ad).removesuffix(".h5ad")
    spatial, pred = annotate(args, name, f"tangram_{args.label}")
    plot_all(spatial, pred, "tangram", args.label, 0.0, name, args.outdir)


def annotate(args, name, key):
    """Run Tangram; returns the Xenium object and the predictions (written to CSV)."""
    ref = read_reference(args.reference, args.label)
    spatial = read_spatial(args.h5ad)
    ref = ref[~ref.obs[args.label].isin(args.drop_labels)].copy()
    ref.obs[args.label] = ref.obs[args.label].cat.remove_unused_categories()
    print(f"reference {ref.shape} (dropped: {', '.join(args.drop_labels) or 'none'}), "
          f"{name} {spatial.shape}")
    print(ref.obs[args.label].value_counts().to_string())

    n_counts = np.asarray(spatial.X.sum(axis=1)).ravel()
    keep = n_counts >= args.min_counts
    print(f"{keep.sum()} of {spatial.n_obs} cells ({keep.mean():.0%}) have "
          f">= {args.min_counts} transcripts and are mapped")
    query = spatial[keep].copy()
    lognorm(ref)
    lognorm(query)

    genes = marker_genes(ref, args.label, args.n_markers)
    probs = run_tangram(ref, query, args.label, genes, args.epochs, args.device, args.seed)

    pred = pd.DataFrame({key: probs.idxmax(axis=1), f"{key}_score": probs.max(axis=1)})
    pred = pd.concat([pred, probs.add_prefix("prob_")], axis=1)
    pred = pred.reindex(spatial.obs_names)
    pred[key] = pred[key].fillna(LOW_COUNTS)
    pred.insert(2, "total_counts", n_counts)
    pred.index.name = "cell_id"
    os.makedirs(os.path.dirname(args.predictions) or ".", exist_ok=True)
    pred.to_csv(args.predictions)
    return spatial, pred


if __name__ == "__main__":
    main()
