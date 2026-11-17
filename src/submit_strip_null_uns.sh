#!/bin/bash
#SBATCH --job-name=strip_null_uns
#SBATCH -A sbsandme_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/strip_null_uns.%j.out
#SBATCH --error=logs/strip_null_uns.%j.err

# Make the per-sample, merged/ and merged_nocorrection/ AnnData readable by
# anndata 0.11: src/strip_null_uns.py deletes the "null"-encoded uns entries
# anndata 0.12 writes (mtimes kept), then every file is read back in full with
# the anndata 0.11.4 env (samui).
# Usage, from the repository root:
#   sbatch src/submit_strip_null_uns.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

PY=/dfs9/ucightf-lab/kstachel/.conda/envs/squidpy_dfs/bin/python
PY011=/dfs9/ucightf-lab/kstachel/.conda/envs/samui/bin/python
H5ADS=(output/xenium/anndata/*.h5ad output/xenium/anndata/merged/*.h5ad
       output/xenium/anndata/merged_nocorrection/*.h5ad)

rc=0
PYTHONNOUSERSITE=1 "$PY" src/strip_null_uns.py "${H5ADS[@]}" || rc=$?

PYTHONNOUSERSITE=1 "$PY011" - "${H5ADS[@]}" <<'PYEOF'
import sys
import anndata as ad
from anndata._io.specs.registry import IORegistryError
failed = 0
for path in sys.argv[1:]:
    try:
        a = ad.read_h5ad(path)
    except (IORegistryError, OSError) as e:
        failed += 1
        print(f"anndata {ad.__version__} FAILED {path}: {type(e).__name__}: {e}", flush=True)
        continue
    cols = [c for c in ("tangram_cell_type", "seurat_celltype") if c in a.obs]
    print(f"anndata {ad.__version__} OK {path}: {a.n_obs} x {a.n_vars}, {cols}", flush=True)
sys.exit(1 if failed else 0)
PYEOF
exit "$rc"
