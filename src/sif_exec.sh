#!/bin/bash
# Run a command inside envs/squidpy/squidpy.sif.
#   bash src/sif_exec.sh python -c "import scanpy; print('ok')"
#   pixi run ipy        (pixi re-splits task arguments and drops their quotes,
#                        so call this script directly for python -c '...')
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
sif="$root/envs/squidpy/squidpy.sif"
if [[ ! -f "$sif" ]]; then
    echo "missing $sif; build it with: pixi run build-sif" >&2
    exit 1
fi

if ! command -v apptainer >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
    module load apptainer/1.4.5
fi

# HPC3 sets SINGULARITY_BINDPATH (/dfs9,/dfs8,/data/homezvol1,/share);
# pass it on under the current name, and add /dfs9 only when nothing binds it
binds="${APPTAINER_BINDPATH:-${SINGULARITY_BINDPATH:-}}"
unset SINGULARITY_BINDPATH
case ",$binds," in
    *,/dfs9,*) ;;
    *) binds="${binds:+$binds,}/dfs9" ;;
esac
export APPTAINER_BINDPATH="$binds"

exec apptainer exec "$sif" "$@"
