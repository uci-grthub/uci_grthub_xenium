#!/bin/bash
#SBATCH --job-name=rsync_xenium
#SBATCH --account=sbsandme_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/rsync_xenium/%x_%j.log

# Mirror a Xenium delivery (SRC_DIR, e.g. <slide>/<run>/output-XETG*/...) into DEST_DIR,
# keeping its directory tree, then verify every file by checksum against the source
# (Xenium deliveries ship no md5sums.txt).
#
# Usage (from the project root; the log dir must exist before sbatch):
#   mkdir -p logs/rsync_xenium
#   sbatch src/rsync_xenium.sh SRC_DIR [DEST_DIR]
# DEST_DIR defaults to data and is created if missing; it must end up group ucightf
# with setgid set.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch src/rsync_xenium.sh SRC_DIR [DEST_DIR]" >&2
    exit 1
fi

SRC=$(realpath "$1")
DEST="${2:-data}"
mkdir -p "${DEST}"
DEST=$(realpath "${DEST}")

if ! find "${SRC}" -name experiment.xenium -print -quit | grep -q .; then
    echo "ERROR: no experiment.xenium under ${SRC}" >&2
    exit 1
fi

# -rlt, not -a: never copy perms/owner/group. Source dirs can be group
# ucightf_lab_share / ccbs_lab_share (1-byte quota on dfs9), and -a with "SRC/" also
# rewrites DEST's own mode and group, dropping its setgid bit. DEST is setgid ucightf,
# so new files inherit ucightf; --chmod gives new entries group-writable, setgid-dir
# modes instead of the source's 755/744.
if [[ "$(stat -c '%G' "${DEST}")" != ucightf || ! -g "${DEST}" ]]; then
    echo "ERROR: ${DEST} must be group ucightf with setgid set" >&2
    exit 1
fi

echo "Copying ${SRC} -> ${DEST}"
rsync -rltv --partial --chmod=D2775,F664 "${SRC}/" "${DEST}/"

echo "Verifying checksums"
DIFFS=$(rsync -rlt --checksum --dry-run --itemize-changes "${SRC}/" "${DEST}/" \
    | grep -E '^[<>ch*]' || true)
if [[ -n "${DIFFS}" ]]; then
    echo "ERROR: files differ from source after copy:" >&2
    echo "${DIFFS}" >&2
    exit 1
fi
echo "All files match source"
