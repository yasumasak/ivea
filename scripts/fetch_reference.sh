#!/usr/bin/env bash
#
# fetch_reference.sh {hg19|hg38} [gencode_release]
#
# Idempotent downloader for the large reference file that IVEA does NOT vendor in
# the repository: the GENCODE annotation GTF (used to build the RSEM reference and,
# for the Gencode-based workflow, gene BED / expression tables).
#
# Everything else IVEA needs is already committed under reference/ (chrom sizes,
# blacklist, collapsed gene bounds, burst sizes) for both hg19 and hg38 -- see
# reference/README.md for the full build map. Only the GTF is fetched on demand
# because it is ~40-50 MB.
#
# GENCODE release:
#   Defaults to release 26, which is the release used in the paper (hg19 vendors
#   gencode.v26lift37; hg38 uses the native GRCh38 gencode.v26). This keeps the two
#   provided reference sets on the same annotation as CollapsedGeneBounds.*.bed and
#   the committed burst sizes.
#
#   A newer release can be requested as the second argument, e.g.:
#       bash scripts/fetch_reference.sh hg38 48
#   IMPORTANT: if you switch to a release other than 26, the committed
#   CollapsedGeneBounds.*.bed and *_burst_sizes.* were built from v26 and will no
#   longer match the GTF's gene set. Regenerate those from the same release for a
#   self-consistent run (see reference/README.md).
#
# Usage:
#   bash scripts/fetch_reference.sh hg38        # native GRCh38, GENCODE v26 (default)
#   bash scripts/fetch_reference.sh hg38 48      # native GRCh38, GENCODE v48
#   bash scripts/fetch_reference.sh hg19        # GRCh37-mapped, GENCODE v26 (vendored)
#
set -euo pipefail

build="${1:-}"
release="${2:-26}"

if ! [[ "$release" =~ ^[0-9]+$ ]]; then
  echo "ERROR: gencode_release must be a number (got '${release}')." >&2
  exit 2
fi

# Resolve reference/ relative to this script so the command works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ref_dir="${script_dir}/../reference"
gencode_base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${release}"

# Known-good byte sizes for the pinned default release, used as a strict check.
# For any other release we fall back to a gzip-integrity test (size is unknown).
declare -A expected_bytes=( ["hg19-26"]=49915673 ["hg38-26"]=37607127 )

case "$build" in
  hg19)
    url="${gencode_base}/GRCh37_mapping/gencode.v${release}lift37.annotation.gtf.gz"
    out="${ref_dir}/gencode.v${release}lift37.annotation.gtf.gz"
    ;;
  hg38)
    url="${gencode_base}/gencode.v${release}.annotation.gtf.gz"
    out="${ref_dir}/gencode.v${release}.annotation.gtf.gz"
    ;;
  *)
    echo "Usage: bash $0 {hg19|hg38} [gencode_release]" >&2
    echo "Downloads the GENCODE GTF (default release 26) for the given build into reference/." >&2
    exit 2
    ;;
esac

want="${expected_bytes[${build}-${release}]:-}"

file_bytes() { stat -c '%s' "$1" 2>/dev/null || stat -f '%z' "$1" 2>/dev/null || echo 0; }

ok_file() {  # a file is "complete" if it matches the known size (when known) or is a valid gzip
  local f="$1"
  [[ -f "$f" ]] || return 1
  if [[ -n "$want" ]]; then
    [[ "$(file_bytes "$f")" == "$want" ]]
  else
    gzip -t "$f" 2>/dev/null
  fi
}

if ok_file "$out"; then
  echo "[fetch_reference] ${build} (GENCODE v${release}): already present and complete -> $out"
  exit 0
fi

echo "[fetch_reference] ${build} (GENCODE v${release}): downloading"
echo "                  from $url"
echo "                  to   $out"
tmp="${out}.part"
if command -v curl >/dev/null 2>&1; then
  curl -fSL --retry 3 -o "$tmp" "$url"
elif command -v wget >/dev/null 2>&1; then
  wget -O "$tmp" "$url"
else
  echo "[fetch_reference] ERROR: neither curl nor wget is available." >&2
  exit 1
fi

if [[ -n "$want" ]]; then
  got_bytes="$(file_bytes "$tmp")"
  if [[ "$got_bytes" != "$want" ]]; then
    echo "[fetch_reference] ERROR: size mismatch (got ${got_bytes}, expected ${want} bytes)." >&2
    echo "                  Leaving partial file at $tmp for inspection." >&2
    exit 1
  fi
elif ! gzip -t "$tmp" 2>/dev/null; then
  echo "[fetch_reference] ERROR: downloaded file is not a valid gzip." >&2
  echo "                  Leaving partial file at $tmp for inspection." >&2
  exit 1
fi

mv "$tmp" "$out"
echo "[fetch_reference] done: $out ($(file_bytes "$out") bytes)"
if [[ "$release" != "26" ]]; then
  echo "[fetch_reference] NOTE: release ${release} != 26. CollapsedGeneBounds.*.bed and the" >&2
  echo "                  committed burst sizes were built from v26; regenerate them from" >&2
  echo "                  release ${release} for a self-consistent run (see reference/README.md)." >&2
fi
