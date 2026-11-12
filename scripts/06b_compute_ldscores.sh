#!/usr/bin/env bash
# ldsc --l2 for the multi-column annotation, chr 1-22. Needs ldsc (py2.7 original) as $LDSC_PY and plink files.
# Usage: LDSC_PY="python /path/to/ldsc/ldsc.py" scripts/06b_compute_ldscores.sh ANNOT_DIR OUT_DIR PLINK_PREFIX HM3_SNPLIST [JOBS]
# Note: no xargs -I (BSD/macOS xargs caps the -I command template at 255 bytes, which long paths exceed);
# values go through the environment and the chromosome arrives as $0 of the child shell.
set -euo pipefail
A="$1"; O="$2"; P="$3"; S="$4"; J="${5:-4}"; mkdir -p "$O"
: "${LDSC_PY:?set LDSC_PY}"
export LDSC_PY A O P S
seq 1 22 | xargs -P "$J" -n1 sh -c '$LDSC_PY --l2 --bfile "${P}$0" --ld-wind-cm 1 --annot "$A/states.$0.annot.gz" --print-snps "$S" --out "$O/states.$0"'
