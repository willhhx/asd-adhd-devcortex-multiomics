#!/usr/bin/env bash
# One S-LDSC regression per trait x state: baselineLD + union + <state>. Needs 08_split_ldscores.py output.
# Usage: LDSC_PY="..." scripts/08_run_sldsc.sh SUMSTATS_DIR SINGLE_DIR BASELINE_PREFIX WEIGHTS_PREFIX FRQ_PREFIX OUT_DIR [JOBS]
set -euo pipefail
SS="$1"; SD="$2"; BL="$3"; W="$4"; FRQ="$5"; O="$6"; J="${7:-4}"; : "${LDSC_PY:?set LDSC_PY}"
mkdir -p "$O"
for f in "$SS"/*.sumstats.gz; do t=$(basename "$f" .sumstats.gz); mkdir -p "$O/$t"
  for s in $(basename -a "$SD"/*.1.l2.ldscore.gz | sed 's/\.1\.l2\.ldscore\.gz$//' | grep -v '^union$'); do
    echo "$t $s"; done; done > "$O/jobs.txt"
xargs -P "$J" -L1 sh -c '[ -f "'"$O"'/$0/$1.results" ] || '"$LDSC_PY"' --h2 "'"$SS"'/$0.sumstats.gz" --ref-ld-chr "'"$BL"'",'"$SD"'/union.,'"$SD"'/$1. --w-ld-chr "'"$W"'" --overlap-annot --frqfile-chr "'"$FRQ"'" --print-coefficients --print-delete-vals --out "'"$O"'/$0/$1"' < "$O/jobs.txt"
