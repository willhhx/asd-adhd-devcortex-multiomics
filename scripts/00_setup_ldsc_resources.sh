#!/usr/bin/env bash
# LDSC reference resources (hg19, 1000G Phase 3 EUR) from the Zenodo mirror of the Price-lab files
# (record 8292725; the old Broad URL now 404s). ~1.1 GB. Resumable: re-run to continue.
set -euo pipefail
D="${1:-data/raw/ldsc}"; mkdir -p "$D"; cd "$D"
Z=https://zenodo.org/records/8292725/files
for f in 1000G_Phase3_plinkfiles.tgz 1000G_Phase3_baselineLD_v2.2_ldscores.tgz 1000G_Phase3_weights_hm3_no_MHC.tgz 1000G_Phase3_frq.tgz; do
  if [ ! -f "$f.done" ]; then
    curl -fL -C - --retry 5 --retry-delay 5 -o "$f" "$Z/$f?download=1"
    tar xzf "$f" && touch "$f.done"
  fi
done
[ -f hm3_no_MHC.list.txt ] || curl -fL -C - --retry 5 -o hm3_no_MHC.list.txt "$Z/hm3_no_MHC.list.txt?download=1"
echo "done:"; ls
