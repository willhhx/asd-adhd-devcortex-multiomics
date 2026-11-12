#!/usr/bin/env bash
# GWAS summary statistics. Both are public but URLs/format change: CHECK the landing pages, then edit the two lines below.
#  ASD  : Grove et al. 2019 Nat Genet (iPSYCH-PGC), PGC download page  https://pgc.unc.edu/for-researchers/download-results/
#  ADHD : Demontis et al. 2023 Nat Genet (ADHD2022 iPSYCH-deCODE-PGC), same PGC page
#  SCZ (positive control, PGC3 2022) and a height or similar negative control (GIANT/Yengo 2022) are on the same site / GWAS Catalog.
set -euo pipefail
D="${1:-data/raw/gwas}"; mkdir -p "$D"
echo "Download the files from the PGC page into $D as: ASD.tsv.gz ADHD.tsv.gz SCZ.tsv.gz HEIGHT.tsv.gz"
echo "Then run: python scripts/01b_standardise_gwas.py --help  (maps column names, writes SNP CHR BP A1 A2 BETA SE P + munge input)"
