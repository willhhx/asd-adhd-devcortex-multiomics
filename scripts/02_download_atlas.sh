#!/usr/bin/env bash
# Human fetal cortex multiome / scATAC atlases. Accessions MUST be verified at download time (see docs/dataset_validation.md).
#   Trevino et al. 2021 Cell: GEO GSE162170 (multiome + scATAC, human cortex ~PCW 17-24 and later)
#   Others listed in config/datasets.tsv with status=TO_VERIFY (Wang 2025 Nature, Zhu/Wang fetal multiome) - fill accession first.
set -euo pipefail
D="${1:-data/raw/atlas}"; mkdir -p "$D"; cd "$D"
GSE="${2:-GSE162170}"
# GEO supplementary files list (ftp over https):
PREFIX="${GSE:0:${#GSE}-3}nnn"
curl -fL "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/suppl/" -o suppl_index.html
echo "Inspect suppl_index.html, then: curl -fL -O https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/suppl/<file>"
echo "Build an h5ad with obs columns configured in config/config.yaml (cell type, age/stage, donor, sample) before step 03."
