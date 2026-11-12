#!/usr/bin/env bash
# Whole pooled class x stage-group analysis (steps 4b-9) for one peak-set size.
# Usage (project root):  LDSC_PY="conda run -n ldsc python $HOME/ldsc/ldsc.py" bash scripts/run_pooled_variant.sh TOP_N OUT_DIR [JOBS]
# e.g.  bash scripts/run_pooled_variant.sh 50000 work/pooled_top50k 6
set -euo pipefail
N="$1"; W="$2"; J="${3:-6}"; : "${LDSC_PY:?set LDSC_PY}"
export LDSC_PY
CFG="$W/config.yaml"; mkdir -p "$W"
sed "s/top_n: [0-9]*/top_n: $N/" config/config_pooled.yaml > "$CFG"
PLINK=data/raw/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.
python scripts/04b_pool_classes.py --config "$CFG" --pseudobulk-dir work/pseudobulk --primary-peakset-dir work/peaksets_atlas_build --out "$W/peaksets_hg38"
python scripts/05_liftover_beds.py --in-dir "$W/peaksets_hg38" --out-dir "$W/peaksets_hg19" --chain data/raw/liftover/hg38ToHg19.over.chain.gz
python scripts/06_make_annotations.py --peakset-dir "$W/peaksets_hg19" --bim-prefix "$PLINK" --out "$W/annot"
bash scripts/06b_compute_ldscores.sh "$W/annot" "$W/ldscores" "$PLINK" data/raw/ldsc/hm3_no_MHC.list.txt 4
[ "$(ls "$W"/ldscores/*.l2.ldscore.gz | wc -l)" -eq 22 ] || { echo "expected 22 LD-score files"; exit 1; }
python scripts/08_split_ldscores.py --ld-prefix "$W/ldscores/states." --names "$W/annot/annot_names.txt" --out "$W/ldscores_single" --annot-pattern "$W/annot/states.{chrom}.annot.gz"
bash scripts/08_run_sldsc.sh work/sumstats_aligned "$W/ldscores_single" data/raw/ldsc/baselineLD. data/raw/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. data/raw/ldsc/1000G_Phase3_frq/1000G.EUR.QC. "$W/sldsc" "$J"
python scripts/09_summarise_sldsc.py --config "$CFG" --sldsc-dir "$W/sldsc" --peakset-summary "$W/peaksets_hg19/peakset_summary.tsv" --out "$W/results"
echo "DONE: $W/results"
