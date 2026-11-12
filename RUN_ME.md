# Commands for you to run (macOS terminal). Nothing here has run on real data yet.

Everything below was written and tested on synthetic data in the cloud sandbox, which has no internet. Real downloads, `ldsc`, `liftover` chain download and `h5ad` reading are the untested parts.

## 0. Environment (once)
```bash
cd ~/<Python Intro>/asd-adhd-devcortex-multiomics
conda env create -f envs/environment.yml && conda activate devcortex
export PYTHONPATH=$PWD/src
pytest -q tests                                   # expect: 28 passed
python synthetic/run_e2e.py --work /tmp/e2e && python synthetic/check_truth.py --work /tmp/e2e   # expect: 0 failure(s)
# ldsc: no conda needed (conda was timing out)
bash scripts/00b_setup_ldsc.sh ~/ldsc      # then paste the two export lines it prints
```
If conda works for you, `envs/environment.yml` is optional: `pip install numpy scipy pandas matplotlib pyyaml pytest pyliftover` is enough.
## 1. Downloads
```bash
bash scripts/00_setup_ldsc_resources.sh data/raw/ldsc        # ~4 GB, hg19 LD reference
bash scripts/01_download_gwas.sh data/raw/gwas                # prints what to fetch from PGC (manual, check terms)
bash scripts/02_download_atlas.sh data/raw/atlas GSE162170    # lists GEO supplementary files; pick files by hand
```
Then: build `data/raw/atlas/atlas.h5ad` (obs: cell_type, stage, donor, sample; var = peaks `chrN:start-end`), edit `config/config.yaml`,
standardise GWAS with `scripts/01b_standardise_gwas.py` (see `--help`) and munge:
```bash
for t in ASD ADHD SCZ HEIGHT; do $MUNGE_PY --sumstats data/gwas/$t.tsv.gz --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --N-cas <Ncases> --N-con <Ncontrols> --out data/munged/$t; done   # quantitative traits: use --N
```
## 2. Pipeline
```bash
python scripts/03_build_pseudobulk.py --config config/config.yaml --format h5ad --input data/raw/atlas/atlas.h5ad --out work/pseudobulk
python scripts/04_define_peak_sets.py --config config/config.yaml --pseudobulk-dir work/pseudobulk --out work/peaksets_atlas_build
python scripts/05_liftover_beds.py --in-dir work/peaksets_atlas_build --out-dir work/peaksets_hg19      # skip build change only if atlas is already hg19: add --assume-same-build
python scripts/06_make_annotations.py --peakset-dir work/peaksets_hg19 --bim-prefix data/raw/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC. --out work/annot
bash scripts/06b_compute_ldscores.sh work/annot work/ldscores data/raw/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC. data/raw/ldsc/hm3_no_MHC.list.txt 4
python scripts/08_split_ldscores.py --ld-prefix work/ldscores/states. --names work/annot/annot_names.txt --out work/ldscores_single
python scripts/07_align_sumstats.py --in ASD=data/munged/ASD.sumstats.gz ADHD=data/munged/ADHD.sumstats.gz SCZ=data/munged/SCZ.sumstats.gz HEIGHT=data/munged/HEIGHT.sumstats.gz --out-dir work/sumstats_aligned
bash scripts/08_run_sldsc.sh work/sumstats_aligned work/ldscores_single data/raw/ldsc/baselineLD. data/raw/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. data/raw/ldsc/1000G_Phase3_frq/1000G.EUR.QC. work/sldsc 4
python scripts/09_summarise_sldsc.py --config config/config.yaml --sldsc-dir work/sldsc --peakset-summary work/peaksets_hg19/peakset_summary.tsv --out work/results
python scripts/12_finemap_and_link.py --config config/config.yaml --trait ADHD --gwas data/gwas/ADHD.tsv.gz --tau-table work/results/tau_star_by_state.tsv --peakset-dir work/peaksets_hg19 --peakset-summary work/peaksets_hg19/peakset_summary.tsv --out work/finemap
python scripts/14_make_figures.py --results work/results --finemap-dir work/finemap --stage-order first_trimester,second_trimester,third_trimester,infancy --out work/figs
```
Check paths of the LD files after untar (names differ slightly); `06b`/`08` print the exact ldsc commands. Expect hours for 06b+08 (≈ states × traits regressions; tune JOBS).
Not written yet / flagged: Fig 1 atlas UMAP, MAGMA/scDRS, SuSiE/PolyFun R, real coloc eQTL harmonisation.
