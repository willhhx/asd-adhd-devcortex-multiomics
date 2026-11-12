#!/usr/bin/env python3
"""Step 6a: BED peak sets -> one multi-column LDSC annotation file per chromosome.

Column 1 is the union-of-all-peaks annotation, the rest are the cell states. Writes
<out>/states.<chr>.annot.gz aligned with <bim-prefix><chr>.bim, plus <out>/annot_names.txt
(column order). Step 6b (06b_compute_ldscores.sh) turns these into LD scores with ldsc.py.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import annotations as an  # noqa: E402
from devcortex import peaks  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peakset-dir", required=True, help="hg19 BED directory (states/ + union_peaks.bed)")
    ap.add_argument("--bim-prefix", required=True, help="e.g. data/raw/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.")
    ap.add_argument("--chroms", default=",".join(str(i) for i in range(1, 23)))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    states_dir = os.path.join(args.peakset_dir, "states")
    beds = {"union": peaks.read_bed(os.path.join(args.peakset_dir, "union_peaks.bed"))}
    for fn in sorted(os.listdir(states_dir)):
        if fn.endswith(".bed"):
            beds[fn[:-4]] = peaks.read_bed(os.path.join(states_dir, fn))
    names = list(beds)
    with open(os.path.join(args.out, "annot_names.txt"), "w") as fh:
        fh.write("\n".join(names) + "\n")

    total = {n: 0 for n in names}
    for chrom in args.chroms.split(","):
        bim = an.read_bim(f"{args.bim_prefix}{chrom}.bim")
        annot = an.make_annot_matrix(bim, beds)
        an.write_annot(annot, bim, os.path.join(args.out, f"states.{chrom}.annot.gz"))
        for n in names:
            total[n] += int(annot[n].sum())
    print("SNPs per annotation (over processed chromosomes):")
    for n in names:
        print(f"  {n}: {total[n]}")
    small = [n for n, v in total.items() if v < 2000 and n != "union"]
    if small:
        print(f"WARNING: {len(small)} annotations have <2,000 SNPs and will be underpowered: {small[:5]} ...")


if __name__ == "__main__":
    main()
