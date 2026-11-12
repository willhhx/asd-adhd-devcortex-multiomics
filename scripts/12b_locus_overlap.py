#!/usr/bin/env python3
"""Step 12b (pre-specified analysis L1): PIP-mass enrichment of fine-mapped loci in each pooled
cell-state group's peaks. Flat-prior PIPs only (no circularity with S-LDSC). One row per group,
BH-FDR across groups within the trait. Output: <out>/<trait>_locus_overlap.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import finemap, jackknife as jk, locus_overlap, peaks  # noqa: E402
from devcortex.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--trait", required=True)
    ap.add_argument("--gwas", required=True, help="SNP CHR BP A1 A2 BETA SE P (hg19)")
    ap.add_argument("--peakset-dir", required=True, help="hg19 BED directory with states/")
    ap.add_argument("--peakset-summary", required=True)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    fm = load_config(args.config)["finemap"]
    os.makedirs(args.out, exist_ok=True)

    g = pd.read_csv(args.gwas, sep="\t")
    g["CHR"] = g["CHR"].astype(str).str.replace("chr", "", regex=False)
    loci = finemap.define_loci(g, p_thresh=fm["p_threshold"], window=fm["window"])
    print(f"{args.trait}: {len(loci)} loci at p < {fm['p_threshold']:g}")
    if loci.empty:
        sys.exit("no genome-wide-significant loci")
    subs = []
    for L in loci.itertuples():
        s = g[(g["CHR"] == str(L.chrom)) & g["BP"].between(L.start, L.end)].sort_values("BP")
        s = s.drop_duplicates("BP")
        if len(s) < 50:
            continue
        s = s.assign(locus=L.locus, pip=finemap.pip(s["BETA"], s["SE"]))
        subs.append(s)
    allv = pd.concat(subs)
    summ = pd.read_csv(args.peakset_summary, sep="\t")
    rows = []
    for state, fname in zip(summ["state"], summ["file"]):
        bed = peaks.read_bed(os.path.join(args.peakset_dir, "states", fname + ".bed"))
        hit = peaks.variant_peak_overlap(allv, bed).notna()
        pips = [s["pip"].to_numpy() for s in subs]
        inpk, off = [], 0
        for s in subs:
            inpk.append(hit.iloc[off:off + len(s)].to_numpy()); off += len(s)
        r = locus_overlap.pip_mass_enrichment(pips, inpk, n_perm=args.n_perm, seed=2026)
        rows.append({"trait": args.trait, "state": state, **r})
    out = pd.DataFrame(rows)
    out["fdr"] = jk.benjamini_hochberg(out["p_one_sided"])
    path = os.path.join(args.out, f"{args.trait}_locus_overlap.tsv")
    out.sort_values("p_one_sided").to_csv(path, sep="\t", index=False)
    print(out.sort_values("p_one_sided").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
