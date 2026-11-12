#!/usr/bin/env python3
"""Step 13: colocalisation of a GWAS locus with an eQTL (or any second trait) via coloc.abf.

Inputs are the harmonised GWAS (SNP CHR BP A1 A2 BETA SE P) and a table with SNP, BETA, SE for
the second trait in the SAME effect-allele orientation (harmonise alleles upstream; the
fetal-brain eQTL availability and allele coding must be checked per resource). The region is
the window around ``--snp`` that exists in both tables.
Writes one line of posterior probabilities (PP.H0..H4) and the top shared-variant candidates.
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import coloc  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--qtl", required=True)
    ap.add_argument("--snp", required=True, help="lead/anchor SNP of the locus")
    ap.add_argument("--window", type=int, default=500_000)
    ap.add_argument("--label", default="locus")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    g = pd.read_csv(args.gwas, sep="\t").drop_duplicates("SNP").set_index("SNP")
    q = pd.read_csv(args.qtl, sep="\t").drop_duplicates("SNP").set_index("SNP")
    if args.snp not in g.index:
        sys.exit(f"anchor SNP {args.snp} not in GWAS")
    chrom, bp = g.loc[args.snp, "CHR"], g.loc[args.snp, "BP"]
    near = g[(g["CHR"] == chrom) & ((g["BP"] - bp).abs() <= args.window)]
    shared = near.index.intersection(q.index)
    if len(shared) < 50:
        sys.exit(f"only {len(shared)} variants shared in the window; coloc needs a dense region")
    r = coloc.coloc_abf(near.loc[shared, "BETA"], near.loc[shared, "SE"], q.loc[shared, "BETA"], q.loc[shared, "SE"])
    snp_pp = pd.Series(r.pop("snp_pp_h4"), index=shared).sort_values(ascending=False)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    row = {"label": args.label, "anchor_snp": args.snp, **{k: v for k, v in r.items()},
           "top_shared_snp": snp_pp.index[0], "top_shared_snp_pp": float(snp_pp.iloc[0])}
    pd.DataFrame([row]).to_csv(args.out, sep="\t", index=False)
    print(pd.DataFrame([row]).T.to_string(header=False))


if __name__ == "__main__":
    main()
