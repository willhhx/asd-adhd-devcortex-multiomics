#!/usr/bin/env python3
"""Prepare a PGC CDG3 (Grotzinger et al.) Genomic-SEM factor GWAS for munge_sumstats.py.

Input columns: SNP CHR BP MAF A1 A2 BETA SE P Q_P (BETA is the unstandardised SNP->factor effect, so
Z = BETA/SE is valid for LD-score regression). The files carry no per-SNP N; the README gives an
implied sample size per factor, passed to munge as a constant ``--N``. tau* and enrichment do not
depend on a constant N (it cancels in tau* = tau*sd/(h2/M)), so a mis-stated N does not bias them.

QSNP rule (pre-specified in docs/analysis_plan_v3.md): drop SNPs whose factor-specific QSNP
heterogeneity p-value is < --qsnp-threshold (default 5e-8), i.e. SNPs that act through a single
indicator rather than through the shared factor (README option 2). Use --qsnp-threshold 0 to keep all.
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--qsnp-threshold", type=float, default=5e-8)
a = ap.parse_args()
df = pd.read_csv(a.inp, sep=r"\s+")
need = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Q_P"}
if not need <= set(df.columns):
    raise SystemExit(f"{a.inp}: missing columns {sorted(need - set(df.columns))}")
n0 = len(df)
df = df.dropna(subset=["SNP", "BETA", "SE", "P"])
df = df[(df["SE"] > 0)].drop_duplicates("SNP")
df["CHR"] = df["CHR"].astype(str).str.replace("chr", "", regex=False)
df = df[df["CHR"].isin([str(i) for i in range(1, 23)])]
n1 = len(df)
if a.qsnp_threshold > 0:
    drop = df["Q_P"] < a.qsnp_threshold
    df = df[~drop]
else:
    drop = pd.Series(False, index=df.index)
df[["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]].to_csv(a.out, sep="\t", index=False)
print(f"{a.inp}: {n0} rows -> {n1} valid autosomal; dropped {int(drop.sum())} with Q_P < {a.qsnp_threshold:g}; wrote {len(df)} -> {a.out}")
