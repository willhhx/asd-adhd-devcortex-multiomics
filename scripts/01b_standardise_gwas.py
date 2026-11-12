#!/usr/bin/env python3
"""Map a PGC/GWAS-catalog table to the pipeline format: SNP CHR BP A1 A2 BETA SE P (hg19).

Example: python scripts/01b_standardise_gwas.py --in raw.tsv.gz --out ASD.tsv.gz \
  --map SNP=SNP CHR=CHR BP=BP A1=A1 A2=A2 BETA=OR SE=SE P=P --or-to-beta
Use --map NEW=OLD (new pipeline name = raw column). If beta is given as OR use --or-to-beta (log OR; SE must be SE of log OR).
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--in", dest="inp", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--map", nargs="+", required=True, metavar="NEW=OLD")
ap.add_argument("--or-to-beta", action="store_true")
ap.add_argument("--sep", default=None, help="default: whitespace/auto")
a = ap.parse_args()
df = pd.read_csv(a.inp, sep=a.sep, engine="python") if a.sep else pd.read_csv(a.inp, sep=r"\s+", comment="#")
m = dict(x.split("=") for x in a.map)
out = pd.DataFrame({k: df[v] for k, v in m.items()})
out["CHR"] = out["CHR"].astype(str).str.replace("chr", "", regex=False)
if a.or_to_beta:
    out["BETA"] = np.log(out["BETA"].astype(float))
out["A1"], out["A2"] = out["A1"].str.upper(), out["A2"].str.upper()
out = out.dropna().drop_duplicates("SNP")
out = out[out["CHR"].isin([str(i) for i in range(1, 23)])]
out.to_csv(a.out, sep="\t", index=False)
print(f"{len(out)} variants -> {a.out}")
