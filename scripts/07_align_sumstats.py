#!/usr/bin/env python3
"""Step 7b: put all traits on the SAME SNP set before S-LDSC.

The paired ASD-vs-ADHD jackknife (devcortex.jackknife.contrast) only works if every trait is
regressed on the same SNPs, so the block partition is identical. This script intersects the
munged ``*.sumstats.gz`` files (columns SNP A1 A2 Z N from munge_sumstats.py), applies one
common chi-square outlier filter, and writes aligned copies.
"""
import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inputs", nargs="+", required=True, metavar="TRAIT=PATH")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-chisq", type=float, default=80.0)
    args = ap.parse_args()

    traits, frames = [], []
    for item in args.inputs:
        trait, path = item.split("=", 1)
        df = pd.read_csv(path, sep=r"\s+")
        need = {"SNP", "A1", "A2", "Z", "N"}
        if not need <= set(df.columns):
            raise SystemExit(f"{path}: expected columns {sorted(need)}, got {list(df.columns)}")
        traits.append(trait)
        frames.append(df.drop_duplicates("SNP").set_index("SNP"))

    common = frames[0].index
    for f in frames[1:]:
        common = common.intersection(f.index)
    ok = pd.Series(True, index=common)
    for f in frames:
        ok &= (f.loc[common, "Z"] ** 2 <= args.max_chisq)
    keep = common[ok.to_numpy()]
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"SNPs: per-trait {[len(f) for f in frames]}, shared {len(common)}, kept after chi2<={args.max_chisq}: {len(keep)}")
    for t, f in zip(traits, frames):
        out = f.loc[keep].reset_index()[["SNP", "A1", "A2", "Z", "N"]]
        out.to_csv(os.path.join(args.out_dir, f"{t}.sumstats.gz"), sep="\t", index=False, compression="gzip")
        print(f"  {t}: mean chi2 = {np.mean(out['Z'] ** 2):.3f}, N median = {out['N'].median():.0f}")


if __name__ == "__main__":
    main()
