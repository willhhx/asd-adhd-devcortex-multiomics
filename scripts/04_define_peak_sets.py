#!/usr/bin/env python3
"""Step 4: pseudobulk accessibility -> equal-size cell-state-specific peak sets (BED).

Writes <out>/states/<state>.bed (BED4, atlas genome build), <out>/union_peaks.bed (every
peak detected in any state), <out>/peakset_summary.tsv and <out>/peakset_jaccard.tsv.
``--method`` and ``--top-n`` allow the sensitivity analyses (zscore specificity, other sizes).
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import peaks, specificity  # noqa: E402
from devcortex.config import ensure_dir, load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--pseudobulk-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--method", default=None, choices=["log2fc", "zscore"])
    ap.add_argument("--top-n", type=int, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ps = cfg["peaksets"]
    method = args.method or ps["method"]
    top_n = args.top_n or ps["top_n"]

    counts = pd.read_csv(os.path.join(args.pseudobulk_dir, "atac_counts.tsv.gz"), sep="\t", index_col=0)
    feat = pd.read_csv(os.path.join(args.pseudobulk_dir, "atac_features.tsv"), sep="\t", index_col=0)
    feat = feat.loc[counts.index]

    union = specificity.union_features(counts, min_cpm=ps["min_cpm"])
    counts_u = counts.loc[union]
    top = specificity.top_specific(counts_u, top_n=top_n, method=method, min_cpm=ps["min_cpm"])

    states_dir = ensure_dir(os.path.join(args.out, "states"))
    rows = []
    for state, ids in top.items():
        name = specificity.safe_name(state)
        bed = feat.loc[ids]
        peaks.write_bed(bed, os.path.join(states_dir, f"{name}.bed"))
        rows.append({"state": state, "file": name, "n_peaks": len(ids),
                     "total_bp": int((bed["end"] - bed["start"]).sum()),
                     "mean_width": float((bed["end"] - bed["start"]).mean())})
    peaks.write_bed(feat.loc[union], os.path.join(args.out, "union_peaks.bed"))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(args.out, "peakset_summary.tsv"), sep="\t", index=False)
    specificity.overlap_matrix(top).to_csv(os.path.join(args.out, "peakset_jaccard.tsv"), sep="\t")
    # size-equality check: unequal total_bp is a known confound, so report it
    cv = summ["total_bp"].std() / summ["total_bp"].mean()
    print(f"{len(summ)} state peak sets of {top_n} peaks; union {len(union)} peaks; "
          f"total-bp CV across states = {cv:.2f}")


if __name__ == "__main__":
    main()
