#!/usr/bin/env python3
"""Step 4b (secondary, post hoc): class x stage-group peak sets from the step-3 pseudobulk.

Writes the same layout as step 4 (<out>/states/*.bed, union_peaks.bed, peakset_summary.tsv) so
steps 5, 6, 6b, 8 and 9 run unchanged. union_peaks.bed is copied from the primary analysis so
the 'union' covariate is identical. Also writes pooled_membership.tsv: which states went into
which pooled group, and which were left out ('Other', e.g. Cajal-Retzius, Vascular, Unknown).
"""
import argparse
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import peaks, pooling, specificity  # noqa: E402
from devcortex.config import ensure_dir, load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--pseudobulk-dir", required=True)
    ap.add_argument("--primary-peakset-dir", required=True, help="step 4 output (hg38), for union_peaks.bed")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ps = load_config(args.config)["peaksets"]

    counts = pd.read_csv(os.path.join(args.pseudobulk_dir, "atac_counts.tsv.gz"), sep="\t", index_col=0)
    feat = pd.read_csv(os.path.join(args.pseudobulk_dir, "atac_features.tsv"), sep="\t", index_col=0).loc[counts.index]
    pooled, member = pooling.pool_counts(counts)
    ensure_dir(args.out)
    member.to_csv(os.path.join(args.out, "pooled_membership.tsv"), sep="\t", index=False)
    left = member[member["pooled"] == ""]
    print(f"{len(member)} states -> {pooled.shape[1]} pooled groups; {len(left)} states left out:")
    print("  " + ", ".join(sorted(set(left["cell_type"]))))
    print(member[member["pooled"] != ""].groupby("pooled").size().to_string())
    union = pd.read_csv(os.path.join(args.primary_peakset_dir, "union_peaks.bed"), sep="\t", header=None)[3]
    keep = pooled.index.intersection(pd.Index(union))
    top = specificity.top_specific(pooled.loc[keep], top_n=ps["top_n"], method=ps["method"], min_cpm=ps["min_cpm"])
    sd = ensure_dir(os.path.join(args.out, "states"))
    rows = []
    for g, ids in top.items():
        name = specificity.safe_name(g)
        bed = feat.loc[ids]
        peaks.write_bed(bed, os.path.join(sd, f"{name}.bed"))
        rows.append({"state": g, "file": name, "n_peaks": len(ids),
                     "total_bp": int((bed["end"] - bed["start"]).sum()),
                     "mean_width": float((bed["end"] - bed["start"]).mean())})
    shutil.copy(os.path.join(args.primary_peakset_dir, "union_peaks.bed"), os.path.join(args.out, "union_peaks.bed"))
    pd.DataFrame(rows).to_csv(os.path.join(args.out, "peakset_summary.tsv"), sep="\t", index=False)
    print(f"wrote {len(rows)} pooled peak sets of up to {ps['top_n']} peaks -> {args.out}")


if __name__ == "__main__":
    main()
