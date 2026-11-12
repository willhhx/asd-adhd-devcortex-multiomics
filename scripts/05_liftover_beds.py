#!/usr/bin/env python3
"""Step 5: lift atlas peak BEDs (hg38) to hg19 so they match the LDSC 1000G Phase 3 panel.

Uses the ``pyliftover`` package with the UCSC chain file (hg38ToHg19.over.chain.gz).
Peaks whose two ends do not both lift to the same chromosome, or whose lifted width changes
by more than ``--max-width-change`` (default 20%), are dropped and counted. Use
``--assume-same-build`` (copy only) if the atlas is already hg19, or for the synthetic test.
"""
import argparse
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import peaks  # noqa: E402


def lift_bed(df, lo):
    keep = []
    for chrom, start, end, pid in df[["chrom", "start", "end", "peak"]].itertuples(index=False):
        a = lo.convert_coordinate(chrom, int(start))
        b = lo.convert_coordinate(chrom, int(end))
        if not a or not b or a[0][0] != b[0][0]:
            continue
        s, e = sorted((a[0][1], b[0][1]))
        w0, w1 = end - start, e - s
        if w0 > 0 and abs(w1 - w0) / w0 > 0.2:
            continue
        keep.append((a[0][0], s, e, pid))
    import pandas as pd
    return pd.DataFrame(keep, columns=["chrom", "start", "end", "peak"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--chain", help="path to hg38ToHg19.over.chain.gz")
    ap.add_argument("--assume-same-build", action="store_true")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out_dir, "states"), exist_ok=True)
    beds = glob.glob(os.path.join(args.in_dir, "states", "*.bed")) + [os.path.join(args.in_dir, "union_peaks.bed")]

    if args.assume_same_build:
        for b in beds:
            dest = os.path.join(args.out_dir, os.path.relpath(b, args.in_dir))
            shutil.copy(b, dest)
        print(f"copied {len(beds)} BED files unchanged")
    else:
        try:
            from pyliftover import LiftOver
        except ImportError:
            sys.exit("pip install pyliftover (or run UCSC liftOver yourself and use --assume-same-build)")
        if not args.chain:
            sys.exit("--chain is required (download hg38ToHg19.over.chain.gz from UCSC goldenPath)")
        lo = LiftOver(args.chain)
        for b in beds:
            df = peaks.read_bed(b)
            out = lift_bed(df, lo)
            dest = os.path.join(args.out_dir, os.path.relpath(b, args.in_dir))
            out.to_csv(dest, sep="\t", header=False, index=False)
            print(f"{os.path.basename(b)}: {len(out)}/{len(df)} peaks lifted")
    for fn in ("peakset_summary.tsv", "peakset_jaccard.tsv"):
        src = os.path.join(args.in_dir, fn)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(args.out_dir, fn))


if __name__ == "__main__":
    main()
