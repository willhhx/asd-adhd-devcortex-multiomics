#!/usr/bin/env python3
"""Step 8a: split the multi-annotation LD scores into single-annotation sets.

Input: <ld-dir>/states.<chr>.l2.ldscore.gz, .l2.M, .l2.M_5_50 from 06b_compute_ldscores.sh
Output: <out>/<annotation>.<chr>.l2.ldscore.gz (+ .M, .M_5_50), one annotation each, so the
regression sees baseline-LD + union + ONE state.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import annotations as an  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ld-prefix", required=True, help="stem before the chromosome number, e.g. data/ldscores/states.")
    ap.add_argument("--names", required=True, help="annot_names.txt from step 6a")
    ap.add_argument("--chroms", default=",".join(str(i) for i in range(1, 23)))
    ap.add_argument("--annot-pattern", help="e.g. work/annot/states.{chrom}.annot.gz; also writes the per-annotation "
                    ".annot.gz files that ldsc --overlap-annot needs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    with open(args.names) as fh:
        names = [ln.strip() for ln in fh if ln.strip()]
    stems = an.split_ldscore_set(args.ld_prefix, names, args.out, args.chroms.split(","))
    if args.annot_pattern:
        an.split_annot_set(args.annot_pattern, names, args.out, args.chroms.split(","))
    print(f"split {len(stems)} annotations into {args.out}")


if __name__ == "__main__":
    main()
