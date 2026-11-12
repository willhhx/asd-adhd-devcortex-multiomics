"""Peak identifiers, BED reading/writing and variant-to-peak overlap."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

_PEAK_RE = re.compile(r"^(?P<chrom>(?:chr)?[0-9XYM]+)[:\-_](?P<start>\d+)[:\-_](?P<end>\d+)$")


def parse_peak_ids(ids) -> pd.DataFrame:
    """Parse 'chr1:1000-2000', 'chr1-1000-2000' or 'chr1_1000_2000' into chrom/start/end."""
    rows = []
    for i in ids:
        m = _PEAK_RE.match(str(i))
        if not m:
            raise ValueError(f"cannot parse peak id {i!r} (expected chrN:start-end)")
        rows.append((m["chrom"], int(m["start"]), int(m["end"])))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"], index=pd.Index(list(ids), name="peak"))
    return df


def norm_chrom(c) -> str:
    return str(c).replace("chr", "")


def write_bed(df: pd.DataFrame, path: str) -> None:
    """BED4: chrom, start, end, peak_id (peak id = index)."""
    out = df[["chrom", "start", "end"]].copy()
    out["peak"] = df.index
    out.to_csv(path, sep="\t", header=False, index=False)


def read_bed(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    names = ["chrom", "start", "end", "peak"][: df.shape[1]]
    df = df.iloc[:, : len(names)]
    df.columns = names
    if "peak" not in df.columns:
        df["peak"] = [f"{c}:{s}-{e}" for c, s, e in zip(df.chrom, df.start, df.end)]
    return df


def variant_peak_overlap(variants: pd.DataFrame, peaks: pd.DataFrame,
                         chrom_col="CHR", pos_col="BP") -> pd.Series:
    """For each variant, the id of the peak containing it (or NaN).

    Peaks are assumed non-overlapping within a BED (true for a common peak set). A 1-based
    variant position p lies in BED interval [start, end) iff start < p <= end.
    """
    out = pd.Series(np.nan, index=variants.index, dtype=object)
    pk = peaks.assign(_c=peaks["chrom"].map(norm_chrom))
    for chrom, sub in pk.groupby("_c"):
        sub = sub.sort_values("start")
        starts, ends, ids = sub["start"].to_numpy(), sub["end"].to_numpy(), sub["peak"].to_numpy()
        mask = variants[chrom_col].map(norm_chrom) == chrom
        pos = variants.loc[mask, pos_col].to_numpy()
        idx = np.searchsorted(starts, pos, side="left") - 1
        ok = idx >= 0
        hit = np.zeros(len(pos), bool)
        hit[ok] = pos[ok] <= ends[idx[ok]]
        res = np.full(len(pos), None, dtype=object)
        res[hit] = ids[idx[hit]]
        out.loc[mask] = res
    return out
