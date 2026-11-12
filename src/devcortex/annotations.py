"""Build LDSC annotation files from peak BED files and split multi-annotation LD scores.

LDSC needs, per chromosome, an ``.annot.gz`` aligned row-for-row with the reference
PLINK ``.bim``. We compute LD scores for ALL state annotations at once (one multi-column
annotation file per chromosome, which shares the expensive LD computation), then split
the resulting ``.l2.ldscore.gz`` / ``.l2.M`` / ``.l2.M_5_50`` into single-annotation sets so
each regression sees: baseline-LD + union-of-peaks + ONE cell state.
"""
from __future__ import annotations

import gzip
import os
from typing import Iterable

import numpy as np
import pandas as pd

BIM_COLS = ["CHR", "SNP", "CM", "BP", "A1", "A2"]


def read_bim(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None, names=BIM_COLS,
                       dtype={"CHR": str}, engine="python")


def snps_in_intervals(bp: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Boolean mask of SNP positions inside any half-open interval [start, end).

    BED is 0-based half-open and bim positions are 1-based, so a SNP at 1-based position p
    lies in a BED interval iff start < p <= end.
    """
    order = np.argsort(starts)
    s, e = starts[order], ends[order]
    # merge overlapping intervals so a single searchsorted is exact
    merged_s, merged_e = [], []
    for a, b in zip(s, e):
        if merged_s and a <= merged_e[-1]:
            merged_e[-1] = max(merged_e[-1], b)
        else:
            merged_s.append(a)
            merged_e.append(b)
    ms, me = np.asarray(merged_s), np.asarray(merged_e)
    idx = np.searchsorted(ms, bp, side="left") - 1  # last interval with start < bp
    ok = idx >= 0
    inside = np.zeros(len(bp), dtype=bool)
    inside[ok] = bp[ok] <= me[idx[ok]]
    return inside


def make_annot_matrix(bim: pd.DataFrame, beds: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """0/1 annotation matrix (rows = bim SNPs, columns = annotation names)."""
    out = {}
    for name, bed in beds.items():
        sub = bed[bed["chrom"].astype(str).str.replace("chr", "", regex=False)
                  == str(bim["CHR"].iloc[0]).replace("chr", "")]
        if sub.empty:
            out[name] = np.zeros(len(bim), dtype=np.int8)
        else:
            out[name] = snps_in_intervals(bim["BP"].to_numpy(), sub["start"].to_numpy(),
                                          sub["end"].to_numpy()).astype(np.int8)
    return pd.DataFrame(out)


def write_annot(annot: pd.DataFrame, bim: pd.DataFrame, path: str) -> None:
    """Write an LDSC annot file: CHR BP SNP CM then annotation columns."""
    df = pd.concat([bim[["CHR", "BP", "SNP", "CM"]].reset_index(drop=True),
                    annot.reset_index(drop=True)], axis=1)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def read_annot(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def write_ldscore(df_ld: pd.DataFrame, bim: pd.DataFrame, path_prefix: str,
                  m: np.ndarray, m_5_50: np.ndarray) -> None:
    """Write ``{prefix}.l2.ldscore.gz`` (+ .l2.M, .l2.M_5_50) in LDSC layout.

    Columns of the ldscore file are CHR SNP BP then one ``<name>L2`` per annotation.
    """
    out = pd.concat([bim[["CHR", "SNP", "BP"]].reset_index(drop=True),
                     df_ld.reset_index(drop=True)], axis=1)
    out.columns = list(out.columns[:3]) + [f"{c}L2" for c in df_ld.columns]
    out.to_csv(path_prefix + ".l2.ldscore.gz", sep="\t", index=False,
               float_format="%.4f", compression="gzip")
    for suffix, vec in ((".l2.M", m), (".l2.M_5_50", m_5_50)):
        with open(path_prefix + suffix, "w") as fh:
            fh.write("\t".join(str(int(x)) for x in vec) + "\n")


def split_ldscore_set(multi_prefix: str, names: Iterable[str], out_dir: str,
                      chroms: Iterable[str]) -> dict[str, str]:
    """Split multi-column LD-score sets (per chromosome) into one set per annotation.

    ``multi_prefix`` is the path stem before the chromosome number (e.g. ``ld/states.``,
    giving ``ld/states.1.l2.ldscore.gz`` ...). Returns {name: single-annotation stem}, where
    the stem is again "stem before chromosome", ready for ``--ref-ld-chr``.
    """
    names = list(names)
    os.makedirs(out_dir, exist_ok=True)
    stems = {n: os.path.join(out_dir, f"{n}.") for n in names}
    for chrom in chroms:
        ld = pd.read_csv(f"{multi_prefix}{chrom}.l2.ldscore.gz", sep="\t")
        meta = ld[["CHR", "SNP", "BP"]]
        cols = [c for c in ld.columns if c.endswith("L2")]
        if len(cols) != len(names):
            raise ValueError(f"chr{chrom}: {len(cols)} ldscore columns vs {len(names)} names")
        with open(f"{multi_prefix}{chrom}.l2.M") as fh:
            m = np.array(fh.read().split(), dtype=np.int64)
        with open(f"{multi_prefix}{chrom}.l2.M_5_50") as fh:
            m550 = np.array(fh.read().split(), dtype=np.int64)
        for i, n in enumerate(names):
            single = pd.concat([meta, ld[[cols[i]]].rename(columns={cols[i]: f"{n}L2"})], axis=1)
            single.to_csv(f"{stems[n]}{chrom}.l2.ldscore.gz", sep="\t", index=False,
                          float_format="%.4f", compression="gzip")
            for suffix, vec in ((".l2.M", m), (".l2.M_5_50", m550)):
                with open(f"{stems[n]}{chrom}{suffix}", "w") as fh:
                    fh.write(str(int(vec[i])) + "\n")
    return stems


def split_annot_set(annot_pattern: str, names: Iterable[str], out_dir: str, chroms: Iterable[str]) -> None:
    """Write one single-annotation ``.annot.gz`` per name and chromosome (``<out>/<name>.<chrom>.annot.gz``).

    ldsc --overlap-annot reads the .annot.gz belonging to every --ref-ld-chr prefix to build the
    Prop._SNPs / Enrichment columns, so each single-annotation LD-score set needs its own annot file.
    ``annot_pattern`` contains ``{chrom}``, e.g. ``work/annot/states.{chrom}.annot.gz`` (CHR BP SNP CM + columns).
    """
    names = list(names)
    os.makedirs(out_dir, exist_ok=True)
    for chrom in chroms:
        df = pd.read_csv(annot_pattern.format(chrom=chrom), sep="\t")
        meta = df[["CHR", "BP", "SNP", "CM"]]
        for n in names:
            if n not in df.columns:
                raise ValueError(f"annotation {n!r} not in {annot_pattern.format(chrom=chrom)}")
            pd.concat([meta, df[[n]]], axis=1).to_csv(
                os.path.join(out_dir, f"{n}.{chrom}.annot.gz"), sep="\t", index=False, compression="gzip")
