"""Cell-state-specific peak (or gene) sets from a pseudobulk matrix.

Input is a features x states matrix of raw pseudobulk counts. Output is, for every state,
the ``top_n`` features most specific to that state, so every state's annotation has the
same number of peaks (this removes annotation size as a hidden driver of tau*).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def to_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return counts.divide(lib, axis=1) * 1e6


def specificity_scores(counts: pd.DataFrame, method: str = "log2fc",
                       pseudocount: float = 1.0) -> pd.DataFrame:
    """Specificity of each feature for each state.

    log2fc : log2((CPM_s + c) / (mean CPM of all OTHER states + c))
    zscore : (log CPM_s - mean over states) / sd over states
    """
    cpm = to_cpm(counts)
    n = cpm.shape[1]
    if n < 3:
        raise ValueError("need at least 3 states to define specificity")
    if method == "log2fc":
        total = cpm.sum(axis=1)
        out = {}
        for s in cpm.columns:
            others = (total - cpm[s]) / (n - 1)
            out[s] = np.log2((cpm[s] + pseudocount) / (others + pseudocount))
        return pd.DataFrame(out, index=cpm.index)
    if method == "zscore":
        lg = np.log2(cpm + pseudocount)
        mu = lg.mean(axis=1)
        sd = lg.std(axis=1, ddof=1).replace(0, np.nan)
        return lg.sub(mu, axis=0).div(sd, axis=0)
    raise ValueError(f"unknown method {method!r}")


def union_features(counts: pd.DataFrame, min_cpm: float = 1.0, min_states: int = 1) -> pd.Index:
    """Features detected (CPM >= min_cpm) in at least ``min_states`` states."""
    cpm = to_cpm(counts)
    keep = (cpm >= min_cpm).sum(axis=1) >= min_states
    return counts.index[keep.to_numpy()]


def top_specific(counts: pd.DataFrame, top_n: int, method: str = "log2fc",
                 min_cpm: float = 1.0) -> dict[str, pd.Index]:
    """Top ``top_n`` specific features per state (must also reach ``min_cpm`` in that state)."""
    cpm = to_cpm(counts)
    spec = specificity_scores(counts, method=method)
    out = {}
    for s in counts.columns:
        ok = cpm[s] >= min_cpm
        sc = spec.loc[ok, s].sort_values(ascending=False, kind="mergesort")
        out[s] = sc.index[:top_n]
    return out


def overlap_matrix(sets: dict[str, pd.Index]) -> pd.DataFrame:
    """Jaccard overlap between state-specific sets (reported as a supplementary diagnostic;
    heavily overlapping states cannot be separated by one-annotation-at-a-time S-LDSC)."""
    names = list(sets)
    mat = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        sa = set(sets[a])
        for j, b in enumerate(names):
            sb = set(sets[b])
            u = len(sa | sb)
            mat[i, j] = len(sa & sb) / u if u else np.nan
    return pd.DataFrame(mat, index=names, columns=names)


def safe_name(state: str) -> str:
    """File-system-safe state label: 'EN-L2_3-IT|second_trimester' -> 'EN-L2_3-IT__second_trimester'."""
    keep = []
    for ch in state:
        keep.append(ch if (ch.isalnum() or ch in "-_.") else "__" if ch == "|" else "_")
    return "".join(keep)


def peaks_to_bed(peak_ids: pd.Index, peak_table: pd.DataFrame) -> pd.DataFrame:
    """peak_table indexed by peak id with columns chrom, start, end -> BED3 for ``peak_ids``."""
    sub = peak_table.loc[peak_ids, ["chrom", "start", "end"]].copy()
    return sub.reset_index(drop=True)
