"""Build one comparable nuclei x peak-tile matrix from per-sample 10x multiome files.

Each sample in the Wang et al. 2025 deposit has its OWN called peak set (ids rarely match
across samples), so raw peak columns cannot be added across samples. We re-express every peak
on a shared grid: the genome is cut into fixed-width tiles (default 500 bp) and each peak is
assigned to the tile containing its midpoint; peak counts landing in one tile are summed.
The consensus feature set is the tiles hit by peaks in >= ``min_samples`` samples.
Limitation (stated in the Methods): peak-level counts, not fragment-level, so a wide peak is
placed in a single tile and sub-tile resolution is lost.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

_PEAK = re.compile(r"^(?:chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
AUTOSOMES = [str(i) for i in range(1, 23)]
_CHROM_BASE = 10**9


def parse_peaks(ids) -> pd.DataFrame:
    """'chr1:10-20' | '1:10-20' | 'KI270713.1:5-9' -> chrom (no 'chr'), start, end; unparsable -> NaN."""
    chrom, start, end = [], [], []
    for i in ids:
        m = _PEAK.match(str(i))
        if m:
            chrom.append(m["chrom"]); start.append(int(m["start"])); end.append(int(m["end"]))
        else:
            chrom.append(None); start.append(-1); end.append(-1)
    return pd.DataFrame({"chrom": chrom, "start": start, "end": end})


def tile_keys(ids, tile: int = 500) -> np.ndarray:
    """Integer key per peak = chrom*1e9 + tile start (midpoint rule); -1 for non-autosomal/unparsable."""
    p = parse_peaks(ids)
    ok = p["chrom"].isin(AUTOSOMES).to_numpy()
    mid = ((p["start"] + p["end"]) // 2).to_numpy()
    chrom = pd.to_numeric(p["chrom"], errors="coerce").fillna(-1).to_numpy().astype(np.int64)
    key = chrom * _CHROM_BASE + (mid // tile) * tile
    key = key.astype(np.int64)
    key[~ok] = -1
    return key


def consensus_tiles(per_sample_keys: dict, min_samples: int = 3) -> np.ndarray:
    """Sorted tile keys present (via any peak) in at least ``min_samples`` samples."""
    counts: dict = {}
    for k in per_sample_keys.values():
        for u in np.unique(k[k >= 0]):
            counts[u] = counts.get(u, 0) + 1
    keep = [k for k, n in counts.items() if n >= min_samples]
    return np.array(sorted(keep), dtype=np.int64)


def tile_ids(keys: np.ndarray, tile: int = 500) -> list:
    return [f"chr{k // _CHROM_BASE}:{k % _CHROM_BASE}-{k % _CHROM_BASE + tile}" for k in keys]


def peak_to_tile_index(peak_keys: np.ndarray, consensus: np.ndarray) -> np.ndarray:
    """Column index into ``consensus`` for every peak, -1 if the peak's tile is not in the consensus."""
    pos = np.searchsorted(consensus, peak_keys)
    pos[pos >= len(consensus)] = 0
    hit = (consensus[pos] == peak_keys) & (peak_keys >= 0)
    return np.where(hit, pos, -1)


def cells_by_tiles(x_peaks_by_cells: sp.spmatrix, tile_idx: np.ndarray, n_tiles: int) -> sp.csr_matrix:
    """(peaks x cells) counts -> (cells x consensus tiles); peaks sharing a tile are summed."""
    rows = np.flatnonzero(tile_idx >= 0)
    m = sp.csr_matrix((np.ones(len(rows)), (rows, tile_idx[rows])), shape=(x_peaks_by_cells.shape[0], n_tiles))
    return (x_peaks_by_cells.tocsc().T @ m).tocsr()


def split_cell_id(index: pd.Series, sample: pd.Series) -> pd.Series:
    """h5ad index '<sample>_<barcode>' -> barcode (raises if a prefix does not match)."""
    out = []
    for cid, s in zip(index, sample):
        pre = f"{s}_"
        if not cid.startswith(pre):
            raise ValueError(f"cell id {cid!r} does not start with its sample {s!r}")
        out.append(cid[len(pre):])
    return pd.Series(out, index=index.index)
