"""Orchestration of the Wang-atlas build with injectable readers (so it is testable without HDF5)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import stagebins, wang_atlas as wa


def select_cells(obs: pd.DataFrame, min_nuclei: int, min_donors: int, n_per_state: int, seed: int):
    """Adequacy filter + equal nuclei per (cell_type, stage); same logic as step 3."""
    tab, kept = stagebins.keep_adequate(obs, min_nuclei=min_nuclei, min_donors=min_donors)
    picked = stagebins.balanced_nuclei_sample(kept, n_per_state, seed=seed)
    return tab, obs.loc[picked]


def build(obs_sel: pd.DataFrame, samples: list, read_features, read_matrix, tile=500, min_samples=3, dry_run=False, log=print):
    """obs_sel: selected nuclei, index '<sample>_<barcode>', columns sample (+labels).
    read_features(sample) -> peak ids (Peaks only, matrix row order); read_matrix(sample, wanted_barcodes) ->
    (barcodes_found, peak_ids, csc peaks x found_cells)."""
    keys = {s: wa.tile_keys(read_features(s), tile) for s in samples}
    cons = wa.consensus_tiles(keys, min_samples)
    log(f"consensus tiles: {len(cons)} (>= {min_samples} samples, {tile} bp)")
    barcode = wa.split_cell_id(pd.Series(obs_sel.index, index=obs_sel.index), obs_sel["sample"])
    rows, report = [], []
    for s in samples:
        sel = obs_sel.index[(obs_sel["sample"] == s).to_numpy()]
        if len(sel) == 0:
            report.append((s, 0, 0)); continue
        want = barcode.loc[sel]
        if dry_run:
            found, _, _ = read_matrix(s, list(want), header_only=True)
            report.append((s, len(sel), len(found))); continue
        found, peak_ids, x = read_matrix(s, list(want))
        pos = {b: i for i, b in enumerate(found)}
        ok = [c for c, b in zip(sel, want.loc[sel]) if b in pos]
        cols = [pos[barcode.loc[c]] for c in ok]
        idx = wa.peak_to_tile_index(wa.tile_keys(peak_ids, tile), cons)
        blk = wa.cells_by_tiles(x[:, cols], idx, len(cons))
        rows.append((ok, blk))
        report.append((s, len(sel), len(ok)))
    rep = pd.DataFrame(report, columns=["sample", "selected", "matched"])
    tot_sel, tot_match = int(rep["selected"].sum()), int(rep["matched"].sum())
    log(f"barcode match: {tot_match}/{tot_sel} selected nuclei")
    if tot_sel and tot_match < 0.99 * tot_sel:
        raise RuntimeError(f"only {tot_match}/{tot_sel} nuclei matched; check barcode format / sample names\n{rep}")
    if dry_run:
        return None, None, tile_ids_or_empty(cons, tile), rep
    cells = [c for ok, _ in rows for c in ok]
    mat = sp.vstack([b for _, b in rows]).tocsr()
    return mat, obs_sel.loc[cells], wa.tile_ids(cons, tile), rep


def tile_ids_or_empty(cons, tile):
    return wa.tile_ids(cons, tile)
