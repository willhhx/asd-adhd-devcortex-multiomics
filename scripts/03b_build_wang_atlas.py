#!/usr/bin/env python3
"""Step 3b: Wang et al. 2025 deposit -> nuclei x consensus-peak-tile matrix in the npz layout of step 3.

Inputs: the CELLxGENE .h5ad (cell labels only; RNA matrix is not read) and the 38 per-sample 10x
filtered_feature_bc_matrix .h5 files from Dryad (peaks + counts).
Output prefix <out>/atlas: .counts.npz (cells x tiles), .obs.tsv, .features.txt, plus build_report.tsv.
Then run step 3 with --format npz --input <out>/atlas.
Use --dry-run first: reads labels + peak lists only and checks barcode matching (a minute or two).
Needs h5py.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import wang_build  # noqa: E402
from devcortex.config import ensure_dir, load_config  # noqa: E402

try:
    import h5py
except ImportError:
    sys.exit("h5py is required: pip install h5py")

COLS = {"cell_type": "Type_updated", "stage": "Group", "donor": "donor_id", "sample": "Sample_ID", "region": "Region"}


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def read_obs(h5ad):
    with h5py.File(h5ad, "r") as f:
        o = f["obs"]
        out = {"cell": _dec(o["_index"][:])}
        for new, old in COLS.items():
            g = o[old]
            out[new] = _dec(g["categories"][:])[g["codes"][:]] if "categories" in g else _dec(g[:])
    return pd.DataFrame(out).set_index("cell")


def make_readers(h5dir):
    path = lambda s: os.path.join(h5dir, f"{s}_filtered_feature_bc_matrix.h5")  # noqa: E731

    def peak_mask(f):
        return _dec(f["matrix"]["features"]["feature_type"][:]) == "Peaks"

    def read_features(s):
        with h5py.File(path(s), "r") as f:
            return list(_dec(f["matrix"]["features"]["id"][:])[peak_mask(f)])

    def read_matrix(s, wanted, header_only=False):
        with h5py.File(path(s), "r") as f:
            m = f["matrix"]
            bc = _dec(m["barcodes"][:])
            idx = {b: i for i, b in enumerate(bc)}
            found = [b for b in wanted if b in idx]
            if header_only:
                return found, None, None
            shape = tuple(int(v) for v in m["shape"][:])             # (features, barcodes)
            x = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=shape)
            mask = peak_mask(f)
            ids = list(_dec(m["features"]["id"][:])[mask])
            cols = [idx[b] for b in found]
            return found, ids, x[mask][:, cols].tocsc()
    return read_features, read_matrix


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--h5-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tile", type=int, default=500)
    ap.add_argument("--min-samples", type=int, default=3)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cfg = load_config(a.config)
    at = cfg["atlas"]
    obs = read_obs(a.h5ad)
    print(f"{len(obs)} labelled nuclei; stages: {obs['stage'].value_counts().to_dict()}")
    tab, sel = wang_build.select_cells(obs, at["min_nuclei_per_state"], at["min_donors_per_state"],
                                       at["nuclei_per_state"], cfg.get("seed", 0))
    print(f"states kept {int(tab['keep'].sum())}/{len(tab)}; selected nuclei {len(sel)}")
    files = sorted(glob.glob(os.path.join(a.h5_dir, "*_filtered_feature_bc_matrix.h5")))
    samples = [os.path.basename(f).replace("_filtered_feature_bc_matrix.h5", "") for f in files]
    missing = sorted(set(sel["sample"]) - set(samples))
    if missing:
        sys.exit(f"no .h5 file for samples: {missing}")
    rf, rm = make_readers(a.h5_dir)
    mat, obs_out, feats, rep = wang_build.build(sel, samples, rf, rm, a.tile, a.min_samples, a.dry_run)
    ensure_dir(a.out)
    rep.to_csv(os.path.join(a.out, "build_report.tsv"), sep="\t", index=False)
    if a.dry_run:
        print(rep.to_string(index=False)); print("dry run OK; nothing else written"); return
    sp.save_npz(os.path.join(a.out, "atlas.counts.npz"), mat)
    obs_out.to_csv(os.path.join(a.out, "atlas.obs.tsv"), sep="\t")
    with open(os.path.join(a.out, "atlas.features.txt"), "w") as fh:
        fh.write("\n".join(feats) + "\n")
    print(f"wrote {mat.shape[0]} nuclei x {mat.shape[1]} tiles -> {a.out}/atlas.*")


if __name__ == "__main__":
    main()
