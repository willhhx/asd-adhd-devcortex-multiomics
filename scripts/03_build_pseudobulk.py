#!/usr/bin/env python3
"""Step 3: nuclei x features count matrix -> balanced pseudobulk per cell state.

Reads either an AnnData .h5ad (``--format h5ad``; needs the ``anndata`` package) or the
plain-file layout used by the synthetic test (``--format npz``: counts.npz [cells x
features, scipy sparse], obs.tsv, features.txt).

Filters thin states (nuclei/donor thresholds from the config), down-samples every state to
the same number of nuclei, sums counts, and writes:
  <out>/<modality>_counts.tsv.gz      features x states raw pseudobulk counts
  <out>/state_adequacy.tsv            nuclei, donors, samples per state and whether kept
  <out>/state_donor_share.tsv         largest single-donor share per kept state
  <out>/<modality>_features.tsv       feature table (chrom,start,end for peaks)
"""
import argparse
import os
import sys

import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import peaks, pseudobulk, stagebins  # noqa: E402
from devcortex.config import ensure_dir, load_config  # noqa: E402


def load_h5ad(path):
    try:
        import anndata as ad
    except ImportError:
        sys.exit("anndata is required for --format h5ad (conda env: envs/devcortex.yml)")
    adata = ad.read_h5ad(path)
    return adata.X, adata.obs.copy(), list(adata.var_names)


def load_npz(prefix):
    x = sp.load_npz(prefix + ".counts.npz")
    obs = pd.read_csv(prefix + ".obs.tsv", sep="\t", index_col=0)
    with open(prefix + ".features.txt") as fh:
        feats = [ln.strip() for ln in fh if ln.strip()]
    return x, obs, feats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--format", choices=["h5ad", "npz"], required=True)
    ap.add_argument("--input", required=True, help=".h5ad path, or the prefix for the npz layout")
    ap.add_argument("--modality", default="atac", choices=["atac", "rna"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    a = cfg["atlas"]
    x, obs, feats = load_h5ad(args.input) if args.format == "h5ad" else load_npz(args.input)
    obs = obs.reset_index(drop=True)
    for col in (a["cell_type_col"], a["stage_col"], a["donor_col"], a["sample_col"]):
        if col not in obs.columns:
            sys.exit(f"obs is missing column {col!r}; available: {list(obs.columns)}")

    obs["cell_type"] = obs[a["cell_type_col"]].astype(str)
    obs["stage"] = obs[a["stage_col"]].astype(str)
    obs["donor"] = obs[a["donor_col"]].astype(str)
    obs["sample"] = obs[a["sample_col"]].astype(str)

    tab, obs_keep = stagebins.keep_adequate(
        obs, min_nuclei=a["min_nuclei_per_state"], min_donors=a["min_donors_per_state"])
    ensure_dir(args.out)
    tab.to_csv(os.path.join(args.out, "state_adequacy.tsv"), sep="\t", index=False)
    n_dropped = int((~tab["keep"]).sum())
    print(f"states: {len(tab)} total, {int(tab['keep'].sum())} kept, {n_dropped} dropped (thin)")

    picked = stagebins.balanced_nuclei_sample(obs_keep, a["nuclei_per_state"], seed=cfg.get("seed", 0))
    stagebins.donor_imbalance(obs.loc[picked]).to_csv(
        os.path.join(args.out, "state_donor_share.tsv"), sep="\t", index=False)

    labels = pseudobulk.state_labels(obs, "cell_type", "stage")
    pb = pseudobulk.aggregate_pseudobulk(x, labels, keep_cells=picked)
    pb.index = pd.Index(feats, name="feature")
    pb.to_csv(os.path.join(args.out, f"{args.modality}_counts.tsv.gz"), sep="\t", compression="gzip")
    if args.modality == "atac":
        peaks.parse_peak_ids(feats).to_csv(os.path.join(args.out, "atac_features.tsv"), sep="\t")
    print(f"wrote pseudobulk {pb.shape[0]} features x {pb.shape[1]} states -> {args.out}")


if __name__ == "__main__":
    main()
