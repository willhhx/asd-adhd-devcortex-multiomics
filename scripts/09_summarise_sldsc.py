#!/usr/bin/env python3
"""Step 9: collect S-LDSC output -> tau*, FDR, paired ASD-vs-ADHD contrasts, stage contrasts.

Expects <sldsc-dir>/<trait>/<state_file>.results and .part_delete (ldsc.py --h2
--print-coefficients --print-delete-vals, last reference annotation = the state).

Analyses (all from the same jackknife delete values):
  1. tau* per trait x state; one-sided z test for tau* > 0; BH-FDR across states within trait.
  2. Paired contrast tau*(A) - tau*(B) for every state; BH-FDR across states.
  3. Stage-shape contrasts within a cell type, for adjacent stages s -> s+1:
     d_trait = tau*(s+1) - tau*(s); contrast d_A - d_B (difference in differences).
  4. Control-trait check: number of states reaching FDR<5% per trait.
  5. Power table: SE, minimum detectable difference at 80% power, effective-N ratio.
Outputs are TSV files in --out.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import jackknife as jk  # noqa: E402
from devcortex import power  # noqa: E402
from devcortex.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sldsc-dir", required=True)
    ap.add_argument("--peakset-summary", required=True, help="peakset_summary.tsv from step 4")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(args.out, exist_ok=True)
    summ = pd.read_csv(args.peakset_summary, sep="\t")
    state_files = dict(zip(summ["state"], summ["file"]))
    stage_order = cfg["atlas"]["stage_order"]
    traits = list(cfg["traits"])

    est, dels, rows = {}, {}, []
    for trait in traits:
        for state, fname in state_files.items():
            res = os.path.join(args.sldsc_dir, trait, fname + ".results")
            dl = os.path.join(args.sldsc_dir, trait, fname + ".part_delete")
            if not (os.path.exists(res) and os.path.exists(dl)):
                print(f"missing: {trait} / {state}", file=sys.stderr)
                continue
            s, d = jk.summarise_one_regression(res, dl)
            ct, stage = state.split("|", 1)
            rows.append({"trait": trait, "state": state, "cell_type": ct, "stage": stage, **s})
            est[(trait, state)] = s["tau_star"]
            dels[(trait, state)] = d
    tau = pd.DataFrame(rows)
    if tau.empty:
        sys.exit("no S-LDSC results found")
    tau["fdr_within_trait"] = np.nan
    for trait, idx in tau.groupby("trait").groups.items():
        tau.loc[idx, "fdr_within_trait"] = jk.benjamini_hochberg(tau.loc[idx, "tau_star_p_one_sided"])
    tau.to_csv(os.path.join(args.out, "tau_star_by_state.tsv"), sep="\t", index=False)
    blocks = tau["n_blocks"].unique()
    if len(blocks) != 1:
        sys.exit(f"traits/states have different numbers of jackknife blocks {blocks}: "
                 "run step 7b so all traits share one SNP set")

    # ---- paired contrasts per state
    pairs = cfg["contrasts"]
    crow = []
    for a, b in pairs:
        for state in state_files:
            if (a, state) in est and (b, state) in est:
                r = jk.contrast(est, dels, {(a, state): 1, (b, state): -1})
                naive = jk.unpaired_difference(est[(a, state)], jk.jackknife_se(dels[(a, state)]),
                                               est[(b, state)], jk.jackknife_se(dels[(b, state)]))
                ct, stage = state.split("|", 1)
                crow.append({"contrast": f"{a}-{b}", "state": state, "cell_type": ct, "stage": stage,
                             f"tau_star_{a}": est[(a, state)], f"tau_star_{b}": est[(b, state)],
                             **r, "p_unpaired_sensitivity": naive["p_two_sided"]})
    con = pd.DataFrame(crow)
    if not con.empty:
        con["fdr"] = np.nan
        for name, idx in con.groupby("contrast").groups.items():
            con.loc[idx, "fdr"] = jk.benjamini_hochberg(con.loc[idx, "p_two_sided"])
        con.to_csv(os.path.join(args.out, "paired_contrasts_by_state.tsv"), sep="\t", index=False)

    # ---- stage-shape contrasts (difference in differences between adjacent stages)
    srow = []
    for a, b in pairs:
        for ct in sorted(tau["cell_type"].unique()):
            for s0, s1 in zip(stage_order[:-1], stage_order[1:]):
                k = [(a, f"{ct}|{s0}"), (a, f"{ct}|{s1}"), (b, f"{ct}|{s0}"), (b, f"{ct}|{s1}")]
                if not all(x in est for x in k):
                    continue
                w = {k[1]: 1, k[0]: -1, k[3]: -1, k[2]: 1}
                r = jk.contrast(est, dels, w)
                srow.append({"contrast": f"{a}-{b}", "cell_type": ct, "from_stage": s0, "to_stage": s1, **r})
    sdf = pd.DataFrame(srow)
    if not sdf.empty:
        sdf["fdr"] = np.nan
        for name, idx in sdf.groupby("contrast").groups.items():
            sdf.loc[idx, "fdr"] = jk.benjamini_hochberg(sdf.loc[idx, "p_two_sided"])
        sdf.to_csv(os.path.join(args.out, "stage_shape_contrasts.tsv"), sep="\t", index=False)

    # ---- controls
    ctrl = (tau.assign(sig=tau["fdr_within_trait"] < 0.05)
               .groupby("trait").agg(n_states=("state", "size"), n_fdr05=("sig", "sum"),
                                     max_tau_star=("tau_star", "max")).reset_index())
    ctrl["role"] = ctrl["trait"].map({t: v.get("role", "") for t, v in cfg["traits"].items()})
    ctrl.to_csv(os.path.join(args.out, "control_summary.tsv"), sep="\t", index=False)

    # ---- power table
    prow = []
    neff = {t: (v["n_eff"] if v.get("n_eff") else power.effective_n(v["n_cases"], v["n_controls"]))
            for t, v in cfg["traits"].items() if v.get("n_eff") or (v.get("n_cases") and v.get("n_controls"))}
    for trait in traits:
        sub = tau[tau["trait"] == trait]
        if sub.empty:
            continue
        se = float(sub["tau_star_se"].median())
        prow.append({"trait": trait, "n_eff": neff.get(trait, np.nan), "median_tau_star_se": se,
                     "min_detectable_tau_star_80pct": power.minimum_detectable_difference(se, alpha=0.05),
                     "min_detectable_tau_star_80pct_bonferroni": power.minimum_detectable_difference(
                         se, alpha=power.bonferroni_alpha(0.05, len(sub)))})
    pd.DataFrame(prow).to_csv(os.path.join(args.out, "power_table.tsv"), sep="\t", index=False)
    print(f"tau* rows: {len(tau)}; contrasts: {len(con)}; stage-shape: {len(sdf)}")
    print(ctrl.to_string(index=False))


if __name__ == "__main__":
    main()
