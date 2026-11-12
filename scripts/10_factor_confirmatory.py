#!/usr/bin/env python3
"""Step 10: pre-specified confirmatory tests (H1, H2a, H2b) and control checks for the factor analysis.

Run AFTER step 9 on the pooled factor regressions:
  python scripts/10_factor_confirmatory.py --config config/config_factors_pooled.yaml \
    --sldsc-dir work/pooled_factors/sldsc --peakset-summary work/pooled/peaksets_hg19/peakset_summary.tsv \
    --tau-table work/pooled_factors/results/tau_star_by_state.tsv --out work/pooled_factors/results
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import factor_tests as ft, jackknife as jk  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sldsc-dir", required=True)
    ap.add_argument("--peakset-summary", required=True)
    ap.add_argument("--tau-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gate-se", type=float, default=0.30)
    a = ap.parse_args()
    summ = pd.read_csv(a.peakset_summary, sep="\t")
    est, dels = {}, {}
    for trait in ("F2", "F3", "F4"):
        for state, fname in zip(summ["state"], summ["file"]):
            res = os.path.join(a.sldsc_dir, trait, fname + ".results")
            dl = os.path.join(a.sldsc_dir, trait, fname + ".part_delete")
            if os.path.exists(res) and os.path.exists(dl):
                s, d = jk.summarise_one_regression(res, dl)
                est[(trait, state)], dels[(trait, state)] = s["tau_star"], d
    out = ft.confirmatory(est, dels, gate_se=a.gate_se)
    os.makedirs(a.out, exist_ok=True)
    out.to_csv(os.path.join(a.out, "confirmatory_tests.tsv"), sep="\t", index=False)
    print(out.round(4).to_string(index=False))

    tau = pd.read_csv(a.tau_table, sep="\t")
    neuronal = tau["cell_type"].isin(["Excitatory", "Inhibitory"])
    ctrl = []
    f2 = tau[(tau.trait == "F2") & neuronal & (tau.stage == "Postnatal")]
    ctrl.append({"check": "F2 postnatal neurons enriched (positive control)",
                 "n_groups": len(f2), "n_fdr05": int((f2["fdr_within_trait"] < 0.05).sum()),
                 "min_fdr": float(f2["fdr_within_trait"].min()) if len(f2) else float("nan")})
    h = tau[(tau.trait == "HEIGHT") & neuronal]
    ctrl.append({"check": "HEIGHT neuronal groups not enriched (specificity control)",
                 "n_groups": len(h), "n_fdr05": int((h["fdr_within_trait"] < 0.05).sum()),
                 "min_fdr": float(h["fdr_within_trait"].min()) if len(h) else float("nan")})
    c = pd.DataFrame(ctrl)
    c.to_csv(os.path.join(a.out, "confirmatory_controls.tsv"), sep="\t", index=False)
    print(c.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
