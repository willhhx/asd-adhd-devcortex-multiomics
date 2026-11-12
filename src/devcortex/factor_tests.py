"""Pre-specified confirmatory contrasts of docs/analysis_plan_v3.md (H1, H2a, H2b) and the power gate.

Groups are the pooled class x stage-group annotations ('Excitatory|Early_fetal', ...). Every contrast
is a linear combination of tau* estimates, tested with the paired block jackknife (jackknife.contrast).
"""
from __future__ import annotations

import pandas as pd
from scipy import stats

from . import jackknife as jk

NEURONAL_FETAL = [f"{c}|{s}" for c in ("Excitatory", "Inhibitory") for s in ("Early_fetal", "Mid_late_fetal")]
NEURONAL_POST = [f"{c}|Postnatal" for c in ("Excitatory", "Inhibitory")]
NEURONAL = NEURONAL_FETAL + NEURONAL_POST
GLIAL = ["Glia|Mid_late_fetal", "Glia|Postnatal"]

# (trait, groups, sign): each term contributes sign * mean(tau*) over its groups
HYPOTHESES = {
    "H1_timing (F3-F2 fetal neurons minus F3-F2 postnatal neurons > 0)": [
        ("F3", NEURONAL_FETAL, +1), ("F2", NEURONAL_FETAL, -1),
        ("F3", NEURONAL_POST, -1), ("F2", NEURONAL_POST, +1)],
    "H2a_class (F2 neuronal minus glial > 0)": [("F2", NEURONAL, +1), ("F2", GLIAL, -1)],
    "H2b_class (F4 glial minus neuronal > 0)": [("F4", GLIAL, +1), ("F4", NEURONAL, -1)],
}


def build_weights(terms) -> dict:
    w: dict = {}
    for trait, groups, sign in terms:
        for g in groups:
            w[(trait, g)] = w.get((trait, g), 0.0) + sign / len(groups)
    return w


def confirmatory(est: dict, dels: dict, gate_se: float = 0.30, alpha: float = 0.05) -> pd.DataFrame:
    """One row per pre-specified hypothesis. One-sided p; Bonferroni alpha/len(HYPOTHESES).
    Verdicts: 'supported', 'not supported (informative)', 'uninformative: underpowered (bound only)'."""
    n = len(HYPOTHESES)
    a = alpha / n
    rows = []
    for name, terms in HYPOTHESES.items():
        w = build_weights(terms)
        missing = [k for k in w if k not in est]
        if missing:
            rows.append({"hypothesis": name, "verdict": f"not computable: missing {missing[:3]}"})
            continue
        r = jk.contrast(est, dels, w)
        p1 = float(stats.norm.sf(r["z"]))
        mdd = float((stats.norm.isf(a) + stats.norm.ppf(0.8)) * r["se"])
        informative = r["se"] <= gate_se
        sig = p1 < a
        verdict = ("supported" if (sig and informative) else
                   "not supported (informative)" if informative else
                   "uninformative: underpowered (bound only)")
        rows.append({"hypothesis": name, "estimate": r["estimate"], "se": r["se"], "z": r["z"],
                     "ci_low": r["ci_low"], "ci_high": r["ci_high"], "p_one_sided": p1,
                     "alpha_bonferroni": a, "gate_se": gate_se, "informative": informative,
                     "min_detectable_80pct": mdd, "verdict": verdict})
    return pd.DataFrame(rows)
