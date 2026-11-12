"""Developmental stage binning and donor/nuclei balance checks.

Primary analysis uses the age groups assigned by the atlas authors (metadata column
``age_group``). This module supplies (i) a fallback assignment from post-conception age,
(ii) a donor/nuclei adequacy table so thin cell-state x stage combinations are dropped
before any enrichment test, and (iii) alternative cut-points for the stage-bin
sensitivity analysis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Post-conception weeks (PCW). Conventional trimester boundaries; the primary analysis
# uses the authors' labels, so these only matter for the fallback / sensitivity runs.
DEFAULT_EDGES_PCW = {
    "first_trimester": (0, 14),
    "second_trimester": (14, 28),
    "third_trimester": (28, 40),
    "infancy": (40, 40 + 52 * 2),        # 0-2 years postnatal
    "childhood_adolescence": (40 + 52 * 2, 40 + 52 * 25),
}


def age_to_pcw(age_value: float, unit: str) -> float:
    """Convert an age to post-conception weeks. unit in {'pcw','gw','years','months','days'}.

    Gestational weeks (GW, from last menstrual period) are ~2 weeks greater than PCW.
    """
    u = unit.lower()
    if u == "pcw":
        return float(age_value)
    if u == "gw":
        return float(age_value) - 2.0
    if u == "years":
        return 40.0 + 52.0 * float(age_value)
    if u == "months":
        return 40.0 + 52.0 * float(age_value) / 12.0
    if u == "days":
        return 40.0 + float(age_value) / 7.0
    raise ValueError(f"unknown unit {unit!r}")


def assign_stage(pcw: float, edges: dict | None = None) -> str:
    edges = edges or DEFAULT_EDGES_PCW
    for name, (lo, hi) in edges.items():
        if lo <= pcw < hi:
            return name
    return "out_of_range"


def shifted_edges(shift_weeks: float, edges: dict | None = None) -> dict:
    """Shift the *fetal* internal boundaries by ``shift_weeks`` (sensitivity analysis).

    Only boundaries between fetal stages move; 0 and the birth boundary (40 PCW) stay put.
    """
    edges = dict(edges or DEFAULT_EDGES_PCW)
    names = list(edges)
    fetal = [n for n in names if edges[n][1] <= 40]
    out = {}
    for i, n in enumerate(names):
        lo, hi = edges[n]
        if n in fetal:
            if i > 0 and names[i - 1] in fetal:
                lo = lo + shift_weeks
            if i + 1 < len(names) and names[i + 1] in fetal:
                hi = hi + shift_weeks
        out[n] = (lo, hi)
    return out


def adequacy_table(obs: pd.DataFrame, state_cols=("cell_type", "stage"),
                   donor_col="donor", sample_col="sample") -> pd.DataFrame:
    """Nuclei, donors and samples per cell-state (x stage)."""
    g = obs.groupby(list(state_cols), observed=True)
    out = g.agg(n_nuclei=(donor_col, "size"), n_donors=(donor_col, "nunique"),
                n_samples=(sample_col, "nunique")).reset_index()
    return out


def keep_adequate(obs: pd.DataFrame, min_nuclei: int = 200, min_donors: int = 2,
                  state_cols=("cell_type", "stage"), donor_col="donor",
                  sample_col="sample") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (adequacy table with a boolean ``keep`` column, filtered obs)."""
    tab = adequacy_table(obs, state_cols, donor_col, sample_col)
    tab["keep"] = (tab["n_nuclei"] >= min_nuclei) & (tab["n_donors"] >= min_donors)
    kept = tab.loc[tab["keep"], list(state_cols)]
    # boolean mask (NOT merge): merge resets the index and would break downstream row selection
    key = pd.MultiIndex.from_frame(obs[list(state_cols)])
    filt = obs[key.isin(pd.MultiIndex.from_frame(kept))]
    return tab, filt


def balanced_nuclei_sample(obs: pd.DataFrame, n_per_state: int, seed: int = 0,
                           state_cols=("cell_type", "stage")) -> pd.Index:
    """Index of nuclei after down-sampling every state to ``n_per_state``.

    Equal nuclei per state removes the dependence of peak detection (and thus of the
    size/composition of specific-peak sets) on how many nuclei a state happens to have.
    States with fewer nuclei than ``n_per_state`` are kept whole (and should have been
    removed by :func:`keep_adequate`).
    """
    rng = np.random.default_rng(seed)
    picks = []
    for _, sub in obs.groupby(list(state_cols), observed=True):
        idx = sub.index.to_numpy()
        if len(idx) > n_per_state:
            idx = rng.choice(idx, size=n_per_state, replace=False)
        picks.append(idx)
    return pd.Index(np.concatenate(picks))


def donor_imbalance(obs: pd.DataFrame, state_cols=("cell_type", "stage"),
                    donor_col="donor") -> pd.DataFrame:
    """Largest single-donor share of nuclei within each state (donor-effect diagnostic)."""
    def share(sub):
        c = sub[donor_col].value_counts(normalize=True)
        return pd.Series({"top_donor": c.index[0], "top_donor_share": float(c.iloc[0])})
    return obs.groupby(list(state_cols), observed=True).apply(share, include_groups=False).reset_index()
