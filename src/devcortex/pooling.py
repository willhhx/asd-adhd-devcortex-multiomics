"""Pool cell states into broader cell classes x stage groups (secondary, post hoc analysis).

Why: per-state annotations are small and the ASD GWAS is underpowered (see power table), so the
primary 94-state analysis cannot resolve ASD vs ADHD. Pooling members' pseudobulk counts gives
fewer, better-powered annotations. Peak sets keep the primary rule and size (top-N log2FC peaks
specific to the pooled group vs the other pooled groups), so annotation size stays constant.
"""
from __future__ import annotations

import pandas as pd

DEFAULT_STAGE_GROUPS = {
    "First_trimester": "Early_fetal",
    "Second_trimester": "Mid_late_fetal",
    "Third_trimester": "Mid_late_fetal",
    "Infancy": "Postnatal",
    "Adolescence": "Postnatal",
}
STAGE_GROUP_ORDER = ["Early_fetal", "Mid_late_fetal", "Postnatal"]


def classify_cell_type(name: str) -> str:
    """Rule-based class from the Wang 2025 `Type_updated` label. Returns 'Other' if unmatched
    (Other is never pooled; the script reports every unmatched label so none is silent)."""
    n = name.lower()
    if n.startswith("unknown"):
        return "Other"
    if "ipc" in n or n.startswith("rg") or "radial" in n or "cycling" in n or "progenitor" in n:
        return "Progenitor"
    if n.startswith("en-"):
        return "Excitatory"
    if n.startswith("in-"):
        return "Inhibitory"
    if any(k in n for k in ("astrocyte", "opc", "oligo", "microglia")):
        return "Glia"
    return "Other"


def pool_counts(counts: pd.DataFrame, stage_groups: dict | None = None):
    """counts: features x states ('cell_type|stage'). Returns (pooled counts, membership table).
    Pooled column = '<Class>|<StageGroup>'; states in class 'Other' or an unmapped stage are
    left out and listed in the membership table with pooled='' ."""
    sg = stage_groups or DEFAULT_STAGE_GROUPS
    rows, groups = [], {}
    for col in counts.columns:
        ct, stage = col.split("|", 1)
        cls = classify_cell_type(ct)
        grp = sg.get(stage)
        pooled = f"{cls}|{grp}" if cls != "Other" and grp else ""
        rows.append({"state": col, "cell_type": ct, "stage": stage, "class": cls,
                     "stage_group": grp or "", "pooled": pooled})
        if pooled:
            groups.setdefault(pooled, []).append(col)
    out = pd.DataFrame({g: counts[cols].sum(axis=1) for g, cols in groups.items()})
    return out, pd.DataFrame(rows)
