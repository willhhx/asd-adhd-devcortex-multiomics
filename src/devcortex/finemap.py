"""Single-causal-variant Bayesian fine-mapping with optional functional priors.

Wakefield approximate Bayes factors (ABF) per variant within a locus, posterior inclusion
probabilities (PIP), and credible sets. A functional prior can up-weight variants that
fall in cell-state-specific peaks, with the up-weight derived from the S-LDSC enrichment
of that state (a simplified, transparent analogue of PolyFun priors). This is the
primary engine for the loci that harbour one dominant signal; the SuSiE/PolyFun route in
scripts/11b_susie_polyfun.R handles multiple signals per locus and needs an LD matrix.

Caveat reported in the paper: the single-causal assumption can under-cover when a locus
carries several independent signals.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import logsumexp

# prior SD of the true effect on the log-odds scale (the coloc default for case-control).
DEFAULT_SD_PRIOR = 0.2


def log_abf(beta, se, sd_prior: float = DEFAULT_SD_PRIOR) -> np.ndarray:
    """Wakefield log ABF in favour of association (alternative vs null)."""
    beta = np.asarray(beta, float)
    v = np.asarray(se, float) ** 2
    w = sd_prior ** 2
    r = w / (v + w)
    z2 = beta ** 2 / v
    return 0.5 * (np.log(1.0 - r) + r * z2)


def functional_prior(in_annot, enrichment: float, cap: float = 10.0) -> np.ndarray:
    """Normalised prior: variants in the annotation get weight ``min(enrichment, cap)``,
    others weight 1. ``enrichment`` < 1 is floored at 1 (we never penalise)."""
    e = float(np.clip(enrichment, 1.0, cap))
    w = np.where(np.asarray(in_annot, bool), e, 1.0)
    return w / w.sum()


def pip(beta, se, prior=None, sd_prior: float = DEFAULT_SD_PRIOR) -> np.ndarray:
    labf = log_abf(beta, se, sd_prior)
    lp = np.zeros_like(labf) if prior is None else np.log(np.asarray(prior, float))
    lp = lp - logsumexp(lp)                      # ensure a proper prior
    post = labf + lp
    return np.exp(post - logsumexp(post))


def credible_set(pips, level: float = 0.95) -> np.ndarray:
    """Indices of the smallest set whose PIPs sum to >= level."""
    order = np.argsort(-np.asarray(pips), kind="mergesort")
    cum = np.cumsum(np.asarray(pips)[order])
    k = int(np.searchsorted(cum, level - 1e-12) + 1)
    return order[:k]


def define_loci(gwas: pd.DataFrame, p_thresh: float = 5e-8, window: int = 500_000,
                chrom_col="CHR", pos_col="BP", p_col="P") -> pd.DataFrame:
    """Greedy independent loci: the most significant SNP, then everything within +-window is
    absorbed, repeat. Returns one row per locus (lead SNP, bounds)."""
    sig = gwas[gwas[p_col] < p_thresh].sort_values(p_col)
    taken = pd.Series(False, index=sig.index)
    loci = []
    for idx, row in sig.iterrows():
        if taken[idx]:
            continue
        near = (sig[chrom_col] == row[chrom_col]) & ((sig[pos_col] - row[pos_col]).abs() <= window)
        taken |= near
        loci.append({"locus": len(loci) + 1, "chrom": row[chrom_col], "lead_pos": int(row[pos_col]),
                     "lead_snp": row.get("SNP", idx), "lead_p": float(row[p_col]),
                     "start": int(row[pos_col] - window), "end": int(row[pos_col] + window)})
    return pd.DataFrame(loci)


def finemap_locus(df: pd.DataFrame, prior=None, level: float = 0.95,
                  sd_prior: float = DEFAULT_SD_PRIOR) -> pd.DataFrame:
    """df needs BETA and SE. Adds ``pip`` and ``in_credible_set``."""
    out = df.copy()
    out["pip"] = pip(df["BETA"], df["SE"], prior, sd_prior)
    cs = credible_set(out["pip"].to_numpy(), level)
    out["in_credible_set"] = False
    out.iloc[cs, out.columns.get_loc("in_credible_set")] = True
    return out
