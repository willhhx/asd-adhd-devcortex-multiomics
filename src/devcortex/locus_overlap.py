"""Locus-level test: do fine-mapped causal-variant posterior probabilities concentrate in the
peaks of a cell-state group more than chance? (Pre-specified analysis L1.)

Statistic: sum over loci of the PIP mass lying in the group's peaks (flat-prior PIPs, so the
test is not circular with the S-LDSC result that selects states). Null: within each locus the
peak indicator is cyclically shifted along the position-sorted SNPs, which keeps the peak
density, clustering and the PIP vector fixed. All cyclic shifts are evaluated exactly with one
FFT per locus; the null then draws one random shift per locus per permutation.
"""
from __future__ import annotations

import numpy as np


def _circ_cross(pip: np.ndarray, inpeak: np.ndarray) -> np.ndarray:
    """c[k] = sum_i pip[i] * inpeak[(i + k) % n] for all k (exact, via FFT)."""
    n = len(pip)
    f = np.fft.rfft(pip[::-1], n) * np.fft.rfft(inpeak.astype(float), n)
    full = np.fft.irfft(f, n)          # circular convolution of reversed pip with inpeak
    return full[(n - 1 + np.arange(n)) % n]   # re-index: c[k] = full[(n-1+k) % n]


def pip_mass_enrichment(pips: list, inpeaks: list, n_perm: int = 5000, seed: int = 0) -> dict:
    """pips[l] (sums to ~1) and inpeaks[l] (bool) are position-sorted vectors for locus l."""
    rng = np.random.default_rng(seed)
    obs = float(sum(np.dot(p, m.astype(float)) for p, m in zip(pips, inpeaks)))
    corr = [_circ_cross(np.asarray(p, float), np.asarray(m, bool)) for p, m in zip(pips, inpeaks)]
    null = np.zeros(n_perm)
    for c in corr:
        null += c[rng.integers(0, len(c), n_perm)]
    mean = float(null.mean())
    p = float((1 + np.sum(null >= obs - 1e-12)) / (n_perm + 1))
    return {"n_loci": len(pips), "observed_pip_mass": obs, "null_mean": mean,
            "ratio": obs / mean if mean > 0 else float("nan"), "p_one_sided": p}
