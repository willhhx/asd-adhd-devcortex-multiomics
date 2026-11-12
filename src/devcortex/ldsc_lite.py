"""A deliberately small LD-score-regression implementation, for VALIDATION ONLY.

The real analysis uses the published LDSC software (scripts/06_run_sldsc.sh). This module
exists so the post-processing code (tau*, paired jackknife contrasts, FDR, figures) can be
checked end to end without downloading 1000 Genomes or GWAS files:

* it simulates GWAS chi-square statistics from a known tau vector;
* it fits weighted LD-score regression with a block jackknife and writes files in the same
  layout as ``ldsc.py --h2 --print-coefficients --print-delete-vals`` (``.results`` and
  ``.delete``), which the production parser then reads.

Simplifications (all stated so nobody mistakes this for LDSC): SNPs are simulated as
independent given their LD scores (no LD-induced correlation of z-scores), the regression
weights are one-step, and the intercept is free. It validates *our* code paths and the
calibration of the paired jackknife contrast, not the biology.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_ld_scores(annot: np.ndarray, chrom: np.ndarray, window: int = 20,
                       rho=0.8) -> np.ndarray:
    """Toy LD scores l_jc = sum_k r2(j,k) a_kc inside a chromosome.

    ``rho`` is a scalar or a per-SNP vector in (0, 1); r2(j, k) = (rho_j * rho_k)^(|j-k|/2).
    A per-SNP vector that varies across regions gives LD scores that vary between SNPs, as
    in real data; with a constant rho the intercept and the all-SNP annotation are nearly
    collinear and the intercept is not identifiable.
    """
    n, k = annot.shape
    rho_v = np.full(n, float(rho)) if np.isscalar(rho) else np.asarray(rho, dtype=float)
    ld = annot.astype(float).copy()
    for d in range(1, window + 1):
        w = (rho_v[:-d] * rho_v[d:]) ** (d / 2.0)          # r2 between j and j+d
        same = (chrom[d:] == chrom[:-d]).astype(float) * w
        ld[:-d] += same[:, None] * annot[d:]               # neighbour at +d contributes to j
        ld[d:] += same[:, None] * annot[:-d]               # neighbour at -d contributes to j
    return ld


def region_rho(n: int, rng: np.random.Generator, segment: int = 100,
               low: float = 0.2, high: float = 0.95) -> np.ndarray:
    """Per-SNP LD strength that is constant inside segments of ``segment`` SNPs."""
    n_seg = -(-n // segment)
    return np.repeat(rng.uniform(low, high, n_seg), segment)[:n]


def simulate_z(ld: np.ndarray, tau: np.ndarray, n_samples: float, intercept: float,
               u: np.ndarray) -> np.ndarray:
    """z_j = sqrt(E[chi2_j]) * u_j,  E[chi2_j] = intercept + N * sum_c tau_c l_jc."""
    expected = intercept + n_samples * (ld @ tau)
    return np.sqrt(np.clip(expected, 1e-6, None)) * u


def correlated_noise(n: int, rho: float, rng: np.random.Generator):
    """Two N(0,1) vectors with correlation ``rho`` (shared controls / sample overlap)."""
    a = rng.standard_normal(n)
    b = rho * a + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
    return a, b


def fit_ldsc_lite(chisq: np.ndarray, ld: np.ndarray, n_samples: float, n_blocks: int = 200):
    """Weighted LD-score regression with block jackknife.

    Returns (tau (K,), intercept, delete (n_blocks, K)).
    """
    n, k = ld.shape
    x = np.column_stack([n_samples * ld, np.ones(n)])
    # step 1: unweighted fit to get fitted values for the heteroskedasticity weights
    beta0, *_ = np.linalg.lstsq(x, chisq, rcond=None)
    fitted = np.clip(x @ beta0, 1.0, None)
    w = 1.0 / fitted ** 2
    sw = np.sqrt(w)
    xw, yw = x * sw[:, None], chisq * sw
    edges = np.floor(np.linspace(0, n, n_blocks + 1)).astype(int)
    xtx_total = xw.T @ xw
    xty_total = xw.T @ yw
    beta = np.linalg.solve(xtx_total, xty_total)
    delete = np.empty((n_blocks, k))
    for b in range(n_blocks):
        lo, hi = edges[b], edges[b + 1]
        xb, yb = xw[lo:hi], yw[lo:hi]
        beta_b = np.linalg.solve(xtx_total - xb.T @ xb, xty_total - xb.T @ yb)
        delete[b] = beta_b[:k]
    return beta[:k], float(beta[k]), delete


def results_table(names, tau, annot: np.ndarray) -> pd.DataFrame:
    """LDSC-style .results table (with --overlap-annot semantics for Prop._h2)."""
    n, k = annot.shape
    m_c = annot.sum(axis=0).astype(float)
    overlap = annot.T.astype(float) @ annot.astype(float)      # O[c, c'] = sum_j a_jc a_jc'
    h2_c = overlap @ tau                                       # h2 attributable to SNPs in c
    h2_total = float((m_c * tau).sum())
    prop_h2 = h2_c / h2_total
    prop_snps = m_c / n
    se_placeholder = np.full(k, np.nan)
    return pd.DataFrame({
        "Category": names,
        "Prop._SNPs": prop_snps,
        "Prop._h2": prop_h2,
        "Prop._h2_std_error": se_placeholder,
        "Enrichment": prop_h2 / prop_snps,
        "Enrichment_std_error": se_placeholder,
        "Enrichment_p": se_placeholder,
        "Coefficient": tau,
        "Coefficient_std_error": se_placeholder,
        "Coefficient_z-score": se_placeholder,
    })


def write_ldsc_like(prefix: str, names, tau, delete, annot: np.ndarray) -> None:
    """Write ``{prefix}.results`` and ``{prefix}.delete`` in LDSC layout."""
    res = results_table(names, tau, annot)
    from .jackknife import jackknife_se
    se = jackknife_se(delete)
    res["Coefficient_std_error"] = se
    res["Coefficient_z-score"] = tau / se
    res.to_csv(prefix + ".results", sep="\t", index=False)
    np.savetxt(prefix + ".delete", delete, delimiter="\t", fmt="%.10e")
