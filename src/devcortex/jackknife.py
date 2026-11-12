"""Block-jackknife post-processing for partitioned LD-score regression (S-LDSC).

LDSC (``ldsc.py --h2 --print-coefficients --print-delete-vals``) writes, for each
regression, a ``.results`` table and a ``.delete`` matrix (n_blocks x n_annotations)
holding the regression coefficients re-estimated with each genomic block removed.
Everything here works on those objects, so it is identical for real LDSC output and
for the synthetic output of :mod:`devcortex.ldsc_lite`.

Key quantities
--------------
tau*   standardised per-annotation effect (Gazal et al. 2017; Finucane et al. 2018):
       tau*_c = tau_c * sd_c / sum_c'(tau_c' * p_c'),  p_c = fraction of reference SNPs in c.
       (M_total cancels because h2 = M_total * sum_c tau_c p_c.)
contrast  a linear combination of tau* values across (trait, state) pairs. Because both
       traits are regressed on the SAME SNPs and the SAME block partition, the jackknife
       delete values can be combined block by block. The resulting SE therefore carries
       the covariance induced by shared annotation structure and shared controls
       (the iPSYCH controls are used by both the ASD and ADHD GWAS).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def jackknife_var(delete_vals) -> np.ndarray:
    """Variance of an estimator from leave-one-block-out values (axis 0 = blocks)."""
    d = np.asarray(delete_vals, dtype=float)
    n = d.shape[0]
    if n < 2:
        raise ValueError("need at least 2 jackknife blocks")
    return (n - 1.0) / n * ((d - d.mean(axis=0)) ** 2).sum(axis=0)


def jackknife_se(delete_vals) -> np.ndarray:
    return np.sqrt(jackknife_var(delete_vals))


def tau_star(tau, p_snps, annot_index, annot_sd):
    """Standardised coefficient for one annotation.

    Parameters
    ----------
    tau : array (..., K) coefficients; the trailing axis indexes annotations. Works for
          a single vector (K,) or the delete matrix (n_blocks, K).
    p_snps : array (K,) fraction of reference SNPs in each annotation (``Prop._SNPs``).
    annot_index : int, column of the annotation of interest (usually the last one).
    annot_sd : float, SD of that annotation across reference SNPs.
    """
    tau = np.asarray(tau, dtype=float)
    p = np.asarray(p_snps, dtype=float)
    h2_over_m = np.sum(tau * p, axis=-1)  # h2 / M_total
    return tau[..., annot_index] * annot_sd / h2_over_m


def binary_annotation_sd(p_snps_of_annot: float) -> float:
    p = float(p_snps_of_annot)
    return float(np.sqrt(p * (1.0 - p)))


def one_sided_p(z):
    """P(Z > z); the convention for 'is the annotation enriched' (Finucane 2018)."""
    return stats.norm.sf(z)


def contrast(estimates: dict, deletes: dict, weights: dict, ci: float = 0.95) -> dict:
    """Paired block-jackknife test of a linear contrast.

    estimates : {key: full-data estimate}
    deletes   : {key: array (n_blocks,) of leave-one-block-out estimates}
    weights   : {key: weight}; the contrast is sum_k w_k * estimate_k.

    Example (ASD minus ADHD in one cell state):
        contrast(est, dele, {("ASD", s): 1, ("ADHD", s): -1})
    Example (difference in differences between two stages):
        weights = {("ASD", s2): 1, ("ASD", s3): -1, ("ADHD", s2): -1, ("ADHD", s3): 1}
    """
    keys = list(weights)
    n_blocks = {len(np.asarray(deletes[k])) for k in keys}
    if len(n_blocks) != 1:
        raise ValueError(
            "traits were regressed on different block partitions "
            f"(n_blocks={sorted(n_blocks)}); align SNP sets before running LDSC"
        )
    est = sum(weights[k] * estimates[k] for k in keys)
    dele = sum(weights[k] * np.asarray(deletes[k], dtype=float) for k in keys)
    se = float(jackknife_se(dele))
    z = est / se if se > 0 else np.nan
    p = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else np.nan
    crit = stats.norm.ppf(0.5 + ci / 2)
    return {"estimate": float(est), "se": se, "z": float(z), "p_two_sided": float(p),
            "ci_low": float(est - crit * se), "ci_high": float(est + crit * se)}


def unpaired_difference(est_a, se_a, est_b, se_b) -> dict:
    """Naive difference that ignores covariance between the two regressions.

    Kept only as a sensitivity analysis / to show what pairing changes."""
    se = float(np.hypot(se_a, se_b))
    est = est_a - est_b
    z = est / se
    return {"estimate": est, "se": se, "z": z, "p_two_sided": float(2 * stats.norm.sf(abs(z)))}


def benjamini_hochberg(pvals) -> np.ndarray:
    """BH-adjusted p-values (q-values); NaNs are passed through."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(ranked, 1.0)
    out[ok] = q
    return out


def parse_ldsc_results(results_path: str) -> pd.DataFrame:
    """Read an LDSC ``.results`` file (tab separated, one row per annotation)."""
    df = pd.read_csv(results_path, sep="\t")
    needed = {"Category", "Prop._SNPs", "Prop._h2", "Enrichment", "Coefficient",
              "Coefficient_std_error"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{results_path}: missing columns {sorted(missing)}")
    return df


def read_delete_matrix(delete_path: str) -> np.ndarray:
    """Read an LDSC ``.delete`` file (no header, tab separated)."""
    arr = np.loadtxt(delete_path, ndmin=2)
    return arr


def summarise_one_regression(results_path: str, delete_path: str, annot_index: int = -1):
    """Return tau*, its delete values, enrichment and basic metadata for one regression.

    The annotation of interest is, by construction in scripts/06_run_sldsc.sh, the last
    reference annotation (baseline-LD + union-of-peaks + one cell state).
    """
    res = parse_ldsc_results(results_path)
    dele = read_delete_matrix(delete_path)
    if dele.shape[1] != len(res):
        raise ValueError(
            f"{delete_path} has {dele.shape[1]} columns but {results_path} has {len(res)} rows"
        )
    k = annot_index % len(res)
    p = res["Prop._SNPs"].to_numpy(float)
    sd = binary_annotation_sd(p[k])
    tau = res["Coefficient"].to_numpy(float)
    est = float(tau_star(tau, p, k, sd))
    dels = tau_star(dele, p, k, sd)
    se = float(jackknife_se(dels))
    row = res.iloc[k]
    return {
        "category": row["Category"],
        "prop_snps": float(p[k]),
        "enrichment": float(row["Enrichment"]),
        "enrichment_se": float(row.get("Enrichment_std_error", np.nan)),
        "tau_star": est,
        "tau_star_se": se,
        "tau_star_z": est / se if se > 0 else np.nan,
        "tau_star_p_one_sided": float(one_sided_p(est / se)) if se > 0 else np.nan,
        "n_blocks": int(dele.shape[0]),
    }, dels
