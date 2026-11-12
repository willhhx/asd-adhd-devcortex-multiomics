"""Power bookkeeping for the ASD-versus-ADHD comparison.

The two GWAS differ ~3x in effective sample size, so a given true difference in tau* is far
easier to detect when the ADHD estimate carries most of the signal. These helpers put that
on the page (Table S, Methods) and are used to label contrasts that are underpowered.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def effective_n(n_cases: float, n_controls: float) -> float:
    """Effective sample size of a case-control GWAS: 4 / (1/Nca + 1/Nco)."""
    return 4.0 / (1.0 / n_cases + 1.0 / n_controls)


def expected_se_ratio(neff_a: float, neff_b: float) -> float:
    """Approximate SE(tau*_A) / SE(tau*_B) at equal heritability: sqrt(Neff_B / Neff_A)."""
    return float(np.sqrt(neff_b / neff_a))


def power_two_sided(delta: float, se: float, alpha: float = 0.05) -> float:
    """Power of a two-sided z-test for a true difference ``delta`` with standard error ``se``."""
    crit = stats.norm.isf(alpha / 2)
    z = delta / se
    return float(stats.norm.sf(crit - z) + stats.norm.cdf(-crit - z))


def minimum_detectable_difference(se: float, alpha: float = 0.05, power: float = 0.8) -> float:
    return float((stats.norm.isf(alpha / 2) + stats.norm.ppf(power)) * se)


def bonferroni_alpha(alpha: float, n_tests: int) -> float:
    return alpha / n_tests
