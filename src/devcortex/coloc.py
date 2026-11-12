"""coloc.abf (Giambartolomei et al. 2014) in numpy, for summary-statistic colocalisation.

Hypotheses: H0 no association; H1 only trait 1; H2 only trait 2; H3 both, different causal
variants; H4 both, one shared causal variant. Priors follow the coloc defaults
(p1 = p2 = 1e-4, p12 = 1e-5). Inputs are per-variant BETA and SE for the two traits on the
SAME variants (merged on SNP id upstream) and in the same allele orientation.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .finemap import log_abf


def coloc_abf(beta1, se1, beta2, se2, p1: float = 1e-4, p2: float = 1e-4, p12: float = 1e-5,
              sd_prior1: float = 0.2, sd_prior2: float = 0.15) -> dict:
    l1 = log_abf(beta1, se1, sd_prior1)
    l2 = log_abf(beta2, se2, sd_prior2)
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    lh0 = 0.0
    lh1 = np.log(p1) + s1
    lh2 = np.log(p2) + s2
    # sum over pairs of DIFFERENT variants: sum_i sum_j!=i ABF1_i ABF2_j = S1*S2 - S12
    ratio = np.exp(s12 - (s1 + s2))
    # ratio can reach 1 (e.g. one variant in the region); then H3 has no admissible pairs
    diff = -np.inf if ratio >= 1.0 else s1 + s2 + np.log1p(-ratio)
    lh3 = np.log(p1) + np.log(p2) + diff
    lh4 = np.log(p12) + s12
    lall = np.array([lh0, lh1, lh2, lh3, lh4])
    pp = np.exp(lall - logsumexp(lall))
    # SNP-level posterior of being the shared causal variant under H4
    snp_pp4 = np.exp((l1 + l2) - s12)
    return {"PP.H0": pp[0], "PP.H1": pp[1], "PP.H2": pp[2], "PP.H3": pp[3], "PP.H4": pp[4],
            "snp_pp_h4": snp_pp4, "n_snps": len(l1)}
