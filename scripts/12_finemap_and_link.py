#!/usr/bin/env python3
"""Step 12: locus fine-mapping with cell-state functional priors and variant-to-gene links.

For one trait:
  1. loci = independent genome-wide-significant (p < threshold) regions +- window;
  2. prior = variants inside the specific peaks of that trait's FDR-significant states get
     weight min(enrichment, cap) (enrichment from S-LDSC), everything else 1;
  3. Wakefield-ABF PIPs with and without the prior; 95% credible sets;
  4. credible-set variants that fall in a state-specific peak are linked to genes through
     (a) peak-gene links (SCENIC+/co-accessibility/Hi-C table) and, failing that, (b) the
     nearest TSS, and the link type is recorded. The paper reports only (a) as
     "variant-to-gene" and treats (b) as proximity.
Output: <out>/<trait>_credible_variants.tsv (all variants in credible sets, both priors) and
<out>/<trait>_table2_candidates.tsv (variants in a significant-state peak, with gene links).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import finemap, peaks  # noqa: E402
from devcortex.config import load_config  # noqa: E402


def read_gwas(path):
    df = pd.read_csv(path, sep="\t")
    need = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"}
    if not need <= set(df.columns):
        sys.exit(f"{path}: need columns {sorted(need)}; got {list(df.columns)}")
    df["CHR"] = df["CHR"].astype(str).str.replace("chr", "", regex=False)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--trait", required=True)
    ap.add_argument("--gwas", required=True, help="harmonised full summary statistics (SNP CHR BP A1 A2 BETA SE P)")
    ap.add_argument("--tau-table", required=True, help="tau_star_by_state.tsv from step 9")
    ap.add_argument("--peakset-dir", required=True, help="hg19 BED directory with states/")
    ap.add_argument("--peakset-summary", required=True)
    ap.add_argument("--peak-gene-links", help="TSV: peak, gene, source, score")
    ap.add_argument("--genes", help="TSV: gene, chrom, tss (hg19) for the nearest-TSS fallback")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fm = cfg["finemap"]
    os.makedirs(args.out, exist_ok=True)

    gwas = read_gwas(args.gwas)
    tau = pd.read_csv(args.tau_table, sep="\t")
    tau = tau[(tau["trait"] == args.trait)]
    sig = tau[tau["fdr_within_trait"] < fm["fdr_for_prior"]]
    summ = pd.read_csv(args.peakset_summary, sep="\t")
    fname = dict(zip(summ["state"], summ["file"]))
    print(f"{args.trait}: {len(sig)} FDR-significant states used for the functional prior")

    # peak id containing each SNP, per significant state
    state_hit = {}
    for state in sig["state"]:
        bed = peaks.read_bed(os.path.join(args.peakset_dir, "states", fname[state] + ".bed"))
        state_hit[state] = peaks.variant_peak_overlap(gwas, bed)
    in_any = pd.Series(False, index=gwas.index)
    for s in state_hit.values():
        in_any |= s.notna()
    enrich = float(np.clip(sig["enrichment"].mean(), 1.0, fm["enrichment_cap"])) if len(sig) else 1.0
    links = pd.read_csv(args.peak_gene_links, sep="\t") if args.peak_gene_links else pd.DataFrame(columns=["peak", "gene"])
    genes = pd.read_csv(args.genes, sep="\t") if args.genes else pd.DataFrame(columns=["gene", "chrom", "tss"])

    loci = finemap.define_loci(gwas, p_thresh=fm["p_threshold"], window=fm["window"])
    print(f"{len(loci)} independent loci at p < {fm['p_threshold']:g}; prior enrichment = {enrich:.2f}")
    cred_rows, cand_rows = [], []
    for L in loci.itertuples():
        mask = (gwas["CHR"] == str(L.chrom)) & gwas["BP"].between(L.start, L.end)
        sub = gwas[mask].copy()
        if sub.empty:
            continue
        prior = finemap.functional_prior(in_any[mask].to_numpy(), enrich, cap=fm["enrichment_cap"])
        flat = finemap.finemap_locus(sub, None, fm["credible_level"])
        fun = finemap.finemap_locus(sub, prior, fm["credible_level"])
        sub["locus"] = L.locus
        sub["pip_flat"] = flat["pip"].to_numpy()
        sub["pip_functional"] = fun["pip"].to_numpy()
        sub["in_cs_flat"] = flat["in_credible_set"].to_numpy()
        sub["in_cs_functional"] = fun["in_credible_set"].to_numpy()
        # first significant state whose peak contains the SNP
        sp_ = pd.Series(np.nan, index=sub.index, dtype=object)
        st_ = pd.Series(np.nan, index=sub.index, dtype=object)
        for state, hit in state_hit.items():
            h = hit.loc[sub.index]
            new = h.notna() & sp_.isna()
            sp_[new] = h[new]
            st_[new] = state
        sub["state_peak"], sub["state"] = sp_, st_
        sub["lead_snp"] = L.lead_snp
        cred = sub[sub["in_cs_functional"] | sub["in_cs_flat"]]
        cred_rows.append(cred)
        for r in cred[cred["state_peak"].notna() & cred["in_cs_functional"]].itertuples():
            lk = links[links["peak"] == r.state_peak]
            if len(lk):
                gene, ltype = ";".join(sorted(set(lk["gene"]))), "peak_gene_link"
            else:
                g = genes[genes["chrom"].astype(str) == str(r.CHR)]
                if len(g):
                    j = (g["tss"] - r.BP).abs().idxmin()
                    gene, ltype = g.loc[j, "gene"], "nearest_tss"
                else:
                    gene, ltype = "", "none"
            cand_rows.append({"trait": args.trait, "locus": r.locus, "lead_snp": r.lead_snp, "snp": r.SNP,
                              "chrom": r.CHR, "pos": r.BP, "p": r.P, "pip_flat": r.pip_flat,
                              "pip_functional": r.pip_functional, "state": r.state, "peak": r.state_peak,
                              "gene": gene, "link_type": ltype})
    cred_df = pd.concat(cred_rows) if cred_rows else pd.DataFrame()
    cred_df.to_csv(os.path.join(args.out, f"{args.trait}_credible_variants.tsv"), sep="\t", index=False)
    cand = pd.DataFrame(cand_rows)
    cand.to_csv(os.path.join(args.out, f"{args.trait}_table2_candidates.tsv"), sep="\t", index=False)
    loci.to_csv(os.path.join(args.out, f"{args.trait}_loci.tsv"), sep="\t", index=False)
    print(f"credible-set variants: {len(cred_df)}; in significant-state peaks: {len(cand)}")


if __name__ == "__main__":
    main()
