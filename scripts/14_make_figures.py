#!/usr/bin/env python3
"""Step 14: Figures 2-6 and Table 2 from step 9/12/13 outputs (Fig 1 atlas: scripts/14a_fig1_atlas.py).

Fig 2 stage x cell-type tau* heatmaps (ASD, ADHD, controls); Fig 3 paired ASD-ADHD difference map;
Fig 4 developmental trajectories with 95% CI; Fig 5 variant-to-gene (Table 2 candidates);
Fig 6 shared vs distinct scatter. Writes PDF+PNG and the plotted numbers as TSV (source data).
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OKABE = {"ASD": "#D55E00", "ADHD": "#0072B2", "SCZ": "#009E73", "HEIGHT": "#999999"}


def grid(df, value, stages, types):
    m = pd.DataFrame(np.nan, index=types, columns=stages)
    for r in df.itertuples():
        if r.cell_type in m.index and r.stage in m.columns:
            m.loc[r.cell_type, r.stage] = getattr(r, value)
    return m


def heat(ax, m, title, cmap, vmin=None, vmax=None, stars=None):
    im = ax.imshow(m.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(m.shape[1]), [s.replace("_", " ") for s in m.columns], rotation=35, ha="right")
    ax.set_yticks(range(m.shape[0]), m.index)
    ax.set_title(title, fontsize=10)
    if stars is not None:
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                if stars.iloc[i, j]:
                    ax.text(j, i, "*", ha="center", va="center", fontsize=14)
    return im


def save(fig, out, name):
    fig.savefig(os.path.join(out, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(out, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--finemap-dir")
    ap.add_argument("--coloc-dir")
    ap.add_argument("--stage-order", required=True, help="comma-separated, developmental order")
    ap.add_argument("--targets", default="ASD,ADHD")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stages = a.stage_order.split(",")
    tau = pd.read_csv(f"{a.results}/tau_star_by_state.tsv", sep="\t")
    pc = pd.read_csv(f"{a.results}/paired_contrasts_by_state.tsv", sep="\t")
    types = sorted(tau["cell_type"].unique())
    traits = list(tau["trait"].unique())
    t1, t2 = a.targets.split(",")

    # ---- Fig 2
    fig, axs = plt.subplots(1, len(traits), figsize=(3.2 * len(traits), 3.2), squeeze=False)
    vmax = float(np.nanpercentile(tau["tau_star"], 99))
    for ax, tr in zip(axs[0], traits):
        d = tau[tau.trait == tr]
        im = heat(ax, grid(d, "tau_star", stages, types), tr, "viridis", 0, vmax,
                  grid(d, "fdr_within_trait", stages, types) < 0.05)
    fig.colorbar(im, ax=axs.ravel().tolist(), label="tau* (FDR<0.05 *)", shrink=0.8)
    save(fig, a.out, "fig2_heritability_enrichment")

    # ---- Fig 3
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    lim = float(np.nanmax(np.abs(pc["estimate"])))
    im = heat(ax, grid(pc, "estimate", stages, types), f"{t1} - {t2} (paired tau* difference)", "RdBu_r", -lim, lim,
              grid(pc, "fdr", stages, types) < 0.05)
    fig.colorbar(im, ax=ax, label="difference in tau*")
    save(fig, a.out, "fig3_difference_map")
    pc.to_csv(f"{a.out}/fig3_source_data.tsv", sep="\t", index=False)

    # ---- Fig 4 trajectories
    fig, axs = plt.subplots(1, len(types), figsize=(2.6 * len(types), 2.8), sharey=True, squeeze=False)
    x = np.arange(len(stages))
    for ax, ct in zip(axs[0], types):
        for tr in (t1, t2):
            d = tau[(tau.trait == tr) & (tau.cell_type == ct)].set_index("stage").reindex(stages)
            ax.errorbar(x, d["tau_star"], yerr=1.96 * d["tau_star_se"], color=OKABE.get(tr, "k"),
                        marker="o", ms=3, capsize=2, label=tr, alpha=0.9)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(x, [s.split("_")[0] for s in stages], rotation=35, ha="right")
        ax.set_title(ct, fontsize=9)
    axs[0][0].set_ylabel("tau* (95% CI)")
    axs[0][0].legend(frameon=False, fontsize=7)
    save(fig, a.out, "fig4_trajectories")

    # ---- Fig 5 variant-to-gene
    if a.finemap_dir:
        f = os.path.join(a.finemap_dir, f"{t2}_table2_candidates.tsv")
        if os.path.exists(f) and os.path.getsize(f) > 5:
            c = pd.read_csv(f, sep="\t").sort_values("pip_functional", ascending=False).head(25)
            fig, ax = plt.subplots(figsize=(5.5, 0.28 * len(c) + 1))
            lab = c["snp"] + "  ->  " + c["gene"].fillna("").astype(str) + " (" + c["link_type"] + ")"
            ax.barh(range(len(c)), c["pip_functional"], color="#0072B2")
            ax.barh(range(len(c)), c["pip_flat"], color="none", edgecolor="k", lw=0.6)
            ax.set_yticks(range(len(c)), lab, fontsize=6)
            ax.invert_yaxis()
            ax.set_xlabel("PIP (filled: with cell-state prior; outline: flat prior)")
            save(fig, a.out, "fig5_variant_to_gene")
            c.to_csv(f"{a.out}/table2_candidates_top.tsv", sep="\t", index=False)

    # ---- Fig 6 shared vs distinct
    w = tau.pivot_table(index="state", columns="trait", values="tau_star")
    se = tau.pivot_table(index="state", columns="trait", values="tau_star_se")
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    ax.errorbar(w[t1], w[t2], xerr=1.96 * se[t1], yerr=1.96 * se[t2], fmt="o", ms=3, lw=0.5, color="0.4")
    top = (w[t1] - w[t2]).abs().sort_values(ascending=False).head(4).index
    for s in top:
        ax.annotate(s.replace("_trimester", ""), (w.loc[s, t1], w.loc[s, t2]), fontsize=6)
    m = float(np.nanmax(w[[t1, t2]].to_numpy()))
    ax.plot([0, m], [0, m], "k--", lw=0.6)
    ax.set_xlabel(f"tau* {t1}")
    ax.set_ylabel(f"tau* {t2}")
    save(fig, a.out, "fig6_shared_vs_distinct")
    print("wrote figures to", a.out)


if __name__ == "__main__":
    main()
