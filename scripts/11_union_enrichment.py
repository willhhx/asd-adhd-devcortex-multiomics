"""Union-of-peaks (pan-cortical open chromatin) enrichment per trait, from EXISTING regressions.

Each regression is baselineLD + union-of-peaks + one state, so the union annotation is the
second-last row of .results. It is located by name ('union', case-insensitive) and falls back
to index -2. Output: one row per trait (median over states of tau*, SE, enrichment) and a
per-regression table. No new LDSC runs are needed."""
import argparse, os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from devcortex import jackknife as jk  # noqa: E402
from scipy import stats  # noqa: E402


def union_index(res_path):
    cats = pd.read_csv(res_path, sep="\t")["Category"].astype(str).tolist()
    hit = [i for i, c in enumerate(cats) if "union" in c.lower()]
    return hit[0] if len(hit) == 1 else -2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sldsc-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for tdir in sorted(glob.glob(os.path.join(a.sldsc_dir, "*"))):
        if not os.path.isdir(tdir):
            continue
        trait = os.path.basename(tdir)
        for res in sorted(glob.glob(os.path.join(tdir, "*.results"))):
            dl = res[:-len(".results")] + ".part_delete"
            if not os.path.exists(dl):
                continue
            s, _ = jk.summarise_one_regression(res, dl, annot_index=union_index(res))
            rows.append({"trait": trait, "file": os.path.basename(res), **s})
    if not rows:
        sys.exit("no results found")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "union_enrichment_per_regression.tsv"), sep="\t", index=False)
    g = df.groupby("trait").agg(n=("tau_star", "size"), category=("category", "first"),
        prop_snps=("prop_snps", "median"), enrichment=("enrichment", "median"),
        enrichment_se=("enrichment_se", "median"), tau_star=("tau_star", "median"),
        tau_star_se=("tau_star_se", "median"), tau_star_min=("tau_star", "min"),
        tau_star_max=("tau_star", "max")).reset_index()
    g["z"] = g["tau_star"] / g["tau_star_se"]
    g["p_one_sided"] = stats.norm.sf(g["z"])
    g.to_csv(os.path.join(a.out, "union_enrichment_by_trait.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(g.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
