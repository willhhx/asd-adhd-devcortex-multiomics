import argparse, os
import pandas as pd
from scipy import stats
ap = argparse.ArgumentParser()
ap.add_argument("--f5-tau", required=True)
ap.add_argument("--sud-tau", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
f5 = pd.read_csv(a.f5_tau, sep="\t"); f5 = f5[f5["trait"] == "F5"]
sud = pd.read_csv(a.sud_tau, sep="\t")
d = pd.concat([f5, sud])
w = d.pivot(index="state", columns="trait", values="tau_star")
se = d.pivot(index="state", columns="trait", values="tau_star_se")
p = d.pivot(index="state", columns="trait", values="tau_star_p_one_sided")
out = pd.concat({"tau": w, "se": se, "p1": p}, axis=1)
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
out.to_csv(a.out, sep="\t")
pd.set_option("display.width", 200)
print(out.round(3).sort_values(("tau", "F5"), ascending=False).to_string())
for t in [c for c in w.columns if c != "F5"]:
    r, pr = stats.spearmanr(w["F5"], w[t])
    print(f"Spearman F5 vs {t}: rho={r:.2f} (p={pr:.2f}, n={len(w)}); nominal p<0.05 groups: {(p[t] < .05).sum()}")
