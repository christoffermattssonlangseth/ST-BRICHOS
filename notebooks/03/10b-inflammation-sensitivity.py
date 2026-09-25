#!/usr/bin/env python3
"""
10b - Sensitivity checks for the microglia inflammation-state result (notebook 10).

Confirms the BRICHOS-vs-PBS state shift (homeostatic UP / DAM:homeostatic DOWN
in microglia-enriched spots) is robust to:
  (1) covariate adjustment  - add microglia fraction as a fixed covariate, so the
      state effect is not just tracking residual abundance differences;
  (2) hotspot threshold     - rerun at z > 1.0 / 1.5 / 2.0.

Runs REMOTELY (sc_py312) on the remote analysis machine, same data assembly as 10.
"""
import warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore"); sc.settings.verbosity = 0

BASE = Path("/Users/christoffer/work/karolinska/development/ST_BRICHOS")
EXPR, ANNO = BASE / "combined_samples.h5ad", BASE / "bundle_annotations.csv"
OUT = BASE / "inflam_out"; FIG, TAB = OUT / "figs", OUT / "tables"
SAMPLE_KEY, TREAT_KEY, MG_PROP = "sample_id", "treatment", "prop_MYE:Microglia"
ORDER = ["WT", "PBS", "BRICHOS"]
SIGS = {
    "Homeostatic":      ["P2ry12","Tmem119","Cx3cr1","Hexb","Siglech","Selplg",
                          "Csf1r","Sall1","Gpr34","Olfml3","Serinc3","Cst3"],
    "DAM":              ["Trem2","Tyrobp","Cst7","Itgax","Apoe","Cd68","Lpl","Ctsb",
                          "Ctsd","Cd9","Spp1","Gpnmb","Lgals3","Clec7a","B2m","Fth1","Axl"],
}

print("[1] assemble + score")
a = sc.read_h5ad(EXPR)
anno = pd.read_csv(ANNO, index_col=0)
common = a.obs_names.intersection(anno.index); a = a[common].copy()
for c in anno.columns: a.obs[c] = anno.loc[a.obs_names, c].values
a = a[a.obs[TREAT_KEY].isin(ORDER)].copy()
a.obs[SAMPLE_KEY] = a.obs[SAMPLE_KEY].astype(str)
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
for name, genes in SIGS.items():
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=f"sig_{name}", ctrl_size=max(50, len(g)))
a.obs["sig_DAMratio"] = a.obs["sig_DAM"] - a.obs["sig_Homeostatic"]
mg = a.obs[MG_PROP].astype(float); a.obs["mg_z"] = (mg - mg.mean()) / mg.std()
a.obs["mg_frac"] = mg.values

def mm(df, y, adjust):
    d = df[df[TREAT_KEY].isin(["PBS", "BRICHOS"])].copy()
    d[TREAT_KEY] = pd.Categorical(d[TREAT_KEY], categories=["PBS", "BRICHOS"])
    f = f"{y} ~ {TREAT_KEY}" + (" + mg_frac" if adjust else "")
    try:
        m = smf.mixedlm(f, d, groups=d[SAMPLE_KEY]).fit(reml=False)
        t = [x for x in m.params.index if "BRICHOS" in x][0]
        ci = m.conf_int().loc[t]
        return dict(beta=float(m.params[t]), se=float(m.bse[t]), p=float(m.pvalues[t]),
                    lo=float(ci[0]), hi=float(ci[1]), n=int(len(d)))
    except Exception as e:
        return dict(beta=np.nan, se=np.nan, p=np.nan, lo=np.nan, hi=np.nan, n=0, err=str(e))

print("[2] threshold x adjustment grid")
rows = []
for z in (1.0, 1.5, 2.0):
    h = a.obs[a.obs["mg_z"] > z]
    npb = int((h[TREAT_KEY] == "PBS").sum()); nbr = int((h[TREAT_KEY] == "BRICHOS").sum())
    for y in ("sig_Homeostatic", "sig_DAM", "sig_DAMratio"):
        for adj in (False, True):
            r = mm(h, y, adj)
            rows.append(dict(z=z, metric=y, adjusted=adj, n_PBS=npb, n_BRICHOS=nbr, **r))
res = pd.DataFrame(rows)
res.to_csv(TAB / "sensitivity_threshold_adjust.csv", index=False)
pd.set_option("display.width", 200)
print(res[["z","metric","adjusted","n_PBS","n_BRICHOS","beta","p","lo","hi"]].to_string(index=False))

print("[3] figure")
metrics = ["sig_Homeostatic", "sig_DAM", "sig_DAMratio"]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
zs = [1.0, 1.5, 2.0]; off = {False: -0.08, True: 0.08}
col = {False: "#4477AA", True: "#EE6677"}
for ax, y in zip(axes, metrics):
    for adj in (False, True):
        sub = res[(res.metric == y) & (res.adjusted == adj)].set_index("z").loc[zs]
        xpos = np.arange(len(zs)) + off[adj]
        ax.errorbar(xpos, sub["beta"], yerr=[sub["beta"]-sub["lo"], sub["hi"]-sub["beta"]],
                    fmt="o", color=col[adj], capsize=3,
                    label="adj. for mg fraction" if adj else "unadjusted")
        for xp, (_, rr) in zip(xpos, sub.iterrows()):
            ax.annotate(f"p={rr['p']:.3g}", (xp, rr["hi"]), fontsize=7, ha="center",
                        va="bottom", color=col[adj])
    ax.axhline(0, color="grey", lw=.6, ls=":")
    ax.set_xticks(range(len(zs))); ax.set_xticklabels([f"z>{z}" for z in zs])
    ax.set_title(y.replace("sig_", "")); ax.set_xlabel("hotspot threshold")
axes[0].set_ylabel("BRICHOS - PBS  (mixed-model beta, 95% CI)")
axes[0].legend(fontsize=8, loc="best")
fig.suptitle("Sensitivity: state effect vs threshold and abundance adjustment", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(FIG / "inflam_sensitivity.png", dpi=140, bbox_inches="tight")
print("[done]")
