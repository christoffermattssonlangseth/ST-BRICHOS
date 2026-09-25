#!/usr/bin/env python3
"""
10c - Follow-ups to the microglia inflammation-state result (notebooks 10/10b).

(A) DECONVOLUTION-WEIGHTED signatures: instead of hard hotspot thresholding, weight
    every spot's signature score by its microglia proportion (deconv fraction) so
    microglia-rich spots dominate continuously -> closer to "pure" microglia state.
    Section-level weighted means, PBS vs BRICHOS (Mann-Whitney; WT descriptive).

(B) REGION-STRATIFIED test: does the BRICHOS state shift vary by region?
    Mixed model  score ~ treatment * region + (1|section)  on hotspot spots
    (Wald test on the treatment:region interaction), plus per-region BRICHOS-PBS
    deltas from section-level pseudobulk.

Runs REMOTELY (sc_py312) on the remote analysis machine.
"""
import warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore"); sc.settings.verbosity = 0

BASE = Path("/Users/christoffer/work/karolinska/development/ST_BRICHOS")
EXPR, ANNO = BASE / "combined_samples.h5ad", BASE / "bundle_annotations.csv"
OUT = BASE / "inflam_out"; FIG, TAB = OUT / "figs", OUT / "tables"
SAMPLE_KEY, TREAT_KEY, MG_PROP, REG = "sample_id", "treatment", "prop_MYE:Microglia", "re_annotation_regions"
ORDER = ["WT", "PBS", "BRICHOS"]; TCOL = {"WT":"#7570b3","PBS":"#d95f02","BRICHOS":"#1b9e77"}
Z_THRESH = 1.5
SIGS = {
    "Homeostatic":["P2ry12","Tmem119","Cx3cr1","Hexb","Siglech","Selplg","Csf1r","Sall1","Gpr34","Olfml3","Serinc3","Cst3"],
    "DAM":["Trem2","Tyrobp","Cst7","Itgax","Apoe","Cd68","Lpl","Ctsb","Ctsd","Cd9","Spp1","Gpnmb","Lgals3","Clec7a","B2m","Fth1","Axl"],
}

print("[1] assemble + score")
a = sc.read_h5ad(EXPR)
anno = pd.read_csv(ANNO, index_col=0)
common = a.obs_names.intersection(anno.index); a = a[common].copy()
for c in anno.columns: a.obs[c] = anno.loc[a.obs_names, c].values
a = a[a.obs[TREAT_KEY].isin(ORDER)].copy(); a.obs[SAMPLE_KEY] = a.obs[SAMPLE_KEY].astype(str)
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
for name, genes in SIGS.items():
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=f"sig_{name}", ctrl_size=max(50, len(g)))
a.obs["sig_DAMratio"] = a.obs["sig_DAM"] - a.obs["sig_Homeostatic"]
mg = a.obs[MG_PROP].astype(float).clip(lower=0)
a.obs["mg_frac"] = mg.values
a.obs["mg_z"] = ((mg - mg.mean()) / mg.std()).values
o = a.obs.copy()
METRICS = ["sig_Homeostatic", "sig_DAM", "sig_DAMratio"]

# =====================================================================
# (A) deconvolution-weighted signatures  (all spots, weight = mg_frac)
# =====================================================================
print("[2] (A) deconv-weighted section means")
def wmean(df, col, wcol="mg_frac"):
    w = df[wcol].values; return np.nan if w.sum() == 0 else float(np.average(df[col].values, weights=w))
recs = []
for (s, t), df in o.groupby([SAMPLE_KEY, TREAT_KEY], observed=True):
    recs.append(dict(sample_id=s, treatment=t, **{m: wmean(df, m) for m in METRICS}))
wsec = pd.DataFrame(recs)
wsec.to_csv(TAB / "weighted_section_means.csv", index=False)
wrows = []
for m in METRICS:
    x = wsec[wsec[TREAT_KEY] == "PBS"][m].values; z = wsec[wsec[TREAT_KEY] == "BRICHOS"][m].values
    wt = wsec[wsec[TREAT_KEY] == "WT"][m].values
    try: _, p = mannwhitneyu(x, z, alternative="two-sided")
    except Exception: p = np.nan
    wrows.append(dict(metric=m, PBS_mean=float(np.mean(x)), BRICHOS_mean=float(np.mean(z)),
                      delta=float(np.mean(z)-np.mean(x)), MWU_p=float(p), WT_mean=float(np.mean(wt))))
wres = pd.DataFrame(wrows); wres.to_csv(TAB / "weighted_summary.csv", index=False)
print(wres.to_string(index=False))

# =====================================================================
# (B) region-stratified interaction + per-region deltas (hotspots)
# =====================================================================
print("[3] (B) region-stratified test on hotspots")
h = o[(o["mg_z"] > Z_THRESH) & o[REG].notna()].copy()
h = h[h[TREAT_KEY].isin(["PBS", "BRICHOS"])].copy()
h[TREAT_KEY] = pd.Categorical(h[TREAT_KEY], categories=["PBS", "BRICHOS"])
# keep regions with both treatments present and enough spots
reg_ok = [r for r, d in h.groupby(REG, observed=True)
          if d[TREAT_KEY].nunique() == 2 and (d[TREAT_KEY] == "PBS").sum() >= 20 and (d[TREAT_KEY] == "BRICHOS").sum() >= 20]
h = h[h[REG].isin(reg_ok)].copy()
h[REG] = h[REG].astype("category")
print("    regions tested:", reg_ok)

interaction_p = {}
for m in METRICS:
    try:
        full = smf.mixedlm(f"{m} ~ C({TREAT_KEY}) * C({REG})", h, groups=h[SAMPLE_KEY]).fit(reml=False)
        red  = smf.mixedlm(f"{m} ~ C({TREAT_KEY}) + C({REG})", h, groups=h[SAMPLE_KEY]).fit(reml=False)
        lr = 2 * (full.llf - red.llf); dfd = int(full.df_modelwc - red.df_modelwc)
        from scipy.stats import chi2
        interaction_p[m] = float(chi2.sf(lr, max(dfd, 1)))
    except Exception as e:
        interaction_p[m] = np.nan; print("    interaction fit failed", m, e)
print("    treatment x region interaction p (LR):", {k: round(v, 4) for k, v in interaction_p.items()})

# per-region section-level deltas (BRICHOS - PBS)
sec = (h.groupby([REG, SAMPLE_KEY, TREAT_KEY], observed=True)[METRICS].mean().reset_index())
rrows = []
for r in reg_ok:
    for m in METRICS:
        sr = sec[sec[REG] == r]
        x = sr[sr[TREAT_KEY] == "PBS"][m].values; z = sr[sr[TREAT_KEY] == "BRICHOS"][m].values
        try: _, p = mannwhitneyu(x, z, alternative="two-sided")
        except Exception: p = np.nan
        rrows.append(dict(region=r, metric=m, PBS_mean=float(np.nanmean(x)), BRICHOS_mean=float(np.nanmean(z)),
                          delta=float(np.nanmean(z)-np.nanmean(x)), n_PBS=len(x), n_BRICHOS=len(z), MWU_p=float(p)))
rres = pd.DataFrame(rrows); rres.to_csv(TAB / "region_stratified.csv", index=False)
print(rres[rres.metric == "sig_DAMratio"].to_string(index=False))

# =====================================================================
# figure
# =====================================================================
print("[4] figure")
fig, axes = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios":[1.1, 1.6]})
# panel A: weighted section means
ax = axes[0]
for i, m in enumerate(METRICS):
    for t in ORDER:
        pts = wsec[wsec[TREAT_KEY] == t][m].values
        ax.scatter(np.full(len(pts), i) + {"WT":-0.22,"PBS":0,"BRICHOS":0.22}[t], pts,
                   color=TCOL[t], edgecolor="k", s=45, label=t if i == 0 else None)
ax.set_xticks(range(len(METRICS))); ax.set_xticklabels([m.replace("sig_","") for m in METRICS])
for i, m in enumerate(METRICS):
    p = wres[wres.metric == m]["MWU_p"].iloc[0]
    ax.annotate(f"MWU p={p:.3g}", (i, ax.get_ylim()[1]), fontsize=8, ha="center", va="bottom")
ax.axhline(0, color="grey", lw=.5, ls=":"); ax.legend(fontsize=8)
ax.set_ylabel("deconv-weighted section mean"); ax.set_title("(A) Microglia-weighted signatures (all spots)")

# panel B: per-region DAMratio delta (BRICHOS - PBS)
ax = axes[1]
d = rres[rres.metric == "sig_DAMratio"].sort_values("delta")
ypos = np.arange(len(d)); cols = ["#1b9e77" if v < 0 else "#d95f02" for v in d["delta"]]
ax.barh(ypos, d["delta"], color=cols)
ax.set_yticks(ypos); ax.set_yticklabels(d["region"])
for y, (_, rr) in zip(ypos, d.iterrows()):
    ax.annotate(f"p={rr['MWU_p']:.2g}", (rr["delta"], y), fontsize=7,
                ha="left" if rr["delta"] >= 0 else "right", va="center")
ax.axvline(0, color="k", lw=.6)
ax.set_xlabel("DAM:homeostatic delta (BRICHOS - PBS);  <0 = less inflammatory in BRICHOS")
ax.set_title(f"(B) Region-stratified state shift\ninteraction p(DAMratio)={interaction_p.get('sig_DAMratio', float('nan')):.3g}")
fig.suptitle("Deconv-weighted signatures + region-stratified BRICHOS effect", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(FIG / "inflam_weighted_regional.png", dpi=140, bbox_inches="tight")
print("[done]")
