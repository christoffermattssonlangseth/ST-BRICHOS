#!/usr/bin/env python3
"""
13 - Within-section contrasts: does BRICHOS change the microglial RESPONSE,
     measured relative to each section's own background?

Runs on the REMOTE analysis machine in env `sc_py312` (same layout as 03/10, 03/12).

Why
---
Treatment is confounded with slide (and sequencing depth; see 03/12). Slide/depth act
mostly as section-wide shifts, so any difference between two spot groups INSIDE the
same section largely cancels them. We compute one such contrast per section x region
and compare the contrasts across treatments.

Contrast
--------
Within each section x region (>= MIN_SPOTS spots, regions shared by all sections):
  hot  = top HOT_Q of microglia deconv fraction (ranked WITHIN that section x region)
  cold = bottom COLD_Q
  diff  = mean score(hot) - mean score(cold)
  slope = OLS slope of score on within-section x region microglia rank (0..1)
Section value = mean over shared regions (equal region weight).
Expression is downsampled to DOWNSAMPLE UMI/spot so noisier scores in shallower
sections can't masquerade as a weaker response.

Controls
--------
  positive : DAM / PIG contrasts should be > 0 in AD sections (plaque niches)
  WT       : no plaques -> small contrast; WT vs BRICHOS is a SAME-SLIDE comparison (P24215)
  negative : neuronal / synapse / oligo / endothelial contrasts should not differ by
             treatment; if they do, slide is leaking into the contrast

Tests (PBS n=6 vs BRICHOS n=3; WT n=1 descriptive)
  MWU on section values; LMM on section x region values (diff ~ treatment + region,
  1|section); within-slide: rank of P28052_401 among the 7 P28052 sections (min p=1/7).
"""
import warnings, sys
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

BASE = Path("/Users/christoffer/work/karolinska/development/ST_BRICHOS")
EXPR = BASE / "combined_samples.h5ad"
ANNO = BASE / "bundle_annotations.csv"
OUT  = BASE / "contrast_out"; FIG = OUT / "figs"
for d in (OUT, FIG): d.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(BASE / "classifier_utils"))
from gene_sets import (PIG_CHEN2020, DAM_KEREN_SHAUL2017, MGND_KRASEMANN2017,
                       NEURONAL_MARKERS, SYNAPSE, OLIGODENDROCYTE_MATURE, ENDOTHELIAL,
                       filter_to_var)

SAMPLE_KEY, TREAT_KEY, REGION_KEY = "sample_id", "treatment", "re_annotation_regions"
MG_PROP = "prop_MYE:Microglia"
DOWNSAMPLE, MIN_SPOTS, HOT_Q, COLD_Q, SEED = 5000, 30, 0.20, 0.50, 0
ORDER = ["WT", "PBS", "BRICHOS"]
TCOL  = {"WT": "#7570b3", "PBS": "#d95f02", "BRICHOS": "#1b9e77"}

# same microglia-state sets as 03/10, plus registry sets and negative controls
SIGS = {
    "Homeostatic":      ["P2ry12","Tmem119","Cx3cr1","Hexb","Siglech","Selplg",
                          "Csf1r","Sall1","Gpr34","Olfml3","Serinc3","Cst3"],
    "DAM":              DAM_KEREN_SHAUL2017,
    "PIG":              PIG_CHEN2020,
    "MGnD":             MGND_KRASEMANN2017,
    "Neuroinflammation":["Il1b","Tnf","Nfkbia","Ccl3","Ccl4","Ccl2","Cxcl10","Nlrp3",
                          "Il6","Tlr2","Ptgs2","Nfkb1","Rela","Il1a","Tnfaip3"],
    "Complement":       ["C1qa","C1qb","C1qc","C3","C4b","C1qbp","Cfh"],
    "ReactiveAstro":    ["Gfap","Vim","Serpina3n","C3","Lcn2","Cp","Steap4"],
    # negative controls
    "neg_Neuronal":     NEURONAL_MARKERS,
    "neg_Synapse":      SYNAPSE,
    "neg_Oligo":        OLIGODENDROCYTE_MATURE,
    "neg_Endothelial":  ENDOTHELIAL,
}

# ================================================================================
print("[1] load, join, downsample")
a = sc.read_h5ad(EXPR)
anno = pd.read_csv(ANNO, index_col=0)
a = a[a.obs_names.intersection(anno.index)].copy()
for c in [SAMPLE_KEY, TREAT_KEY, REGION_KEY, MG_PROP]:
    a.obs[c] = anno.loc[a.obs_names, c].values
a.obs[SAMPLE_KEY] = a.obs[SAMPLE_KEY].astype(str)
a = a[np.asarray(a.X.sum(1)).ravel() >= DOWNSAMPLE].copy()
sc.pp.downsample_counts(a, counts_per_cell=DOWNSAMPLE, random_state=SEED)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

ct = pd.crosstab(a.obs[REGION_KEY], a.obs[SAMPLE_KEY])
shared = list(ct.index[(ct >= MIN_SPOTS).all(1)])
a = a[a.obs[REGION_KEY].isin(shared)].copy()
print(f"    {a.n_obs} spots, shared regions: {shared}")

print("[2] scoring signatures")
for name, genes in SIGS.items():
    g = filter_to_var(genes, a.var_names)
    sc.tl.score_genes(a, g, score_name=name, ctrl_size=max(50, len(g)), random_state=SEED)
    print(f"    {name}: {len(g)}/{len(genes)}")
a.obs["DAMratio"] = a.obs["DAM"] - a.obs["Homeostatic"]
SCORES = ["DAMratio"] + list(SIGS)
o = a.obs[[SAMPLE_KEY, TREAT_KEY, REGION_KEY, MG_PROP] + SCORES].copy()
o["slide"] = o[SAMPLE_KEY].str[:6]

# ================================================================================
print("[3] within-section x region contrasts")
rows = []
for (s, r), d in o.groupby([SAMPLE_KEY, REGION_KEY], observed=True):
    rank = d[MG_PROP].rank(pct=True).values
    hot, cold = rank > 1 - HOT_Q, rank <= COLD_Q
    row = dict(section=s, region=r, treatment=d[TREAT_KEY].iloc[0], slide=s[:6], n=len(d),
               mg_hot=d[MG_PROP].values[hot].mean(), mg_cold=d[MG_PROP].values[cold].mean())
    for y in SCORES:
        v = d[y].values
        row[f"diff_{y}"] = v[hot].mean() - v[cold].mean()
        row[f"slope_{y}"] = np.polyfit(rank, v, 1)[0]
    rows.append(row)
sr = pd.DataFrame(rows)
sr.to_csv(OUT / "contrasts_section_region.csv", index=False)
sec = sr.groupby(["section", "treatment", "slide"]).mean(numeric_only=True).reset_index()
sec.to_csv(OUT / "contrasts_section.csv", index=False)

# ================================================================================
print("[4] stats")
def lmm(y):
    d = sr[sr.treatment.isin(["PBS", "BRICHOS"])].copy()
    d["treatment"] = pd.Categorical(d.treatment, ["PBS", "BRICHOS"])
    try:
        m = smf.mixedlm(f"{y} ~ treatment + region", d, groups=d.section).fit(reml=False)
        t = [k for k in m.params.index if "BRICHOS" in k][0]
        return m.params[t], m.pvalues[t]
    except Exception:
        return np.nan, np.nan

stats = []
for kind in ("diff", "slope"):
    for y in SCORES:
        col = f"{kind}_{y}"
        v = sec.set_index("section")[col]
        pbs = sec[sec.treatment == "PBS"][col]; bri = sec[sec.treatment == "BRICHOS"][col]
        beta, p_lmm = lmm(col)
        p28 = sec[sec.slide == "P28052"].set_index("section")[col]
        # 401 rank among P28052 sections, in the direction of the BRICHOS-vs-PBS effect
        sign = np.sign(bri.mean() - pbs.mean()) or 1
        r401 = int((sign * p28).rank(ascending=False)["P28052_401"])
        bri24 = sec[(sec.treatment == "BRICHOS") & (sec.slide == "P24215")][col].mean()
        stats.append(dict(kind=kind, score=y, negative_control=y.startswith("neg_"),
                          WT=v.get("P24215_101", np.nan), PBS_mean=pbs.mean(), BRICHOS_mean=bri.mean(),
                          delta=bri.mean() - pbs.mean(), MWU_p=mannwhitneyu(pbs, bri).pvalue,
                          lmm_beta=beta, lmm_p=p_lmm,
                          P28052_401=v["P28052_401"], P28052_PBS_min=p28.drop("P28052_401").min(),
                          P28052_PBS_max=p28.drop("P28052_401").max(),
                          rank401_in_P28052=r401, p401=r401 / len(p28),
                          P24215_BRI_minus_WT=bri24 - v.get("P24215_101", np.nan)))
st = pd.DataFrame(stats)
st.to_csv(OUT / "contrast_stats.csv", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 30):
    print(st[st.kind == "diff"].drop(columns=["kind"]).round(3).to_string(index=False))
    print(st[st.kind == "slope"][["score", "delta", "MWU_p", "lmm_p", "rank401_in_P28052"]].round(3).to_string(index=False))

# ================================================================================
print("[5] figure")
show = ["DAMratio", "Homeostatic", "DAM", "PIG", "MGnD", "Complement", "Neuroinflammation",
        "ReactiveAstro", "neg_Neuronal", "neg_Synapse", "neg_Oligo", "neg_Endothelial"]
fig, axes = plt.subplots(2, 6, figsize=(22, 8), sharex=True)
rng = np.random.default_rng(SEED)
for ax, y in zip(axes.ravel(), show):
    col = f"diff_{y}"
    for i, t in enumerate(ORDER):
        d = sec[sec.treatment == t]
        x = i + rng.uniform(-.12, .12, len(d))
        for xi, (_, row) in zip(x, d.iterrows()):
            ax.scatter(xi, row[col], s=70, color=TCOL[t],
                       marker="o" if row.slide == "P28052" else "^",
                       edgecolor="k" if row.section == "P28052_401" else "none", linewidth=1.8, zorder=3)
        if t != "WT": ax.hlines(d[col].mean(), i - .25, i + .25, color=TCOL[t], lw=2)
    ax.axhline(0, color="grey", lw=.6, ls=":")
    s_ = st[(st.kind == "diff") & (st.score == y)].iloc[0]
    ax.set_title(f"{y}\nMWU p={s_.MWU_p:.2f}  LMM p={s_.lmm_p:.2f}  401 rank {s_.rank401_in_P28052}/7",
                 fontsize=9, color="dimgrey" if y.startswith("neg_") else "k")
    ax.set_xticks(range(3)); ax.set_xticklabels(ORDER)
axes[0, 0].set_ylabel("score(microglia-hot) - score(cold)\nwithin section x region")
axes[1, 0].set_ylabel("same")
fig.suptitle("Within-section microglia-niche contrasts (depth-matched). "
             "o = slide P28052, ^ = P24215, black edge = P28052_401 (only BRICHOS on the PBS slide). "
             "Bottom-right four = negative controls.", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "within_section_contrasts.png", dpi=140, bbox_inches="tight")
print("done")
