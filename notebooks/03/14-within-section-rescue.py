#!/usr/bin/env python3
"""
14 - Within-section (difference-in-differences) version of the regional rescue tests (03/02, 03/05).

Runs LOCALLY on the (section, region) pseudobulk from 03/01.

Why
---
Treatment is confounded with slide (all PBS on P28052, 2/3 BRICHOS on P24215). In 03/02 the
disease axis is PBS(P28052) - WT(P24215) and the rescue axis is BRICHOS(mostly P24215) - PBS(P28052),
so a pure slide effect enters both with opposite sign and by itself produces "concordance" and a
negative "rescue slope". 03/02 cannot distinguish that from rescue.

Within-section design
---------------------
For every section s, region r, gene g:   y[s,r,g] = log2CPM[s,r,g] - mean_r' log2CPM[s,r',g]
(centring over the 7 regions every section has). Any section/slide/depth effect that is the same
in all regions of a section (additive on log scale) cancels. The cost: the GLOBAL rescue (same in
every region) cancels too. What survives is the regional PROFILE of rescue: does BRICHOS rescue
region r MORE than the section-average region? Centred rescue LFCs sum to ~0 over regions by
construction, so expect some regions above and some below zero.

Assumption: slide effect does not differ by region (no slide x region interaction).
Checked in [4] using BRICHOS sections on both slides (P24215_201/301 vs P28052_401).

Disease axes
------------
  PIG       : Chen 2020 literature set (slide-free; primary)
  data_up/dn: 03/02 data-driven PBS-vs-WT genes (slide-contaminated selection; secondary)

Tests per region (PBS n=6 vs BRICHOS n=3)
  concordance : fraction of disease genes moved in the rescue direction
  slope       : rescue LFC on disease LFC over disease-up genes
  perm p      : 84 reassignments of 3 BRICHOS among the 9 AD sections (gene-level binomial
                p-values ignore gene correlation and are not reported as evidence)
  within-slide: P28052_401 vs the 6 PBS sections on P28052, uncentred (no assumption; n=1)
"""
import sys, warnings
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from utils.gene_sets import PIG_CHEN2020

TBL = ROOT / "results" / "tables" / "attenuation"
OUT = ROOT / "results" / "within_section_rescue"; OUT.mkdir(parents=True, exist_ok=True)
MIN_COUNT, PSEUDO, N_BOOT, SEED = 10, 1.0, 1000, 0
LFC_T, PADJ_T = 0.5, 0.1                      # same data-driven thresholds as 03/02
TCOL = {"WT": "#7570b3", "PBS": "#d95f02", "BRICHOS": "#1b9e77"}

# ================================================================================
print("[1] load pseudobulk, keep shared regions")
counts = pd.read_csv(TBL / "pseudobulk_counts.tsv", sep="\t", index_col=0)
meta = pd.read_csv(TBL / "pseudobulk_meta.tsv", sep="\t", index_col=0)
ct = pd.crosstab(meta.region, meta["sample"])
REGIONS = sorted(ct.index[(ct > 0).all(1)])
meta = meta[meta.region.isin(REGIONS)]
counts = counts.loc[meta.index]
counts = counts.loc[:, counts.sum() >= MIN_COUNT * len(meta) / 10]
sec_meta = meta.groupby("sample").treatment.first().to_frame()
sec_meta["slide"] = sec_meta.index.str[:6]
SECTIONS = list(sec_meta.index)
print(f"    {len(SECTIONS)} sections x {len(REGIONS)} regions, {counts.shape[1]} genes")
print(f"    regions: {REGIONS}")

lib = counts.sum(1).values[:, None]
logcpm = pd.DataFrame(np.log2(counts.values / lib * 1e6 + PSEUDO), index=counts.index, columns=counts.columns)
# array [section, region, gene]
X = np.stack([np.stack([logcpm.loc[f"{s}__{r}"].values for r in REGIONS]) for s in SECTIONS])
Xc = X - X.mean(axis=1, keepdims=True)        # section-centred
genes = logcpm.columns

# ================================================================================
print("[2] disease gene sets")
dis = pd.read_csv(TBL / "disease_signature.tsv", sep="\t", index_col=0).reindex(genes)
data_up = genes[(dis.lfc > LFC_T) & (dis.padj < PADJ_T)]
data_dn = genes[(dis.lfc < -LFC_T) & (dis.padj < PADJ_T)]
pig = genes.intersection(PIG_CHEN2020)
SETS = {"PIG": (pig, pd.Index([])), "data": (data_up.union(pig), data_dn)}
for k, (u, d) in SETS.items():
    print(f"    {k}: up {len(u)}, dn {len(d)}")

def group_lfc(A, idx_a, idx_b):
    """[region, gene] mean difference between section groups."""
    return A[idx_a].mean(0) - A[idx_b].mean(0)

def rescue_stats(lfc_rg, up, dn):
    """per-region concordance + slope (on disease LFC over disease-up genes)."""
    iu, idn = genes.get_indexer(up), genes.get_indexer(dn)
    dl = dis.lfc.values
    out = []
    for j in range(len(REGIONS)):
        signed = np.concatenate([-lfc_rg[j, iu], lfc_rg[j, idn]])
        ok = np.isfinite(dl[iu])
        slope = np.polyfit(dl[iu][ok], lfc_rg[j, iu][ok], 1)[0]
        out.append((np.mean(signed > 0), np.median(signed), slope))
    return np.array(out)                       # [region, (conc, median_signed, slope)]

pbs = [i for i, s in enumerate(SECTIONS) if sec_meta.treatment[s] == "PBS"]
bri = [i for i, s in enumerate(SECTIONS) if sec_meta.treatment[s] == "BRICHOS"]
ad = pbs + bri
i401 = SECTIONS.index("P28052_401")
p28_pbs = [i for i in pbs if SECTIONS[i].startswith("P28052")]
bri24 = [i for i in bri if SECTIONS[i].startswith("P24215")]

# ================================================================================
print("[3] observed + 84-permutation null")
DESIGNS = {
    "uncentred (03/02 analogue)": (X, bri, pbs),
    "within-section (centred)":   (Xc, bri, pbs),
    "within-slide 401 vs P28052 PBS": (X, [i401], p28_pbs),
    "within-slide 401, centred":  (Xc, [i401], p28_pbs),
}
rows = []
perms = [list(c) for c in combinations(ad, len(bri))]
for set_name, (up, dn) in SETS.items():
    for dname, (A, ia, ib) in DESIGNS.items():
        obs = rescue_stats(group_lfc(A, ia, ib), up, dn)
        if dname.startswith("within-slide"):
            # only one BRICHOS on P28052: null = each of the 7 P28052 sections taking the BRICHOS role
            pool = p28_pbs + [i401]
            null = np.stack([rescue_stats(group_lfc(A, [k], [i for i in pool if i != k]), up, dn) for k in pool])
        else:
            null = np.stack([rescue_stats(group_lfc(A, p, [i for i in ad if i not in p]), up, dn) for p in perms])
        for j, r in enumerate(REGIONS):
            rows.append(dict(gene_set=set_name, design=dname, region=r,
                             concordance=obs[j, 0], median_signed_lfc=obs[j, 1], slope=obs[j, 2],
                             perm_p_conc=np.mean(null[:, j, 0] >= obs[j, 0]),
                             perm_p_slope=np.mean(null[:, j, 2] <= obs[j, 2]),
                             n_perm=len(null)))
res = pd.DataFrame(rows)
res.to_csv(OUT / "rescue_within_section.csv", index=False)
with pd.option_context("display.width", 200, "display.max_rows", 200):
    for set_name in SETS:
        print(f"\n--- gene set: {set_name}")
        print(res[res.gene_set == set_name].drop(columns=["gene_set", "n_perm"])
              .round(3).to_string(index=False))

# ================================================================================
print("\n[4] assumption check: is the slide effect region-uniform?")
# both slides carry BRICHOS: slide effect estimate = mean(P24215 BRICHOS) - P28052_401, per region x gene.
# Expressed genes only; compare the region-uniform part (gene mean over regions) with the
# region-specific part (residual). Animal-to-animal noise inflates the residual, so calibrate
# against the same split among PBS sections on one slide (no slide effect at all).
expr = logcpm.mean(0).values > 3
def split_var(D):
    D = D[:, expr]
    main = D.mean(0, keepdims=True)
    return np.var(main) , np.var(D - main)
slide_eff = X[bri24].mean(0) - X[i401]
v_main, v_inter = split_var(slide_eff)
calib = []
for k in p28_pbs:                               # 2-vs-1 splits of PBS sections, same slide
    for pair in combinations([i for i in p28_pbs if i != k], 2):
        calib.append(split_var(X[list(pair)].mean(0) - X[k]))
calib = np.array(calib)
chk = pd.DataFrame({
    "contrast": ["slide (P24215 BRI - P28052_401)", "noise (PBS 2-vs-1, same slide), median"],
    "var_region_uniform": [v_main, np.median(calib[:, 0])],
    "var_region_specific": [v_inter, np.median(calib[:, 1])],
})
chk["frac_region_uniform"] = chk.var_region_uniform / (chk.var_region_uniform + chk.var_region_specific)
pct = np.mean(calib[:, 1] >= v_inter)
print(chk.round(4).to_string(index=False))
print(f"    region-specific variance of slide contrast exceeds {100 * (1 - pct):.0f}% of same-slide noise splits")
chk.to_csv(OUT / "slide_uniformity_check.csv", index=False)

# per-region: how well does the slide contrast follow its own gene mean?
r_reg = [np.corrcoef(slide_eff[j, expr], slide_eff[:, expr].mean(0))[0, 1] for j in range(len(REGIONS))]
print("    corr(region slide effect, region-average slide effect):",
      dict(zip(REGIONS, np.round(r_reg, 2))))

# ================================================================================
print("[5] figure")
fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)
order = (res[(res.gene_set == "PIG") & (res.design == "uncentred (03/02 analogue)")]
         .sort_values("slope").region.tolist())
for ax, metric, ref, lab in [(axes[0], "concordance", 0.5, "PIG concordance (frac. moved toward rescue)"),
                             (axes[1], "slope", 0.0, "PIG rescue slope (neg = rescue)"),
                             (axes[2], "perm_p_slope", None, "perm p (slope)")]:
    for k, (dname, mk, c) in enumerate([("uncentred (03/02 analogue)", "o", "#999999"),
                                        ("within-section (centred)", "s", "#1b9e77"),
                                        ("within-slide 401 vs P28052 PBS", "^", "#d95f02")]):
        d = res[(res.gene_set == "PIG") & (res.design == dname)].set_index("region").loc[order]
        ax.scatter(d[metric], np.arange(len(order)) + (k - 1) * 0.2, marker=mk, color=c, s=55, label=dname)
    if ref is not None: ax.axvline(ref, color="grey", ls=":", lw=.8)
    else: ax.axvline(0.05, color="red", ls="--", lw=.6); ax.set_xscale("log")
    ax.set_xlabel(lab)
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels(order)
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle("Regional rescue: original (slide-confounded) vs within-section (slide-cancelled; global part removed) "
             "vs within-slide (401 only)", fontsize=11)
fig.tight_layout(); fig.savefig(OUT / "rescue_within_section.png", dpi=140, bbox_inches="tight")
print("done")
