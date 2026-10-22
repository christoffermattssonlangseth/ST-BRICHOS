#!/usr/bin/env python3
"""
15 - Model-based slide removal, then re-run the regional rescue tests (03/02, 03/14).

Runs LOCALLY on the (section, region) pseudobulk from 03/01.

Treatment and slide are confounded, so slide can only be removed under an assumption. Three
methods, each with a different one:

  RUVs  (replicate-based, Risso 2014): unwanted factors = top singular vectors of within-group
        residuals (group = treatment x region). Only BRICHOS spans both slides, so the slide
        direction is learnt from P24215 BRICHOS vs P28052_401, pooled over genes.
        Assumes: slide acts along a few gene-wide factors.
  RUVg  (control genes): factors = top singular vectors of region-centred expression of
        "boring" genes (expressed, stable across regions and PBS sections, not PIG/disease).
        Assumes: BRICHOS does not change those genes.
  tech  (CQN-like): slide contrast among BRICHOS sections (P24215 BRI - 401) modelled as a
        function of gene GC, length, expression and mito/ribo membership (gradient boosting,
        gene-wise CV). Predicted slide effect subtracted from P24215 sections.
        Assumes: the slide effect is technical (driven by gene properties); biology is not.

After each correction: disease LFC = PBS - WT (region mean), per-region rescue concordance and
slope (as 03/14), null = 84 reassignments of 3 BRICHOS among the 9 AD sections (RUVs refits its
factors per permutation because its groups depend on the labels). Residual slide = variance of
the BRICHOS-vs-BRICHOS slide contrast, relative to same-slide PBS 2-vs-1 noise splits.
"""
import sys, warnings
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold, cross_val_predict

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from utils.gene_sets import PIG_CHEN2020

TBL = ROOT / "results" / "tables" / "attenuation"
FEAT = ROOT / "data" / "external" / "gene_features" / "mm_gene_features.tsv"   # Ensembl BioMart
OUT = ROOT / "results" / "slide_removal"; OUT.mkdir(parents=True, exist_ok=True)
MIN_COUNT, PSEUDO, SEED = 10, 1.0, 0
LFC_T, PADJ_T = 0.5, 0.1
KS = (1, 2)
N_CTRL = 1000

# ================================================================================
print("[1] load pseudobulk")
counts = pd.read_csv(TBL / "pseudobulk_counts.tsv", sep="\t", index_col=0)
meta = pd.read_csv(TBL / "pseudobulk_meta.tsv", sep="\t", index_col=0)
ct = pd.crosstab(meta.region, meta["sample"])
REGIONS = sorted(ct.index[(ct > 0).all(1)])
meta = meta[meta.region.isin(REGIONS)].copy()
counts = counts.loc[meta.index]
counts = counts.loc[:, counts.sum() >= MIN_COUNT * len(meta) / 10]
meta["slide"] = meta["sample"].str[:6]
lib = counts.sum(1).values[:, None]
Y = np.log2(counts.values / lib * 1e6 + PSEUDO)          # [sample, gene]
genes = counts.columns
SECTIONS = sorted(meta["sample"].unique())
sec_t = meta.groupby("sample").treatment.first()
sec_of = meta["sample"].values; reg_of = meta.region.values
print(f"    {len(meta)} samples ({len(SECTIONS)} sections x {len(REGIONS)} regions), {len(genes)} genes")

dis_tab = pd.read_csv(TBL / "disease_signature.tsv", sep="\t", index_col=0).reindex(genes)
pig = genes.intersection(PIG_CHEN2020)
data_up = genes[(dis_tab.lfc > LFC_T) & (dis_tab.padj < PADJ_T)].union(pig)
data_dn = genes[(dis_tab.lfc < -LFC_T) & (dis_tab.padj < PADJ_T)]
SETS = {"PIG": (pig, pd.Index([])), "data": (data_up, data_dn)}

PBS = [s for s in SECTIONS if sec_t[s] == "PBS"]
BRI = [s for s in SECTIONS if sec_t[s] == "BRICHOS"]
AD = PBS + BRI
BRI24 = [s for s in BRI if s.startswith("P24215")]

def means(Yc, secs):
    """[region, gene] mean over sections."""
    return np.stack([Yc[np.isin(sec_of, secs) & (reg_of == r)].mean(0) for r in REGIONS])

def rescue(Yc, bri):
    pbs = [s for s in AD if s not in bri]
    resc = means(Yc, bri) - means(Yc, pbs)
    dis = (means(Yc, PBS) - means(Yc, ["P24215_101"])).mean(0)     # disease axis on corrected data
    out = {}
    for name, (up, dn) in SETS.items():
        iu, idn = genes.get_indexer(up), genes.get_indexer(dn)
        out[name] = np.array([(np.mean(np.concatenate([-resc[j, iu], resc[j, idn]]) > 0),
                               np.polyfit(dis[iu], resc[j, iu], 1)[0]) for j in range(len(REGIONS))])
    return out

def labels(bri):
    return np.where(np.isin(sec_of, bri), "BRICHOS", np.where(sec_of == "P24215_101", "WT", "PBS")).astype(object)

def regress_out(W, bri=None):
    """remove W from Y, fitted jointly with the treatment x region design."""
    t = labels(BRI if bri is None else bri)
    grp = pd.get_dummies(pd.Series(t + "|" + reg_of.astype(object))).values.astype(float)
    D = np.hstack([grp, W])
    B = np.linalg.lstsq(D, Y, rcond=None)[0]
    return Y - W @ B[grp.shape[1]:]

expr = Y.mean(0) > 3

# ================================================================================
print("[2] corrections")
def ruvs(k, bri):
    g = pd.Series(labels(bri) + "|" + reg_of.astype(object))
    E = Y.copy()
    for _, idx in g.groupby(g).groups.items():
        E[idx] -= Y[idx].mean(0)
    _, _, Vt = np.linalg.svd(E[:, expr], full_matrices=False)
    alpha = Vt[:k]
    W = (Y[:, expr] - Y[:, expr].mean(0)) @ alpha.T
    return W

# RUVg controls: expressed, low region variance and low PBS between-section variance, not disease genes
pbs_rows = np.isin(sec_of, PBS)
reg_var = pd.DataFrame(Y[pbs_rows]).groupby(reg_of[pbs_rows]).mean().var().values
Yr = Y.copy()
for r in REGIONS:
    Yr[reg_of == r] -= Y[reg_of == r].mean(0)            # region-centred
sec_var = Yr[pbs_rows].var(0)
cand = expr & ~genes.isin(data_up.union(data_dn))
score = pd.Series(np.argsort(np.argsort(reg_var)) + np.argsort(np.argsort(sec_var)), index=genes)[cand]
ctrl = genes.get_indexer(score.nsmallest(N_CTRL).index)
pd.Series(genes[ctrl]).to_csv(OUT / "ruvg_control_genes.csv", index=False, header=False)
_, _, Vt_g = np.linalg.svd(Yr[:, ctrl], full_matrices=False)
def ruvg(k):
    return Yr[:, ctrl] @ Vt_g[:k].T

# technical: gene features
f = pd.read_csv(FEAT, sep="\t")
f.columns = ["gene", "biotype", "chrom", "start", "end", "gc", "tx_len"]
f = (f.groupby("gene").agg(gc=("gc", "first"), gene_len=("end", "first"), start=("start", "first"),
                          tx_len=("tx_len", "median"), chrom=("chrom", "first")))
f["gene_len"] = np.log10(f.gene_len - f.start + 1); f["tx_len"] = np.log10(f.tx_len)
F = pd.DataFrame(index=genes)
F = F.join(f[["gc", "gene_len", "tx_len"]])
F["log_mean"] = Y.mean(0)
F["mt"] = genes.str.startswith("mt-").astype(int)
F["ribo"] = genes.str.match(r"^Rp[ls]\d").astype(int)
print(f"    gene features found for {F.gc.notna().mean():.0%} of genes")

def tech_fit(diff_rg, label):
    """predict [region, gene] contrast from gene features (+region); gene-wise 5-fold CV R2."""
    rows = []
    for j in range(len(REGIONS)):
        Xf = F.copy(); Xf["region"] = j; Xf["y"] = diff_rg[j]; rows.append(Xf)
    d = pd.concat(rows)
    gene_ids = np.tile(np.arange(len(genes)), len(REGIONS))
    model = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=SEED)
    folds = KFold(5, shuffle=True, random_state=SEED)
    gene_fold = np.empty(len(genes), int)
    for i, (_, te) in enumerate(folds.split(np.arange(len(genes)))): gene_fold[te] = i
    cv = np.zeros(len(d))
    for i in range(5):
        tr, te = gene_fold[gene_ids] != i, gene_fold[gene_ids] == i
        cv[te] = model.fit(d.drop(columns="y")[tr], d.y[tr]).predict(d.drop(columns="y")[te])
    w = np.tile(expr, len(REGIONS))
    r2 = 1 - np.var(d.y.values[w] - cv[w]) / np.var(d.y.values[w])
    print(f"    tech model CV R2 ({label}): {r2:.3f}")
    return cv.reshape(len(REGIONS), len(genes)), r2

slide_diff = means(Y, BRI24) - means(Y, ["P28052_401"])
s_hat, r2_slide = tech_fit(slide_diff, "slide: P24215 BRI - 401")
p28_pbs = [s for s in PBS if s.startswith("P28052")]
_, r2_noise = tech_fit(means(Y, p28_pbs[:2]) - means(Y, [p28_pbs[2]]), "noise: PBS 2-vs-1 same slide")
Y_tech = Y.copy()
for j, r in enumerate(REGIONS):
    Y_tech[(meta.slide.values == "P24215") & (reg_of == r)] -= s_hat[j]

CORR = {"none": Y, "tech": Y_tech}
for k in KS:
    CORR[f"RUVs_k{k}"] = regress_out(ruvs(k, BRI))
    CORR[f"RUVg_k{k}"] = regress_out(ruvg(k))

# ================================================================================
print("[3] residual slide + factor diagnostics")
def split_var(D):
    D = D[:, expr]; m = D.mean(0, keepdims=True)
    return np.var(m), np.var(D - m)
noise = np.array([split_var(means(Y, list(p)) - means(Y, [k_]))
                  for k_ in p28_pbs for p in combinations([s for s in p28_pbs if s != k_], 2)])
diag = []
for name, Yc in CORR.items():
    vm, vi = split_var(means(Yc, BRI24) - means(Yc, ["P28052_401"]))
    diag.append(dict(method=name, slide_var_uniform=vm, slide_var_regional=vi,
                     slide_uniform_vs_noise=vm / np.median(noise[:, 0]),
                     slide_regional_vs_noise=vi / np.median(noise[:, 1])))
for k in KS:
    for nm, W in [(f"RUVs_k{k}", ruvs(k, BRI)), (f"RUVg_k{k}", ruvg(k))]:
        sw = pd.DataFrame(W).groupby(sec_of).mean()
        tr = (sec_t.reindex(sw.index) == "BRICHOS").astype(float).values
        sl = (sw.index.str[:6] == "P24215").astype(float)
        r_t = max(abs(np.corrcoef(sw[c], tr)[0, 1]) for c in sw)
        r_s = max(abs(np.corrcoef(sw[c], sl)[0, 1]) for c in sw)
        d = [x for x in diag if x["method"] == nm][0]
        d.update(max_abs_r_factor_treatment=r_t, max_abs_r_factor_slide=r_s)
diag = pd.DataFrame(diag)
diag["tech_cv_r2_slide"] = r2_slide; diag["tech_cv_r2_noise"] = r2_noise
diag.to_csv(OUT / "correction_diagnostics.csv", index=False)
print(diag.round(3).to_string(index=False))

# ================================================================================
print("[4] rescue tests + permutation null")
perms = [list(c) for c in combinations(AD, len(BRI))]
rows = []
for name, Yc in CORR.items():
    obs = rescue(Yc, BRI)
    nulls = []
    for p in perms:
        if name.startswith("RUVs"):
            Yp = regress_out(ruvs(int(name[-1]), p), p)
        elif name.startswith("RUVg"):
            Yp = regress_out(ruvg(int(name[-1])), p)
        else:
            Yp = Yc
        nulls.append(rescue(Yp, p))
    for set_name in SETS:
        null = np.stack([n[set_name] for n in nulls])
        for j, r in enumerate(REGIONS):
            rows.append(dict(method=name, gene_set=set_name, region=r,
                             concordance=obs[set_name][j, 0], slope=obs[set_name][j, 1],
                             perm_p_conc=np.mean(null[:, j, 0] >= obs[set_name][j, 0]),
                             perm_p_slope=np.mean(null[:, j, 1] <= obs[set_name][j, 1])))
    print(f"    {name} done")
res = pd.DataFrame(rows)
res.to_csv(OUT / "rescue_after_slide_removal.csv", index=False)
with pd.option_context("display.width", 200, "display.max_rows", 300):
    for set_name in SETS:
        print(f"\n--- gene set: {set_name}")
        print(res[res.gene_set == set_name].pivot(index="region", columns="method",
              values=["concordance", "perm_p_slope"]).round(3).to_string())

# ================================================================================
print("[5] figure")
methods = list(CORR)
fig, axes = plt.subplots(1, 2, figsize=(15, 5))
for ax, val, ref in [(axes[0], "concordance", 0.5), (axes[1], "slope", 0.0)]:
    m = res[res.gene_set == "PIG"].pivot(index="region", columns="method", values=val)[methods]
    p = res[res.gene_set == "PIG"].pivot(index="region", columns="method", values="perm_p_slope")[methods]
    lim = np.nanmax(np.abs(m.values - ref))
    im = ax.imshow(m.values - ref, cmap="RdBu" if val == "concordance" else "RdBu_r",
                   vmin=-lim, vmax=lim, aspect="auto")
    for (i, j), v in np.ndenumerate(m.values):
        ax.text(j, i, f"{v:.2f}" + ("*" if p.values[i, j] < 0.05 else ""), ha="center", va="center", fontsize=8)
    ax.set_xticks(range(len(methods))); ax.set_xticklabels(methods, rotation=40, ha="right")
    ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index)
    ax.set_title(f"PIG rescue {val} after slide removal (* perm p<0.05, slope)")
    plt.colorbar(im, ax=ax, shrink=.7)
fig.tight_layout(); fig.savefig(OUT / "rescue_after_slide_removal.png", dpi=140, bbox_inches="tight")
print("done")
