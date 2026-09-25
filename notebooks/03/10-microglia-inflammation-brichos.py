#!/usr/bin/env python3
"""
10 - Microglia inflammatory-state analysis: does BRICHOS calm the microglia?

Runs on the REMOTE analysis machine in env `sc_py312`,
because that is where the full ST expression object lives:
    /Users/christoffer/work/karolinska/development/ST_BRICHOS/combined_samples.h5ad
Annotations (treatment / regions / deconv props) are joined in from the local
deconvolution bundle, exported to bundle_annotations.csv and pushed up.

Two axes of "reducing inflammation":
  (a) ABUNDANCE  - fewer microglia recruited          (deconv fraction)
  (b) STATE      - microglia shifted DAM -> homeostatic (expression signatures)

Statistics: WT (n=1) is descriptive only. The testable contrast is
PBS (n=6) vs BRICHOS (n=3). Primary test = linear mixed model with section as a
random effect (spot-level); sensitivity = section-level pseudobulk Mann-Whitney.
"""
import warnings, json
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

BASE = Path("/Users/christoffer/work/karolinska/development/ST_BRICHOS")
EXPR = BASE / "combined_samples.h5ad"
ANNO = BASE / "bundle_annotations.csv"
OUT  = BASE / "inflam_out"; FIG = OUT / "figs"; TAB = OUT / "tables"
for d in (OUT, FIG, TAB): d.mkdir(parents=True, exist_ok=True)

SAMPLE_KEY, TREAT_KEY = "sample_id", "treatment"
MG_PROP = "prop_MYE:Microglia"
Z_THRESH = 1.5                       # microglia-enriched ("hotspot") cutoff
ORDER = ["WT", "PBS", "BRICHOS"]
TCOL  = {"WT": "#7570b3", "PBS": "#d95f02", "BRICHOS": "#1b9e77"}

# ---- signatures (mouse symbols) -------------------------------------------------
SIGS = {
    "Homeostatic":      ["P2ry12","Tmem119","Cx3cr1","Hexb","Siglech","Selplg",
                          "Csf1r","Sall1","Gpr34","Olfml3","Serinc3","Cst3"],
    "DAM":              ["Trem2","Tyrobp","Cst7","Itgax","Apoe","Cd68","Lpl","Ctsb",
                          "Ctsd","Cd9","Spp1","Gpnmb","Lgals3","Clec7a","B2m","Fth1","Axl"],
    "Neuroinflammation":["Il1b","Tnf","Nfkbia","Ccl3","Ccl4","Ccl2","Cxcl10","Nlrp3",
                          "Il6","Tlr2","Ptgs2","Nfkb1","Rela","Il1a","Tnfaip3"],
    "Complement":       ["C1qa","C1qb","C1qc","C3","C4b","C1qbp","Cfh"],
    "ReactiveAstro":    ["Gfap","Vim","Serpina3n","C3","Lcn2","Cp","Steap4"],
}

# ================================================================================
print("[1] loading expression + annotations")
a = sc.read_h5ad(EXPR)
anno = pd.read_csv(ANNO, index_col=0)
common = a.obs_names.intersection(anno.index)
print(f"    expr {a.n_obs} spots, anno {len(anno)} -> {len(common)} joined")
a = a[common].copy()
for c in anno.columns:
    a.obs[c] = anno.loc[a.obs_names, c].values
a = a[a.obs[TREAT_KEY].isin(ORDER)].copy()
a.obs[SAMPLE_KEY] = a.obs[SAMPLE_KEY].astype(str)
print("    spots/treatment:\n", a.obs[TREAT_KEY].value_counts())
print("    sections/treatment:")
s2t = a.obs[[SAMPLE_KEY, TREAT_KEY]].drop_duplicates()
print(s2t.groupby(TREAT_KEY)[SAMPLE_KEY].apply(lambda s: sorted(s)).to_dict())

# ---- normalize (keep raw counts for pseudobulk DE) -----------------------------
a.layers["counts"] = a.X.copy()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

# ---- score signatures ----------------------------------------------------------
print("[2] scoring signatures")
present = {}
for name, genes in SIGS.items():
    g = [x for x in genes if x in a.var_names]
    present[name] = g
    sc.tl.score_genes(a, g, score_name=f"sig_{name}", ctrl_size=max(50, len(g)))
    print(f"    {name}: {len(g)}/{len(genes)} genes  ({g})")
a.obs["sig_DAMratio"] = a.obs["sig_DAM"] - a.obs["sig_Homeostatic"]

# ---- microglia hotspots (global z on deconv fraction) --------------------------
mg = a.obs[MG_PROP].astype(float)
a.obs["mg_z"] = (mg - mg.mean()) / mg.std()
a.obs["hotspot"] = a.obs["mg_z"] > Z_THRESH
print(f"[3] hotspots (mg_z > {Z_THRESH}): {int(a.obs['hotspot'].sum())} spots")
print(a.obs.loc[a.obs["hotspot"], TREAT_KEY].value_counts())

SIG_COLS = [f"sig_{n}" for n in SIGS] + ["sig_DAMratio"]
h = a.obs.loc[a.obs["hotspot"]].copy()          # hotspot spots only

# ================================================================================
# [4] STATE: mixed model (BRICHOS vs PBS) + pseudobulk sensitivity, in hotspots
# ================================================================================
print("[4] mixed models on signature scores (hotspots, BRICHOS vs PBS)")
pb = h.groupby([SAMPLE_KEY, TREAT_KEY], observed=True)[SIG_COLS + ["mg_z", MG_PROP]].mean().reset_index()

def mixed_test(df, y, groupcol=SAMPLE_KEY):
    d = df[df[TREAT_KEY].isin(["PBS", "BRICHOS"])].copy()
    d[TREAT_KEY] = pd.Categorical(d[TREAT_KEY], categories=["PBS", "BRICHOS"])
    try:
        m = smf.mixedlm(f"{y} ~ {TREAT_KEY}", d, groups=d[groupcol]).fit(reml=False)
        term = [t for t in m.params.index if "BRICHOS" in t][0]
        return dict(beta=float(m.params[term]), se=float(m.bse[term]),
                    p=float(m.pvalues[term]), n=int(len(d)))
    except Exception as e:
        return dict(beta=np.nan, se=np.nan, p=np.nan, n=0, err=str(e))

def pseudobulk_test(y):
    x = pb[pb[TREAT_KEY] == "PBS"][y].values
    z = pb[pb[TREAT_KEY] == "BRICHOS"][y].values
    try:  u, p = mannwhitneyu(x, z, alternative="two-sided")
    except Exception: u, p = np.nan, np.nan
    return dict(PBS_mean=float(np.mean(x)), BRICHOS_mean=float(np.mean(z)),
                delta=float(np.mean(z) - np.mean(x)), MWU_p=float(p),
                n_PBS=len(x), n_BRICHOS=len(z))

rows = []
for y in SIG_COLS + ["mg_z", MG_PROP]:
    mm = mixed_test(h, y)                 # h = all hotspot spots, has every column
    pt = pseudobulk_test(y)
    # WT descriptive
    wt = pb[pb[TREAT_KEY] == "WT"][y]
    rows.append({"metric": y, **{f"mm_{k}": v for k, v in mm.items()},
                 **pt, "WT_mean": float(wt.mean()) if len(wt) else np.nan})
res = pd.DataFrame(rows)
res.to_csv(TAB / "state_mixedmodel_pseudobulk.csv", index=False)
print(res[["metric","mm_beta","mm_p","PBS_mean","BRICHOS_mean","delta","MWU_p","WT_mean"]].to_string(index=False))
pb.to_csv(TAB / "state_section_pseudobulk_means.csv", index=False)

# ---- figure: signature scores per treatment (hotspots), section means overlaid --
print("[5] figure: state by treatment")
panels = ["sig_Homeostatic","sig_DAM","sig_DAMratio","sig_Neuroinflammation","sig_Complement"]
fig, axes = plt.subplots(1, len(panels), figsize=(4.0*len(panels), 4.6))
for ax, y in zip(axes, panels):
    data = [h[h[TREAT_KEY] == t][y].values for t in ORDER]
    parts = ax.violinplot(data, showmeans=False, showextrema=False)
    for pc, t in zip(parts["bodies"], ORDER):
        pc.set_facecolor(TCOL[t]); pc.set_alpha(.35)
    for i, t in enumerate(ORDER, 1):
        pts = pb[pb[TREAT_KEY] == t][y].values
        ax.scatter(np.full(len(pts), i), pts, color=TCOL[t], edgecolor="k",
                   s=45, zorder=3)
    ax.set_xticks(range(1, len(ORDER)+1)); ax.set_xticklabels(ORDER, rotation=30)
    row = res[res.metric == y].iloc[0]
    ax.set_title(f"{y.replace('sig_','')}\nMM p(BRvsPBS)={row.mm_p:.3g}", fontsize=10)
    ax.axhline(0, color="grey", lw=.5, ls=":")
fig.suptitle(f"Microglia state in enriched spots (mg_z>{Z_THRESH})  -  points = section means",
             fontsize=13)
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(FIG / "inflam_state_by_treatment.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# ---- figure: regional breakdown of DAM:homeostatic ratio -----------------------
print("[6] figure: regional DAM:homeostatic ratio")
if "re_annotation_regions" in h.columns:
    reg = (h.groupby(["re_annotation_regions", TREAT_KEY], observed=True)["sig_DAMratio"]
             .mean().reset_index())
    regions = sorted(h["re_annotation_regions"].dropna().unique())
    fig, ax = plt.subplots(figsize=(max(8, 0.7*len(regions)), 5))
    w = 0.25
    for i, t in enumerate(ORDER):
        vals = [reg[(reg.re_annotation_regions == r) & (reg[TREAT_KEY] == t)]["sig_DAMratio"]
                for r in regions]
        vals = [v.iloc[0] if len(v) else np.nan for v in vals]
        ax.bar(np.arange(len(regions)) + (i-1)*w, vals, w, label=t, color=TCOL[t])
    ax.set_xticks(range(len(regions))); ax.set_xticklabels(regions, rotation=45, ha="right")
    ax.axhline(0, color="k", lw=.6); ax.set_ylabel("DAM : homeostatic ratio (mean)")
    ax.set_title("Regional microglia inflammatory state (hotspots)"); ax.legend()
    fig.tight_layout(); fig.savefig(FIG / "inflam_regional_DAMratio.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    reg.to_csv(TAB / "regional_DAMratio.csv", index=False)

# ================================================================================
# [7] DE: pseudobulk hotspot counts, BRICHOS vs PBS (pydeseq2) + GSEA prerank
# ================================================================================
print("[7] pseudobulk DE (pydeseq2) BRICHOS vs PBS in hotspots")
import scipy.sparse as sp
hot = a[a.obs["hotspot"].values & a.obs[TREAT_KEY].isin(["PBS","BRICHOS"]).values]
counts = hot.layers["counts"]
counts = counts.toarray() if sp.issparse(counts) else np.asarray(counts)
cdf = pd.DataFrame(counts, index=hot.obs_names, columns=hot.var_names)
cdf[SAMPLE_KEY] = hot.obs[SAMPLE_KEY].values
pbc = cdf.groupby(SAMPLE_KEY).sum().round().astype(int)                 # section x gene
meta = hot.obs[[SAMPLE_KEY, TREAT_KEY]].drop_duplicates().set_index(SAMPLE_KEY)
meta = meta.loc[pbc.index]
keep = (pbc > 5).sum(axis=0) >= 3                                       # expressed filter
pbc = pbc.loc[:, keep]
print(f"    pseudobulk matrix {pbc.shape}, sections: {list(pbc.index)}")

de_ok = False
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    meta[TREAT_KEY] = pd.Categorical(meta[TREAT_KEY], categories=["PBS","BRICHOS"])
    try:
        dds = DeseqDataSet(counts=pbc, metadata=meta, design="~" + TREAT_KEY)
    except TypeError:
        dds = DeseqDataSet(counts=pbc, metadata=meta, design_factors=TREAT_KEY)
    dds.deseq2()
    st = DeseqStats(dds, contrast=[TREAT_KEY, "BRICHOS", "PBS"]); st.summary()
    de = st.results_df.sort_values("stat")
    de.to_csv(TAB / "de_hotspot_BRICHOS_vs_PBS.csv")
    de_ok = True
    sig = de[(de.padj < 0.1)]
    print(f"    DE genes padj<0.1: {len(sig)}  (down in BRICHOS: {(sig.stat<0).sum()}, up: {(sig.stat>0).sum()})")
    # volcano
    fig, ax = plt.subplots(figsize=(6,6))
    x = de["log2FoldChange"]; y = -np.log10(de["pvalue"].clip(lower=1e-300))
    ax.scatter(x, y, s=6, c="lightgrey")
    m = de.padj < 0.1
    ax.scatter(x[m], y[m], s=10, c="crimson")
    for g in de.index[m][:20]:
        ax.annotate(g, (de.loc[g,"log2FoldChange"], -np.log10(max(de.loc[g,"pvalue"],1e-300))),
                    fontsize=7)
    ax.axvline(0, color="grey", lw=.5); ax.set_xlabel("log2FC (BRICHOS / PBS)")
    ax.set_ylabel("-log10 p"); ax.set_title("Hotspot pseudobulk DE: BRICHOS vs PBS")
    fig.tight_layout(); fig.savefig(FIG / "inflam_de_volcano.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
except Exception as e:
    print("    DE failed:", e)

# ---- GSEA prerank against curated inflammatory sets (offline) ------------------
if de_ok:
    print("[8] GSEA prerank on DE ranking vs curated inflammatory sets")
    try:
        import gseapy as gp
        rnk = de["stat"].dropna().sort_values(ascending=False)
        gene_sets = {k: [g for g in v if g in a.var_names] for k, v in SIGS.items()}
        pre = gp.prerank(rnk=rnk.reset_index().rename(columns={"index":0,"stat":1}),
                         gene_sets=gene_sets, min_size=3, max_size=500,
                         permutation_num=1000, seed=0, no_plot=True,
                         outdir=str(OUT / "gsea"))
        gres = pre.res2d
        gres.to_csv(TAB / "gsea_prerank_inflammatory.csv", index=False)
        print(gres[["Term","NES","NOM p-val","FDR q-val"]].to_string(index=False))
    except Exception as e:
        print("    GSEA failed:", e)

# ---- export per-spot scores (so spatial maps can be redrawn locally on H&E) -----
obs_export = a.obs[[SAMPLE_KEY, TREAT_KEY, "re_annotation_regions", MG_PROP,
                    "mg_z", "hotspot"] + SIG_COLS].copy()
obs_export.to_csv(TAB / "per_spot_scores.csv")

# ================================================================================
# [9] spatial maps of DAM:homeostatic ratio, all sections, grouped by treatment
# ================================================================================
print("[9] spatial DAM:homeostatic ratio maps")
try:
    samples = []
    for t in ORDER:
        samples += sorted(s2t.loc[s2t[TREAT_KEY] == t, SAMPLE_KEY].astype(str))
    n = len(samples); ncol = 5; nrow = int(np.ceil(n/ncol)); vmax = 1.0
    have_img = "spatial" in a.uns
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2*ncol, 4.4*nrow), squeeze=False)
    for ax, s in zip(np.ravel(axes), samples):
        m = (a.obs[SAMPLE_KEY] == s).values
        xy = a.obsm["spatial"][m]
        if have_img and s in a.uns["spatial"]:
            sp_ = a.uns["spatial"][s]; sf = sp_["scalefactors"]["tissue_hires_scalef"]
            ax.imshow(sp_["images"]["hires"], origin="upper", alpha=0.45); xy = xy*sf
        sc = ax.scatter(xy[:,0], xy[:,1], s=8, c=a.obs["sig_DAMratio"].values[m],
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax, edgecolors="none")
        t = a.obs.loc[m, TREAT_KEY].iloc[0]
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_title(f"{s} [{t}]", fontsize=10, color=TCOL.get(t,"k"), fontweight="bold")
    for ax in np.ravel(axes)[n:]: ax.axis("off")
    cb = fig.colorbar(sc, ax=np.ravel(axes).tolist(), fraction=0.02, pad=0.01)
    cb.set_label("DAM : homeostatic ratio")
    fig.suptitle("Microglia inflammatory state (DAM:homeostatic) across sections", fontsize=14)
    fig.savefig(FIG / "inflam_spatial_DAMratio_grid.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
except Exception as e:
    print("    spatial map failed:", e)

# ---- write a small run summary -------------------------------------------------
summary = dict(
    z_thresh=Z_THRESH,
    n_spots=int(a.n_obs),
    n_hotspots=int(a.obs["hotspot"].sum()),
    hotspots_per_treatment=a.obs.loc[a.obs["hotspot"], TREAT_KEY].value_counts().to_dict(),
    signatures={k: v for k, v in present.items()},
    de_ok=de_ok,
)
(OUT / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print("[done] outputs in", OUT)
