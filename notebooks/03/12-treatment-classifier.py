#!/usr/bin/env python3
"""
12 - Can we predict a section's treatment from its spots?  (and is it just the slide?)

Runs on the REMOTE analysis machine in env `sc_py312` (full expression lives there),
same layout as notebooks 03/10*. Annotations joined from bundle_annotations.csv,
gene sets from utils/gene_sets.py (pushed alongside).

Design
------
* Unit of evaluation = SECTION, not spot. Spots of one section share slide, depth,
  RNA quality, so a random spot split leaks section identity. We use
  leave-one-section-out (LOSO): train on all spots of the other sections, predict
  every spot of the held-out section, section score = mean p(BRICHOS).
* Unsupervised preprocessing (HVG selection, scaling, PCA) is re-fit per fold on
  training sections only. It uses no labels, so the per-fold matrices are cached
  and reused across tasks and label permutations.
* Sample weights: every section contributes equal total weight and both classes
  contribute equal total weight (otherwise LOSO with 6 v 3 just predicts PBS).
* Feature sets, most to least constrained:
    qc   - log counts, log genes, %mt, %ribo         (technical negative control)
    sig  - curated signature scores (gene_sets.SIGNATURE_REGISTRY)
    pca  - 30 PCs of per-fold HVGs
    hvg  - 2000 per-fold HVGs, scaled, L2 logistic
* Tasks:
    T1 treatment  PBS vs BRICHOS (9 sections, WT excluded)
    T2 slide      P24215 vs P28052 (10 sections) - batch baseline
    WT projection: T1 model trained on all 9 AD sections, applied to WT.
* Confound tests (treatment ~ slide, see memory batch-confound):
    perm  - all C(9,3)=84 relabelings of which sections are BRICHOS; section AUC null.
    401   - slide-conditional test: keep P24215_201/301 as BRICHOS and let each of the
            7 P28052 sections in turn play the 3rd BRICHOS. Is the real one (401)
            more recognisable than a mislabelled PBS section? Min p = 1/7.

Run 1 (no flags) showed the classifier mostly tracks sequencing depth (P24215 and
P28052_401 are the shallowest sections), cutting plane (hypothalamic neuropeptides)
and immediate-early genes. The --downsample / --drop-ieg / --match-regions /
--drop-anatomy flags
remove those; outputs go to classifier_out_<tag>.
"""
import warnings, sys, itertools, json
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

BASE = Path("/Users/christoffer/work/karolinska/development/ST_BRICHOS")
EXPR = BASE / "combined_samples.h5ad"
ANNO = BASE / "bundle_annotations.csv"
sys.path.insert(0, str(BASE / "classifier_utils"))
from gene_sets import SIGNATURE_REGISTRY, filter_to_var

# ---- confound removal (all off = the original run) ------------------------------
#   --downsample N   : binomially downsample every spot to exactly N UMIs, drop spots < N
#                      (removes sequencing depth by construction)
#   --drop-ieg       : drop immediate-early / handling-stress genes before everything
#   --match-regions  : keep only regions with >= MIN_REGION_SPOTS in every section, and
#                      give each region equal weight within a section (training weights
#                      AND held-out section score), so cutting plane / region mix can't
#                      separate sections
#   --drop-anatomy   : drop the Hypothalamus region and plane-sensitive marker genes
#                      (hypothalamic neuropeptides / TFs, choroid plexus) - sub-region
#                      cutting plane survives --match-regions
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--downsample", type=int, default=0)
ap.add_argument("--drop-ieg", action="store_true")
ap.add_argument("--match-regions", action="store_true")
ap.add_argument("--drop-anatomy", action="store_true")
ap.add_argument("--tag", default="")
ARGS = ap.parse_args()
OUT  = BASE / ("classifier_out" + (f"_{ARGS.tag}" if ARGS.tag else "")); OUT.mkdir(exist_ok=True)
MIN_REGION_SPOTS = 30
IEG = ["Fos", "Fosb", "Fosl2", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3", "Egr4", "Arc",
       "Nr4a1", "Nr4a2", "Nr4a3", "Ier2", "Ier3", "Ier5", "Dusp1", "Dusp6", "Klf2", "Klf4",
       "Klf6", "Atf3", "Npas4", "Btg2", "Cyr61", "Ccn1", "Ppp1r15a", "Gadd45b", "Gadd45g",
       "Zfp36", "Hspa1a", "Hspa1b", "Hspb1", "Socs3", "Cebpb", "Cebpd"]
ANATOMY_REGIONS = ["Hypothalamus"]
ANATOMY_GENES = [  # hypothalamic neuropeptides / hormones / TFs, choroid plexus
    "Oxt", "Avp", "Pmch", "Hcrt", "Pomc", "Agrp", "Npy", "Cartpt", "Trh", "Crh", "Gal",
    "Ghrh", "Gnrh1", "Kiss1", "Nts", "Tac2", "Sst", "Prl", "Gh", "Pnoc", "Nmb", "Penk",
    "Sim1", "Otp", "Fezf1", "Nr5a1", "Lhx6", "Lhx8", "Shox2", "Six3", "Foxb1", "Pitx2",
    "Ttr", "Kl", "Folr1", "Clic6", "Kcnj13", "Prlr", "Sostdc1", "Enpp2", "Calml4"]

SAMPLE_KEY, TREAT_KEY, REGION_KEY = "sample_id", "treatment", "re_annotation_regions"
N_TRAIN_PER_SECTION = 800          # spots are pseudo-replicates; subsample for speed
N_HVG, N_PC, C_REG = 2000, 30, 0.05
SEED = 0
FEATURE_SETS = ["qc", "sig", "pca", "hvg"]
print("config:", vars(ARGS))

# ================================================================================
print("[1] load + join")
a = sc.read_h5ad(EXPR)
anno = pd.read_csv(ANNO, index_col=0)
a = a[a.obs_names.intersection(anno.index)].copy()
for c in [SAMPLE_KEY, TREAT_KEY, REGION_KEY]:
    a.obs[c] = anno.loc[a.obs_names, c].astype(str).values
a.obs["slide"] = a.obs[SAMPLE_KEY].str[:6]
print(pd.crosstab(a.obs[SAMPLE_KEY], [a.obs["slide"], a.obs[TREAT_KEY]]))

if ARGS.drop_ieg:
    drop = a.var_names.isin(IEG)
    print(f"    dropping {drop.sum()} IEG/stress genes")
    a = a[:, ~drop].copy()
if ARGS.drop_anatomy:
    drop = a.var_names.isin(ANATOMY_GENES)
    print(f"    dropping {drop.sum()} anatomy genes + regions {ANATOMY_REGIONS}")
    a = a[~a.obs[REGION_KEY].isin(ANATOMY_REGIONS)][:, ~drop].copy()
if ARGS.downsample:
    n0 = a.n_obs
    a = a[np.asarray(a.X.sum(1)).ravel() >= ARGS.downsample].copy()
    sc.pp.downsample_counts(a, counts_per_cell=ARGS.downsample, random_state=SEED)
    print(f"    downsampled to {ARGS.downsample} UMI/spot; kept {a.n_obs}/{n0} spots")
    print("    kept per section:", a.obs[SAMPLE_KEY].value_counts().sort_index().to_dict())
if ARGS.match_regions:
    ct = pd.crosstab(a.obs[REGION_KEY], a.obs[SAMPLE_KEY])
    shared = ct.index[(ct >= MIN_REGION_SPOTS).all(1)]
    dropped = sorted(set(ct.index) - set(shared))
    a = a[a.obs[REGION_KEY].isin(shared)].copy()
    print(f"    shared regions ({len(shared)}): {list(shared)}")
    print(f"    dropped regions: {dropped}; {a.n_obs} spots left")

a.var["mt"] = a.var_names.str.startswith("mt-")
a.var["ribo"] = a.var_names.str.match(r"^Rp[sl]\d")
sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None, inplace=True)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

obs = a.obs
sections = sorted(obs[SAMPLE_KEY].unique())
sec_of = obs[SAMPLE_KEY].values
reg_of = obs[REGION_KEY].values

def section_score(p, idx):
    """held-out section score: plain mean, or region-balanced mean with --match-regions"""
    if not ARGS.match_regions: return p.mean()
    return pd.Series(p).groupby(reg_of[idx]).mean().mean()

# fixed training subsample (same spots in every fold)
rng = np.random.default_rng(SEED)
train_pool = np.zeros(a.n_obs, bool)
for s in sections:
    idx = np.where(sec_of == s)[0]
    train_pool[rng.choice(idx, min(N_TRAIN_PER_SECTION, len(idx)), replace=False)] = True

# ---- fold-independent features ---------------------------------------------------
print("[2] signature + QC features")
qc = np.column_stack([np.log1p(obs.total_counts), np.log1p(obs.n_genes_by_counts),
                      obs.pct_counts_mt, obs.pct_counts_ribo]).astype(np.float32)
sig_names = []
for name, genes in SIGNATURE_REGISTRY.items():
    g = filter_to_var(genes, a.var_names)
    if len(g) >= 3:
        sc.tl.score_genes(a, g, score_name=f"sig_{name}", ctrl_size=max(50, len(g)), random_state=SEED)
        sig_names.append(f"sig_{name}")
sig = obs[sig_names].values.astype(np.float32)
print(f"    {len(sig_names)} signatures")

# ---- per-fold features (fit on training sections only) ---------------------------
def build_fold(held):
    """Return dict feature_set -> (X_train_pool_rows, X_heldout_rows) for held-out section
    (held=None -> fit on all sections, X_heldout = all spots)."""
    tr = train_pool & (sec_of != held) if held else train_pool
    te = (sec_of == held) if held else np.ones(a.n_obs, bool)
    sub = a[tr]
    hv = sc.pp.highly_variable_genes(sub, n_top_genes=N_HVG, flavor="seurat", inplace=False)
    genes = sub.var_names[hv["highly_variable"].values]
    def dense(m): return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float32)
    Xtr, Xte = dense(a[tr, genes].X), dense(a[te, genes].X)
    scl = StandardScaler().fit(Xtr)
    Xtr, Xte = np.clip(scl.transform(Xtr), -10, 10), np.clip(scl.transform(Xte), -10, 10)
    pca = PCA(N_PC, random_state=SEED).fit(Xtr)
    feats = {"hvg": (Xtr, Xte), "pca": (pca.transform(Xtr), pca.transform(Xte))}
    for k, M in (("qc", qc), ("sig", sig)):
        s = StandardScaler().fit(M[tr]); feats[k] = (s.transform(M[tr]), s.transform(M[te]))
    return dict(tr=np.where(tr)[0], te=np.where(te)[0], feats=feats, genes=list(genes))

print("[3] building per-fold features")
FOLDS = {s: build_fold(s) for s in sections}
FOLDS["ALL"] = build_fold(None)

# ---- core fit/predict -------------------------------------------------------------
def weights(sec, y, reg=None):
    """equal total weight per section (and per region within a section if reg given),
    equal total weight per class"""
    w = np.ones(len(y))
    for s in np.unique(sec):
        m = sec == s
        if reg is None: w[m] = 1.0 / m.sum(); continue
        rs = np.unique(reg[m])
        for r in rs: w[m & (reg == r)] = 1.0 / ((m & (reg == r)).sum() * len(rs))
    for c in (0, 1):
        n_sec = len(np.unique(sec[y == c])); w[y == c] /= max(n_sec, 1)
    return w * len(y) / w.sum()

def fit(fold, fs, label_of_section, train_sections):
    tr_idx = fold["tr"]; Xtr = fold["feats"][fs][0]
    keep = np.isin(sec_of[tr_idx], train_sections)
    y = np.array([label_of_section[s] for s in sec_of[tr_idx][keep]])
    m = LogisticRegression(C=C_REG, max_iter=3000)
    reg = reg_of[tr_idx][keep] if ARGS.match_regions else None
    m.fit(Xtr[keep], y, sample_weight=weights(sec_of[tr_idx][keep], y, reg))
    return m

def loso(fs, label_of_section, task_sections):
    """held-out section -> (mean p, frac spots > .5, per-spot p)"""
    out = {}
    for s in task_sections:
        f = FOLDS[s]
        m = fit(f, fs, label_of_section, [t for t in task_sections if t != s])
        p = m.predict_proba(f["feats"][fs][1])[:, 1]
        out[s] = (section_score(p, f["te"]), (p > .5).mean(), p)
    return out

def section_auc(res, label_of_section):
    secs = list(res); y = [label_of_section[s] for s in secs]
    return roc_auc_score(y, [res[s][0] for s in secs])

AD = [s for s in sections if obs.loc[sec_of == s, TREAT_KEY].iloc[0] != "WT"]
WT = [s for s in sections if s not in AD]
treat = {s: int(obs.loc[sec_of == s, TREAT_KEY].iloc[0] == "BRICHOS") for s in sections}
slide = {s: int(s.startswith("P24215")) for s in sections}

# ================================================================================
print("[4] LOSO: treatment (T1) and slide (T2)")
rows, spot_rows, summ = [], [], []
for fs in FEATURE_SETS:
    for task, lab, secs in (("treatment", treat, AD), ("slide", slide, sections)):
        res = loso(fs, lab, secs)
        acc = np.mean([(res[s][0] > .5) == lab[s] for s in secs])
        auc = section_auc(res, lab)
        summ.append(dict(feature_set=fs, task=task, n_sections=len(secs), section_acc=acc, section_auc=auc))
        print(f"    {fs:4s} {task:9s} acc={acc:.2f} auc={auc:.2f}")
        for s in secs:
            rows.append(dict(feature_set=fs, task=task, section=s, slide=s[:6],
                             treatment=obs.loc[sec_of == s, TREAT_KEY].iloc[0],
                             label=lab[s], mean_p=res[s][0], frac_gt_05=res[s][1]))
            if task == "treatment":
                spot_rows.append(pd.DataFrame({"barcode": a.obs_names[FOLDS[s]["te"]],
                                               "feature_set": fs, "p_BRICHOS": res[s][2]}))
    # WT projection: train on all 9 AD sections (WT fold = features fit without WT)
    for w in WT:
        f = FOLDS[w]; m = fit(f, fs, treat, AD)
        p = m.predict_proba(f["feats"][fs][1])[:, 1]
        rows.append(dict(feature_set=fs, task="WT_projection", section=w, slide=w[:6],
                         treatment="WT", label=np.nan, mean_p=section_score(p, f["te"]), frac_gt_05=(p > .5).mean()))
        spot_rows.append(pd.DataFrame({"barcode": a.obs_names[f["te"]], "feature_set": fs, "p_BRICHOS": p}))
pd.DataFrame(rows).to_csv(OUT / "section_predictions.csv", index=False)
pd.DataFrame(summ).to_csv(OUT / "loso_summary.csv", index=False)
pd.concat(spot_rows).to_csv(OUT / "spot_predictions_heldout.csv.gz", index=False)

# ================================================================================
print("[5] permutation null: all C(9,3) BRICHOS labelings")
labelings = list(itertools.combinations(AD, 3))
def perm_one(fs, combo):
    lab = {s: int(s in combo) for s in AD}
    return dict(feature_set=fs, brichos=";".join(combo),
                n_on_P24215=sum(s.startswith("P24215") for s in combo),
                is_true=set(combo) == {s for s in AD if treat[s]},
                section_auc=section_auc(loso(fs, lab, AD), lab))
perm = pd.DataFrame(Parallel(n_jobs=8, prefer="threads")(
    delayed(perm_one)(fs, c) for fs in FEATURE_SETS for c in labelings))
perm.to_csv(OUT / "permutation_null.csv", index=False)
for fs in FEATURE_SETS:
    d = perm[perm.feature_set == fs]; t = d[d.is_true].section_auc.iloc[0]
    print(f"    {fs:4s} true AUC={t:.2f}  p_perm={(d.section_auc >= t).mean():.3f}  "
          f"(slide-aligned labelings mean AUC={d[d.n_on_P24215 == 2].section_auc.mean():.2f})")

# ================================================================================
print("[6] slide-conditional test for P28052_401")
p24_bri = [s for s in AD if s.startswith("P24215") and treat[s]]
cands = [s for s in AD if s.startswith("P28052")]
cond = []
for fs in FEATURE_SETS:
    for c in cands:
        lab = {s: int(s in p24_bri or s == c) for s in AD}
        f = FOLDS[c]; m = fit(f, fs, lab, [t for t in AD if t != c])
        cond.append(dict(feature_set=fs, candidate=c, is_true=treat[c] == 1,
                         heldout_mean_p=section_score(m.predict_proba(f["feats"][fs][1])[:, 1], f["te"])))
cond = pd.DataFrame(cond)
cond["rank"] = cond.groupby("feature_set").heldout_mean_p.rank(ascending=False)
cond.to_csv(OUT / "slide_conditional_401.csv", index=False)
for fs in FEATURE_SETS:
    d = cond[cond.feature_set == fs]; r = d[d.is_true]["rank"].iloc[0]
    print(f"    {fs:4s} 401 rank {int(r)}/{len(d)}  p={r/len(d):.2f}")

# ================================================================================
print("[7] what does each model look at? (ALL-fold hvg, treatment vs slide)")
f = FOLDS["ALL"]
coef = {}
for task, lab, secs in (("treatment", treat, AD), ("slide", slide, sections)):
    coef[task] = fit(f, "hvg", lab, secs).coef_[0]
cf = pd.DataFrame(coef, index=f["genes"])
cf.to_csv(OUT / "hvg_coefficients.csv")
r = np.corrcoef(cf.treatment, cf.slide)[0, 1]
print(f"    corr(treatment coef, slide coef) = {r:.2f}")
print("    top BRICHOS-ward genes:", list(cf.treatment.nlargest(15).index))
print("    top PBS-ward genes:    ", list(cf.treatment.nsmallest(15).index))
json.dump(dict(coef_corr_treatment_vs_slide=r, n_train_per_section=N_TRAIN_PER_SECTION,
               n_hvg=N_HVG, n_pc=N_PC, C=C_REG), open(OUT / "meta.json", "w"), indent=1)
print("done")
