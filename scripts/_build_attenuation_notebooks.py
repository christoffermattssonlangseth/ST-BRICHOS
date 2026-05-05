"""Generate notebooks/03/*.ipynb implementing the attenuation stats roadmap."""
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "notebooks" / "03"
OUT_DIR.mkdir(parents=True, exist_ok=True)


COMMON_HEADER_CODE = """\
from __future__ import annotations
import sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

ROOT = Path.cwd().resolve()
while not (ROOT / 'utils').exists():
    if ROOT.parent == ROOT:
        raise RuntimeError('could not locate project root')
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))

# regional-annotation h5ad lives on the processing volume
BASEDIR = Path('/Volumes/processing2/ST_BRICHOS/data')
H5AD_ORIENTED = BASEDIR / 'ST_BRICHOS_region_subcluster_oriented.h5ad'
H5AD_BASE     = BASEDIR / 'ST_BRICHOS_region_subcluster.h5ad'
H5AD = H5AD_ORIENTED if H5AD_ORIENTED.exists() else H5AD_BASE
COUNT_LAYER = 'counts'   # raw integer counts live here, not in .X

TBL = ROOT / 'results' / 'tables' / 'attenuation'
TBL.mkdir(parents=True, exist_ok=True)
FIG = ROOT / 'results' / 'figures' / 'manuscript'
FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_KEY    = 'sample_id'           # change to 'library_id' if obs uses that
REGION_KEY    = 'anatomical_region'   # adjust to your obs column for regions
TREATMENT_KEY = 'treatment'
print('h5ad        :', H5AD)
print('count layer :', COUNT_LAYER)
"""


def code(src: str):
    return nbf.v4.new_code_cell(src.rstrip() + "\n")


def md(src: str):
    return nbf.v4.new_markdown_cell(src.rstrip())


def write(name: str, cells):
    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python"},
    }
    out = OUT_DIR / name
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"wrote {out}")


# -----------------------------------------------------------------------------
# 01 — pseudobulk
# -----------------------------------------------------------------------------

write("01-pseudobulk.ipynb", [
    md(
        "# 03/01 — Pseudobulk per (mouse, region)\n"
        "\n"
        "Collapse spot-level counts to one row per (sample x region). "
        "Everything downstream (DGE, attenuation slope, mixed model) "
        "consumes this table, so the unit of replication is the *mouse* "
        "rather than the spot.\n"
        "\n"
        "Outputs:\n"
        "* `results/tables/attenuation/pseudobulk_counts.tsv`\n"
        "* `results/tables/attenuation/pseudobulk_meta.tsv`"
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from utils.attenuation import make_pseudobulk\n"
        "\n"
        "adata = sc.read_h5ad(H5AD)\n"
        "print(adata)\n"
        "for k in (SAMPLE_KEY, REGION_KEY, TREATMENT_KEY):\n"
        "    if k not in adata.obs.columns:\n"
        "        cands = [c for c in adata.obs.columns if k.split('_')[0] in c]\n"
        "        print(f'!! {k!r} missing — candidates: {cands}')\n"
    ),
    md(
        "### Sum spot counts to (sample, region)\n"
        "\n"
        "Raw counts live in `adata.layers['counts']`. `make_pseudobulk` "
        "reads that layer directly — no need to mutate `.X`."
    ),
    code(
        "counts, meta = make_pseudobulk(\n"
        "    adata,\n"
        "    sample_key=SAMPLE_KEY,\n"
        "    region_key=REGION_KEY,\n"
        "    treatment_key=TREATMENT_KEY,\n"
        "    min_spots=20,\n"
        "    layer=COUNT_LAYER,\n"
        ")\n"
        "# cast to integer counts (matches DESeq2/edgeR expectations downstream)\n"
        "counts = counts.round().astype(int)\n"
        "print(counts.shape, meta.shape)\n"
        "meta.head()\n"
    ),
    md("### Sanity"),
    code(
        "print('mice per treatment:')\n"
        "print(meta.groupby('treatment')['sample'].nunique())\n"
        "print()\n"
        "print('regions per treatment:')\n"
        "print(meta.groupby(['treatment','region']).size().unstack(fill_value=0))\n"
    ),
    md("### Persist"),
    code(
        "counts.to_csv(TBL / 'pseudobulk_counts.tsv', sep='\\t')\n"
        "meta.to_csv  (TBL / 'pseudobulk_meta.tsv',   sep='\\t')\n"
        "print('wrote', TBL / 'pseudobulk_counts.tsv')\n"
        "print('wrote', TBL / 'pseudobulk_meta.tsv')\n"
    ),
])


# -----------------------------------------------------------------------------
# 02 — attenuation stats (T1, T2, T3)
# -----------------------------------------------------------------------------

write("02-attenuation-stats.ipynb", [
    md(
        "# 03/02 — Attenuation statistics per region\n"
        "\n"
        "Three tests on pseudobulk:\n"
        "\n"
        "* **T1 — interaction LFC**  : per-gene Welch on log-CPM, contrast "
        "BRI vs PBS.\n"
        "* **T2 — sign concordance** : of disease-DEGs, fraction that "
        "reverse direction under BRICHOS (binomial test vs 50%).\n"
        "* **T3 — attenuation slope**: regress LFC_residual on "
        "LFC_disease over disease-DEGs; slope <1 = rescue, with bootstrap CI.\n"
        "\n"
        "Outputs `results/tables/attenuation/regional_attenuation_stats.tsv`."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from utils.attenuation import (\n"
        "    pseudobulk_lfc, sign_concordance, attenuation_slope\n"
        ")\n"
        "\n"
        "counts = pd.read_csv(TBL / 'pseudobulk_counts.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "meta   = pd.read_csv(TBL / 'pseudobulk_meta.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "regions = sorted(meta['region'].dropna().unique())\n"
        "print(len(regions), 'regions')\n"
    ),
    md("### Per-region tests"),
    code(
        "rows = []\n"
        "padj_thresh = 0.05\n"
        "for region in regions:\n"
        "    try:\n"
        "        pbs_v_wt = pseudobulk_lfc(counts, meta,\n"
        "                                  group_a='PBS', group_b='WT',\n"
        "                                  region=region)\n"
        "        bri_v_wt = pseudobulk_lfc(counts, meta,\n"
        "                                  group_a='BRICHOS', group_b='WT',\n"
        "                                  region=region)\n"
        "        bri_v_pbs = pseudobulk_lfc(counts, meta,\n"
        "                                   group_a='BRICHOS', group_b='PBS',\n"
        "                                   region=region)\n"
        "    except ValueError as e:\n"
        "        print(f'  {region}: skip ({e})')\n"
        "        continue\n"
        "\n"
        "    common = (pbs_v_wt.lfc.index\n"
        "              .intersection(bri_v_wt.lfc.index)\n"
        "              .intersection(bri_v_pbs.lfc.index))\n"
        "    sig = pbs_v_wt.padj.loc[common] < padj_thresh\n"
        "\n"
        "    sc_res = sign_concordance(pbs_v_wt.lfc.loc[common],\n"
        "                              bri_v_pbs.lfc.loc[common],\n"
        "                              sig)\n"
        "    slope_res = attenuation_slope(pbs_v_wt.lfc.loc[common],\n"
        "                                  bri_v_wt.lfc.loc[common],\n"
        "                                  sig, n_boot=1000)\n"
        "\n"
        "    rows.append(dict(\n"
        "        region=region,\n"
        "        n_disease_DEG=int(sig.sum()),\n"
        "        n_reverse=sc_res['n_reverse'],\n"
        "        frac_reverse=sc_res['fraction'],\n"
        "        binomial_p=sc_res['p'],\n"
        "        slope=slope_res['slope'],\n"
        "        slope_ci_low=slope_res['ci_low'],\n"
        "        slope_ci_high=slope_res['ci_high'],\n"
        "        intercept=slope_res['intercept'],\n"
        "        n_genes_in_slope=slope_res['n'],\n"
        "    ))\n"
        "stats_df = pd.DataFrame(rows).set_index('region').sort_values('slope')\n"
        "stats_df\n"
    ),
    md("### Persist + plot"),
    code(
        "stats_df.to_csv(TBL / 'regional_attenuation_stats.tsv', sep='\\t')\n"
        "\n"
        "fig, ax = plt.subplots(figsize=(5.5, 0.4 * len(stats_df) + 1.0))\n"
        "y = np.arange(len(stats_df))\n"
        "ax.errorbar(stats_df['slope'], y,\n"
        "            xerr=[stats_df['slope'] - stats_df['slope_ci_low'],\n"
        "                  stats_df['slope_ci_high'] - stats_df['slope']],\n"
        "            fmt='o', color='#222', ecolor='#888', capsize=2)\n"
        "ax.axvline(1, color='#bbb', linewidth=0.6, linestyle='--')\n"
        "ax.axvline(0, color='#bbb', linewidth=0.6, linestyle=':')\n"
        "ax.set_yticks(y); ax.set_yticklabels(stats_df.index, fontsize=8)\n"
        "ax.set_xlabel('attenuation slope (LFC$_{BRI-WT}$ on LFC$_{PBS-WT}$)')\n"
        "ax.set_title('1.0 = no rescue   |   0.0 = full rescue', fontsize=8,\n"
        "             loc='left', color='#555')\n"
        "for s in ('top', 'right'):\n"
        "    ax.spines[s].set_visible(False)\n"
        "fig.tight_layout()\n"
        "fig.savefig(FIG / 'regional_attenuation_slope.svg',\n"
        "            bbox_inches='tight')\n"
        "fig.savefig(FIG / 'regional_attenuation_slope.png',\n"
        "            bbox_inches='tight', dpi=200)\n"
        "plt.show()\n"
    ),
])


# -----------------------------------------------------------------------------
# 03 — signature mixed model (T4)
# -----------------------------------------------------------------------------

write("03-signature-mixed-model.ipynb", [
    md(
        "# 03/03 — Mixed-effects model on signature scores\n"
        "\n"
        "Per-spot signature scores (PIG, microglial activation, "
        "endocytosis) tested with random effects for mouse and section. "
        "One p-value per signature per region instead of thousands.\n"
        "\n"
        "    score ~ treatment + region + treatment:region + (1|mouse) + (1|section)\n"
        "\n"
        "Uses `statsmodels.MixedLM`. The relevant contrast is\n"
        "`treatment[BRICHOS] - treatment[PBS]` per region."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "import statsmodels.formula.api as smf\n"
        "from utils.attenuation import responder_zshift\n"
        "\n"
        "adata = sc.read_h5ad(H5AD)\n"
        "obs = adata.obs.copy()\n"
        "obs['mouse']   = obs[SAMPLE_KEY]\n"
        "obs['section'] = obs[SAMPLE_KEY]   # adjust if section != sample\n"
        "obs['region']  = obs[REGION_KEY]\n"
        "obs['treatment'] = obs[TREATMENT_KEY]\n"
    ),
    md(
        "### Compute signature scores if not already present\n"
        "\n"
        "Replace `PIG_GENES` etc. with the gene lists you actually use. "
        "`sc.tl.score_genes` adds a column to `adata.obs` and assumes "
        "**log-normalised data in `.X`**. Do *not* overwrite `.X` with "
        "counts before running this — only do that for pseudobulk."
    ),
    code(
        "PIG_GENES = ['Cst7', 'Itgax', 'Apoe', 'Trem2', 'Tyrobp', 'Lpl',\n"
        "             'Csf1', 'Spp1', 'Clec7a', 'Ccl3', 'Ctsb', 'Ctsd']\n"
        "MICROGLIA_GENES = ['Iba1', 'Aif1', 'Cx3cr1', 'P2ry12', 'Tmem119']\n"
        "ENDOCYTOSIS_GENES = ['Rab5a', 'Rab7a', 'Eea1', 'Lamp1', 'Lamp2']\n"
        "\n"
        "for name, genes in [('PIG', PIG_GENES),\n"
        "                    ('microglia', MICROGLIA_GENES),\n"
        "                    ('endocytosis', ENDOCYTOSIS_GENES)]:\n"
        "    present = [g for g in genes if g in adata.var_names]\n"
        "    if not present:\n"
        "        print(f'!! {name}: no genes match var_names; skipping')\n"
        "        continue\n"
        "    sc.tl.score_genes(adata, present, score_name=f'{name}_score',\n"
        "                      use_raw=False)\n"
        "    obs[f'{name}_score'] = adata.obs[f'{name}_score'].values\n"
    ),
    md("### Mixed model per signature"),
    code(
        "results = []\n"
        "for sig in ['PIG_score', 'microglia_score', 'endocytosis_score']:\n"
        "    if sig not in obs.columns:\n"
        "        continue\n"
        "    df = obs[[sig, 'treatment', 'region', 'mouse']].dropna()\n"
        "    df = df.rename(columns={sig: 'score'})\n"
        "    md_ = smf.mixedlm('score ~ C(treatment, Treatment(\"WT\"))'\n"
        "                      ' * C(region)',\n"
        "                      data=df, groups=df['mouse'])\n"
        "    fit = md_.fit(method='lbfgs', reml=False, maxiter=200)\n"
        "    res = fit.summary().tables[1]\n"
        "    res = pd.DataFrame(res.data[1:], columns=res.data[0])\n"
        "    res.insert(0, 'signature', sig)\n"
        "    results.append(res)\n"
        "\n"
        "if results:\n"
        "    mixed_df = pd.concat(results, ignore_index=True)\n"
        "    mixed_df.to_csv(TBL / 'signature_mixed_model.tsv', sep='\\t', index=False)\n"
        "    mixed_df.head(20)\n"
    ),
    md(
        "### Notes\n"
        "\n"
        "* If convergence warnings appear, try `reml=True` and "
        "`method='powell'`.\n"
        "* For section-level random effect, add `re_formula='1'` and a "
        "`vc_formula={'section': '0 + C(section)'}`. statsmodels can be "
        "fussy; if it fights you, drop section and keep only mouse.\n"
        "* Confounders (slide batch, sex, age) go in the fixed-effect "
        "formula — `+ C(slide_batch) + sex + age`."
    ),
])


# -----------------------------------------------------------------------------
# 04 — responder spots
# -----------------------------------------------------------------------------

write("04-responder-spots.ipynb", [
    md(
        "# 03/04 — Responder spot map\n"
        "\n"
        "For each spot, compute z-shift of a pathology signature relative "
        "to that region's WT distribution. Spots with strongly negative "
        "z (closer to WT than expected) under BRICHOS = *responder spots*. "
        "Map to tissue, ask whether responders cluster regionally or "
        "follow a meningeal-to-parenchymal axis.\n"
        "\n"
        "Output: an updated AnnData with a `PIG_zshift` obs column and a "
        "`responder` boolean."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from utils.attenuation import responder_zshift\n"
        "\n"
        "adata = sc.read_h5ad(H5AD)\n"
        "if 'PIG_score' not in adata.obs.columns:\n"
        "    PIG_GENES = ['Cst7', 'Itgax', 'Apoe', 'Trem2', 'Tyrobp', 'Lpl',\n"
        "                 'Csf1', 'Spp1', 'Clec7a', 'Ccl3', 'Ctsb', 'Ctsd']\n"
        "    present = [g for g in PIG_GENES if g in adata.var_names]\n"
        "    sc.tl.score_genes(adata, present, score_name='PIG_score',\n"
        "                      use_raw=False)\n"
    ),
    md("### z-shift relative to regional WT"),
    code(
        "z = responder_zshift(adata,\n"
        "                     score_key='PIG_score',\n"
        "                     treatment_key=TREATMENT_KEY,\n"
        "                     region_key=REGION_KEY,\n"
        "                     sample_key=SAMPLE_KEY,\n"
        "                     ref='WT')\n"
        "adata.obs['PIG_zshift'] = z\n"
        "adata.obs['responder']  = adata.obs['PIG_zshift'] < -1.0\n"
        "adata.obs.groupby([TREATMENT_KEY, REGION_KEY])['responder']\\\n"
        "         .mean().unstack(fill_value=0).round(3)\n"
    ),
    md("### Spatial plots — per-treatment composite"),
    code(
        "import matplotlib.pyplot as plt\n"
        "for treat in ['WT', 'PBS', 'BRICHOS']:\n"
        "    sub = adata[adata.obs[TREATMENT_KEY] == treat]\n"
        "    if sub.n_obs == 0:\n"
        "        continue\n"
        "    libs = sub.obs[SAMPLE_KEY].unique()\n"
        "    fig, axes = plt.subplots(1, len(libs),\n"
        "                             figsize=(3 * len(libs), 3.2),\n"
        "                             squeeze=False)\n"
        "    for ax, lib in zip(axes.flat, libs):\n"
        "        sub_l = sub[sub.obs[SAMPLE_KEY] == lib]\n"
        "        sp = adata.uns['spatial'][lib]\n"
        "        sf = sp['scalefactors']['tissue_hires_scalef']\n"
        "        ax.imshow(sp['images']['hires'], origin='upper')\n"
        "        xy = sub_l.obsm['spatial'] * sf\n"
        "        c = sub_l.obs['PIG_zshift'].clip(-3, 3).values\n"
        "        ax.scatter(xy[:, 0], xy[:, 1], c=c, s=0.6,\n"
        "                   cmap='coolwarm', vmin=-3, vmax=3, edgecolors='none')\n"
        "        ax.set_title(f'{lib} [{treat}]', fontsize=8, loc='left')\n"
        "        ax.set_xticks([]); ax.set_yticks([])\n"
        "        for s in ax.spines.values():\n"
        "            s.set_visible(False)\n"
        "    fig.suptitle(f'{treat}: PIG z-shift relative to regional WT',\n"
        "                 fontsize=10)\n"
        "    fig.tight_layout()\n"
        "    fig.savefig(FIG / f'responders_{treat}.svg', bbox_inches='tight')\n"
        "    fig.savefig(FIG / f'responders_{treat}.png',\n"
        "                bbox_inches='tight', dpi=200)\n"
        "    plt.show()\n"
    ),
    md(
        "### Quantify: responder fraction per region\n"
        "\n"
        "Compare BRICHOS vs PBS responder fraction with a per-mouse "
        "Mann-Whitney for robustness."
    ),
    code(
        "from scipy.stats import mannwhitneyu\n"
        "rows = []\n"
        "for region in adata.obs[REGION_KEY].dropna().unique():\n"
        "    by_mouse = (adata.obs[adata.obs[REGION_KEY] == region]\n"
        "                .groupby(SAMPLE_KEY)\n"
        "                .agg(treatment=(TREATMENT_KEY, 'first'),\n"
        "                     responder_frac=('responder', 'mean')))\n"
        "    pbs = by_mouse.loc[by_mouse['treatment'] == 'PBS',\n"
        "                       'responder_frac'].values\n"
        "    bri = by_mouse.loc[by_mouse['treatment'] == 'BRICHOS',\n"
        "                       'responder_frac'].values\n"
        "    if len(pbs) >= 2 and len(bri) >= 2:\n"
        "        u, p = mannwhitneyu(bri, pbs, alternative='greater')\n"
        "    else:\n"
        "        u, p = (np.nan, np.nan)\n"
        "    rows.append(dict(region=region, mean_pbs=np.mean(pbs) if len(pbs) else np.nan,\n"
        "                     mean_bri=np.mean(bri) if len(bri) else np.nan,\n"
        "                     mwu_p=p))\n"
        "resp_df = pd.DataFrame(rows).set_index('region').sort_values('mwu_p')\n"
        "resp_df.to_csv(TBL / 'responder_fraction_by_region.tsv', sep='\\t')\n"
        "resp_df\n"
    ),
])


# -----------------------------------------------------------------------------
# 05 — permutation null
# -----------------------------------------------------------------------------

write("05-permutation-null.ipynb", [
    md(
        "# 03/05 — Permutation null for attenuation slope\n"
        "\n"
        "Shuffle treatment labels at the *mouse* level and recompute the "
        "regional attenuation slope. Empirical p per region vs the "
        "observed slope. With 3 vs 3 mice the permutation space is small "
        "(~20 unique permutations) so this is a sanity check, not a "
        "high-power test — but it falsifies obviously spurious slopes."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from itertools import permutations\n"
        "from utils.attenuation import pseudobulk_lfc, attenuation_slope\n"
        "\n"
        "counts = pd.read_csv(TBL / 'pseudobulk_counts.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "meta   = pd.read_csv(TBL / 'pseudobulk_meta.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "\n"
        "obs_stats = pd.read_csv(TBL / 'regional_attenuation_stats.tsv',\n"
        "                        sep='\\t', index_col=0)\n"
    ),
    code(
        "def slope_for(meta_local, region):\n"
        "    pbs_v_wt = pseudobulk_lfc(counts, meta_local,\n"
        "                              group_a='PBS', group_b='WT',\n"
        "                              region=region)\n"
        "    bri_v_wt = pseudobulk_lfc(counts, meta_local,\n"
        "                              group_a='BRICHOS', group_b='WT',\n"
        "                              region=region)\n"
        "    common = pbs_v_wt.lfc.index.intersection(bri_v_wt.lfc.index)\n"
        "    sig = pbs_v_wt.padj.loc[common] < 0.05\n"
        "    return attenuation_slope(pbs_v_wt.lfc.loc[common],\n"
        "                             bri_v_wt.lfc.loc[common],\n"
        "                             sig, n_boot=0)['slope']\n"
        "\n"
        "mice = meta[['sample', 'treatment']].drop_duplicates()\n"
        "labels = mice['treatment'].values\n"
        "samples = mice['sample'].values\n"
        "\n"
        "perms = list(set(permutations(labels)))\n"
        "print(f'{len(perms)} unique permutations')\n"
    ),
    code(
        "rng_perms = perms[:200]   # cap for speed\n"
        "rows = []\n"
        "for region in obs_stats.index:\n"
        "    obs_slope = obs_stats.loc[region, 'slope']\n"
        "    null_slopes = []\n"
        "    for perm in rng_perms:\n"
        "        local = meta.copy()\n"
        "        mapping = dict(zip(samples, perm))\n"
        "        local['treatment'] = local['sample'].map(mapping)\n"
        "        try:\n"
        "            null_slopes.append(slope_for(local, region))\n"
        "        except Exception:\n"
        "            continue\n"
        "    null_slopes = np.array(null_slopes)\n"
        "    p_left = (null_slopes <= obs_slope).mean()  # rescue test\n"
        "    rows.append(dict(region=region, obs_slope=obs_slope,\n"
        "                     null_median=np.nanmedian(null_slopes),\n"
        "                     p_one_sided=p_left,\n"
        "                     n_perms=len(null_slopes)))\n"
        "perm_df = pd.DataFrame(rows).set_index('region')\n"
        "perm_df.to_csv(TBL / 'permutation_null.tsv', sep='\\t')\n"
        "perm_df\n"
    ),
])


# -----------------------------------------------------------------------------
# 06 — confounder + composition robustness
# -----------------------------------------------------------------------------

write("06-confounder-checks.ipynb", [
    md(
        "# 03/06 — Confounder + composition robustness\n"
        "\n"
        "Sanity layers reviewers will ask about:\n"
        "\n"
        "1. **Section QC** — drop low-quality sections, re-run T1/T3.\n"
        "2. **Slide batch** — re-fit signature mixed model with batch as "
        "fixed effect; check estimates do not collapse.\n"
        "3. **Cell composition** — regress out tangram cell-type "
        "fractions per spot before computing signature scores; rerun "
        "responder map.\n"
        "4. **Meninges contamination** — compare attenuation in "
        "hippocampus with vs without spots within N pixels of the rim."
    ),
    code(COMMON_HEADER_CODE),
    md("### 1. Section QC drop-and-rerun"),
    code(
        "adata = sc.read_h5ad(H5AD)\n"
        "qc = (adata.obs.groupby(SAMPLE_KEY)\n"
        "      .agg(n_spots=('total_counts','size'),\n"
        "           median_counts=('total_counts','median'),\n"
        "           pct_mito=('pct_counts_mt','mean'))\n"
        "      .sort_values('median_counts'))\n"
        "qc\n"
    ),
    code(
        "BAD_SAMPLES = []  # fill in based on QC: e.g. ['P28052_103']\n"
        "if BAD_SAMPLES:\n"
        "    keep = ~adata.obs[SAMPLE_KEY].isin(BAD_SAMPLES)\n"
        "    adata_qc = adata[keep].copy()\n"
        "else:\n"
        "    adata_qc = adata\n"
        "print(adata_qc)\n"
    ),
    md("### 2. Refit attenuation on QC-filtered cohort"),
    code(
        "from utils.attenuation import (\n"
        "    make_pseudobulk, pseudobulk_lfc, attenuation_slope,\n"
        "    sign_concordance,\n"
        ")\n"
        "counts_qc, meta_qc = make_pseudobulk(adata_qc,\n"
        "    sample_key=SAMPLE_KEY, region_key=REGION_KEY,\n"
        "    treatment_key=TREATMENT_KEY, min_spots=20,\n"
        "    layer=COUNT_LAYER)\n"
        "counts_qc = counts_qc.round().astype(int)\n"
        "rows = []\n"
        "for region in sorted(meta_qc['region'].dropna().unique()):\n"
        "    try:\n"
        "        pw = pseudobulk_lfc(counts_qc, meta_qc, 'PBS', 'WT', region)\n"
        "        bw = pseudobulk_lfc(counts_qc, meta_qc, 'BRICHOS', 'WT', region)\n"
        "    except ValueError:\n"
        "        continue\n"
        "    common = pw.lfc.index.intersection(bw.lfc.index)\n"
        "    sig = pw.padj.loc[common] < 0.05\n"
        "    res = attenuation_slope(pw.lfc.loc[common],\n"
        "                            bw.lfc.loc[common], sig)\n"
        "    rows.append(dict(region=region, slope_qc=res['slope'],\n"
        "                     ci_low=res['ci_low'], ci_high=res['ci_high']))\n"
        "qc_df = pd.DataFrame(rows).set_index('region')\n"
        "qc_df.to_csv(TBL / 'qc_filtered_attenuation.tsv', sep='\\t')\n"
        "qc_df\n"
    ),
    md(
        "### 3. Composition adjustment\n"
        "\n"
        "Requires a `cell_type_fractions` obsm slot from tangram (or "
        "another deconvolution). For each signature score, regress out "
        "fractions and rerun the mixed model on residuals.\n"
        "\n"
        "Skeleton below — fill in once `obsm['cell_type_fractions']` is "
        "populated."
    ),
    code(
        "if 'cell_type_fractions' not in adata.obsm:\n"
        "    print('skip — no obsm[cell_type_fractions] yet')\n"
        "else:\n"
        "    from sklearn.linear_model import LinearRegression\n"
        "    X = adata.obsm['cell_type_fractions']\n"
        "    for sig in ['PIG_score', 'microglia_score']:\n"
        "        if sig not in adata.obs.columns:\n"
        "            continue\n"
        "        y = adata.obs[sig].values\n"
        "        keep = np.isfinite(y).all() if y.ndim > 1 else np.isfinite(y)\n"
        "        lr = LinearRegression().fit(X[keep], y[keep])\n"
        "        adata.obs[f'{sig}_resid'] = y - lr.predict(X)\n"
        "        print(f'{sig}: residualised against cell types')\n"
    ),
    md(
        "### 4. Meninges-rim sensitivity\n"
        "\n"
        "Recompute hippocampal attenuation slope with and without spots "
        "within `D` pixels of the rim mask. If the slope is sensitive "
        "to D, the regional rescue may be partly meningeal contamination."
    ),
    code(
        "# Pseudo-code: assumes obs has 'distance_to_rim_px' from nb 07.\n"
        "if 'distance_to_rim_px' not in adata.obs.columns:\n"
        "    print('skip — distance_to_rim_px not yet computed')\n"
        "else:\n"
        "    rows = []\n"
        "    for d in [0, 50, 100, 200, 400]:\n"
        "        mask = adata.obs['distance_to_rim_px'] >= d\n"
        "        a = adata[mask].copy()\n"
        "        c, m = make_pseudobulk(a, sample_key=SAMPLE_KEY,\n"
        "                                region_key=REGION_KEY,\n"
        "                                treatment_key=TREATMENT_KEY,\n"
        "                                min_spots=10,\n"
        "                                layer=COUNT_LAYER)\n"
        "        c = c.round().astype(int)\n"
        "        try:\n"
        "            pw = pseudobulk_lfc(c, m, 'PBS', 'WT',\n"
        "                                 region='Hippocampal_formation')\n"
        "            bw = pseudobulk_lfc(c, m, 'BRICHOS', 'WT',\n"
        "                                 region='Hippocampal_formation')\n"
        "        except (KeyError, ValueError):\n"
        "            continue\n"
        "        common = pw.lfc.index.intersection(bw.lfc.index)\n"
        "        sig = pw.padj.loc[common] < 0.05\n"
        "        res = attenuation_slope(pw.lfc.loc[common],\n"
        "                                bw.lfc.loc[common], sig)\n"
        "        rows.append(dict(min_distance_px=d, slope=res['slope'],\n"
        "                         n=res['n']))\n"
        "    pd.DataFrame(rows)\n"
    ),
])


# -----------------------------------------------------------------------------
# README
# -----------------------------------------------------------------------------

readme = OUT_DIR / "README.md"
readme.write_text(
    "# notebooks/03 — Statistical attenuation analysis\n"
    "\n"
    "Run order:\n"
    "\n"
    "1. `01-pseudobulk.ipynb` — collapse spot counts to (mouse, region).\n"
    "2. `02-attenuation-stats.ipynb` — T1 (LFC), T2 (sign concordance), "
    "T3 (slope + bootstrap CI). Produces `regional_attenuation_stats.tsv`.\n"
    "3. `03-signature-mixed-model.ipynb` — T4 mixed-effects on PIG / "
    "microglia / endocytosis scores.\n"
    "4. `04-responder-spots.ipynb` — per-spot z-shift, responder maps, "
    "regional fraction test.\n"
    "5. `05-permutation-null.ipynb` — mouse-level label permutation, "
    "empirical p per region.\n"
    "6. `06-confounder-checks.ipynb` — QC filter, slide-batch, "
    "composition adjustment, meninges-rim sensitivity.\n"
    "\n"
    "Outputs land in `results/tables/attenuation/` and "
    "`results/figures/manuscript/`.\n"
    "\n"
    "Top-of-notebook variables to set on first run:\n"
    "\n"
    "* `SAMPLE_KEY` — obs column for mouse/section id (try `sample_id` "
    "or `library_id`).\n"
    "* `REGION_KEY` — obs column for region label "
    "(`anatomical_region`, `region`, etc.).\n"
    "* `TREATMENT_KEY` — usually `treatment` with values WT / PBS / BRICHOS.\n"
)
print(f"wrote {readme}")
