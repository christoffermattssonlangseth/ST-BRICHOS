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
        "    min_spots=10,            # was 20 — relax to keep small regions\n"
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
        "# 03/02 — Disease signature + per-region rescue\n"
        "\n"
        "With only 1 WT mouse, per-region WT-anchored t-tests are not "
        "valid. Pivot:\n"
        "\n"
        "1. **Disease signature (global, once).** Whole-tissue PBS vs WT "
        "pseudobulk (n=7 vs n=1). Pooled-variance fallback in "
        "`pseudobulk_lfc(min_n=1)`. Rank disease-up / disease-down "
        "genes. Cross-check against literature PIG/DAM lists.\n"
        "2. **Per-region rescue (well-powered).** BRI vs PBS per region "
        "(n=3 vs n=2-7). Welch on log-CPM.\n"
        "3. **T2' sign concordance.** Of disease-up genes, fraction "
        "*lowered* by BRICHOS per region (binomial vs 0.5).\n"
        "4. **T3' rescue slope.** Regress per-region LFC_BRI-vs-PBS on "
        "global LFC_disease over disease-DEGs. Slope < 0 = rescue.\n"
        "\n"
        "Outputs `regional_rescue_stats.tsv`, `disease_signature.tsv`."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from utils.attenuation import (\n"
        "    pseudobulk_lfc, signed_rescue, rescue_slope\n"
        ")\n"
        "\n"
        "counts = pd.read_csv(TBL / 'pseudobulk_counts.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "meta   = pd.read_csv(TBL / 'pseudobulk_meta.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "regions = sorted(meta['region'].dropna().unique())\n"
        "print(len(regions), 'regions')\n"
    ),
    md(
        "### 1. Disease signature — whole-tissue PBS vs WT\n"
        "\n"
        "Sum across regions per (sample). With n=1 WT we use "
        "`min_n=1`: pooled variance from PBS replicates. **Significance "
        "values here are inflated** — use the LFC ranking and treat "
        "padj as descriptive."
    ),
    code(
        "global_counts = counts.copy()\n"
        "global_counts['__sample__'] = meta['sample']\n"
        "global_counts = global_counts.groupby('__sample__').sum()\n"
        "global_meta = (meta.groupby('sample')\n"
        "               .agg(treatment=('treatment','first'))\n"
        "               .reset_index().set_index('sample'))\n"
        "global_meta.index.name = None\n"
        "global_meta = global_meta.rename_axis(None)\n"
        "# emulate the (sample, region) row-id format expected by helpers\n"
        "global_counts.index = global_counts.index.astype(str)\n"
        "global_meta.index   = global_meta.index.astype(str)\n"
        "global_meta['region'] = '__all__'\n"
        "global_meta['sample'] = global_meta.index\n"
        "\n"
        "disease = pseudobulk_lfc(global_counts, global_meta,\n"
        "                         group_a='PBS', group_b='WT',\n"
        "                         region='__all__', min_n=1)\n"
        "disease_df = disease.to_frame()\n"
        "disease_df.to_csv(TBL / 'disease_signature.tsv', sep='\\t')\n"
        "print(f'  n_a={disease.n_a}, n_b={disease.n_b}')\n"
        "disease_df.sort_values('lfc', ascending=False).head()\n"
    ),
    md(
        "### Pick disease gene set\n"
        "\n"
        "Two ways, pick one or report both:\n"
        "\n"
        "* **Data-driven**: top genes by |LFC| in PBS-vs-WT (descriptive "
        "padj filter).\n"
        "* **Literature**: known PIG / DAM gene lists.\n"
        "\n"
        "Default: union, with provenance tag."
    ),
    code(
        "# Data-driven set\n"
        "lfc_thresh, padj_thresh = 0.5, 0.1\n"
        "data_up = disease_df.index[(disease_df.lfc >  lfc_thresh) &\n"
        "                            (disease_df.padj < padj_thresh)]\n"
        "data_dn = disease_df.index[(disease_df.lfc < -lfc_thresh) &\n"
        "                            (disease_df.padj < padj_thresh)]\n"
        "\n"
        "# Literature PIG (Chen et al, Cell 2020) — full curated list\n"
        "# from utils/gene_sets.py.\n"
        "from utils.gene_sets import PIG_CHEN2020, filter_to_var\n"
        "lit_up = filter_to_var(PIG_CHEN2020, disease_df.index)\n"
        "print(f'PIG genes in disease table: {len(lit_up)} / {len(PIG_CHEN2020)}')\n"
        "\n"
        "disease_up = pd.Index(set(list(data_up) + lit_up))\n"
        "disease_dn = pd.Index(data_dn)\n"
        "print(f'data-driven up: {len(data_up)}, dn: {len(data_dn)}')\n"
        "print(f'literature PIG present: {len(lit_up)}')\n"
        "print(f'union disease-up: {len(disease_up)}')\n"
    ),
    md("### 2. Per-region BRI vs PBS LFC"),
    code(
        "rescue_lfc = {}\n"
        "for region in regions:\n"
        "    try:\n"
        "        res = pseudobulk_lfc(counts, meta,\n"
        "                             group_a='BRICHOS', group_b='PBS',\n"
        "                             region=region, min_n=2)\n"
        "    except ValueError as e:\n"
        "        print(f'  {region}: skip ({e})')\n"
        "        continue\n"
        "    rescue_lfc[region] = res\n"
        "    print(f'  {region}: n_BRI={res.n_a}, n_PBS={res.n_b}, '\n"
        "          f'n_genes={len(res.lfc)}')\n"
    ),
    md("### 3. + 4. Per-region rescue stats"),
    code(
        "rows = []\n"
        "for region, res in rescue_lfc.items():\n"
        "    sc_res = signed_rescue(res.lfc, disease_up, disease_dn)\n"
        "    sl_res = rescue_slope(disease_df.lfc, res.lfc,\n"
        "                          sig_mask=list(disease_up),\n"
        "                          n_boot=1000)\n"
        "    rows.append(dict(\n"
        "        region=region,\n"
        "        n_BRI=res.n_a, n_PBS=res.n_b,\n"
        "        n_signature=sc_res['n'],\n"
        "        frac_concordant=sc_res['frac_concordant'],\n"
        "        binom_p=sc_res['binom_p'],\n"
        "        wilcoxon_p=sc_res['wilcoxon_p'],\n"
        "        median_signed_lfc=sc_res['median_signed_lfc'],\n"
        "        slope=sl_res['slope'],\n"
        "        slope_ci_low=sl_res['ci_low'],\n"
        "        slope_ci_high=sl_res['ci_high'],\n"
        "    ))\n"
        "stats_df = (pd.DataFrame(rows).set_index('region')\n"
        "            .sort_values('slope'))\n"
        "stats_df.to_csv(TBL / 'regional_rescue_stats.tsv', sep='\\t')\n"
        "stats_df\n"
    ),
    md(
        "### Plot — rescue slope per region\n"
        "\n"
        "Negative slope = treatment opposes disease. Bootstrap CI shown."
    ),
    code(
        "fig, ax = plt.subplots(figsize=(5.5, 0.4 * len(stats_df) + 1.0))\n"
        "y = np.arange(len(stats_df))\n"
        "ax.errorbar(stats_df['slope'], y,\n"
        "            xerr=[stats_df['slope'] - stats_df['slope_ci_low'],\n"
        "                  stats_df['slope_ci_high'] - stats_df['slope']],\n"
        "            fmt='o', color='#222', ecolor='#888', capsize=2)\n"
        "ax.axvline(0, color='#bbb', linewidth=0.6, linestyle=':')\n"
        "ax.axvline(-1, color='#d62728', linewidth=0.6, linestyle='--',\n"
        "           alpha=0.4)\n"
        "ax.set_yticks(y); ax.set_yticklabels(stats_df.index, fontsize=8)\n"
        "ax.set_xlabel('rescue slope (LFC$_{BRI-PBS}$ on LFC$_{disease}$)')\n"
        "ax.set_title('0 = no effect   |   -1 = full rescue', fontsize=8,\n"
        "             loc='left', color='#555')\n"
        "for s in ('top', 'right'):\n"
        "    ax.spines[s].set_visible(False)\n"
        "fig.tight_layout()\n"
        "fig.savefig(FIG / 'regional_rescue_slope.svg', bbox_inches='tight')\n"
        "fig.savefig(FIG / 'regional_rescue_slope.png',\n"
        "            bbox_inches='tight', dpi=200)\n"
        "plt.show()\n"
    ),
    md(
        "### Plot — sign concordance fraction\n"
        "\n"
        "Of disease-signature genes, fraction moved in the rescue "
        "direction by BRICHOS. 0.5 is the chance level."
    ),
    code(
        "fig, ax = plt.subplots(figsize=(5.5, 0.4 * len(stats_df) + 1.0))\n"
        "y = np.arange(len(stats_df))\n"
        "ax.barh(y, stats_df['frac_concordant'] - 0.5, left=0.5,\n"
        "        color='#1f77b4', alpha=0.7, edgecolor='#222')\n"
        "for i, p in enumerate(stats_df['binom_p'].values):\n"
        "    star = '***' if p < 1e-3 else ('**' if p < 1e-2 else\n"
        "           ('*' if p < 5e-2 else ''))\n"
        "    ax.text(stats_df['frac_concordant'].iloc[i] + 0.01, i,\n"
        "            star, fontsize=8, va='center')\n"
        "ax.axvline(0.5, color='#bbb', linewidth=0.6, linestyle=':')\n"
        "ax.set_yticks(y); ax.set_yticklabels(stats_df.index, fontsize=8)\n"
        "ax.set_xlim(0.3, 1.0)\n"
        "ax.set_xlabel('fraction of signature genes moved by BRICHOS '\n"
        "              'in rescue direction')\n"
        "for s in ('top', 'right'):\n"
        "    ax.spines[s].set_visible(False)\n"
        "fig.tight_layout()\n"
        "fig.savefig(FIG / 'regional_sign_concordance.svg',\n"
        "            bbox_inches='tight')\n"
        "fig.savefig(FIG / 'regional_sign_concordance.png',\n"
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
        "from utils.gene_sets import (\n"
        "    PIG_CHEN2020, MICROGLIA_HOMEOSTATIC, ENDOCYTOSIS,\n"
        "    DAM_KEREN_SHAUL2017, ARM_SALA_FRIGERIO2019,\n"
        "    filter_to_var,\n"
        ")\n"
        "\n"
        "signatures = {\n"
        "    'PIG':         PIG_CHEN2020,\n"
        "    'microglia':   MICROGLIA_HOMEOSTATIC,\n"
        "    'endocytosis': ENDOCYTOSIS,\n"
        "    'DAM':         DAM_KEREN_SHAUL2017,\n"
        "    'ARM':         ARM_SALA_FRIGERIO2019,\n"
        "}\n"
        "for name, genes in signatures.items():\n"
        "    present = filter_to_var(genes, adata.var_names)\n"
        "    if not present:\n"
        "        print(f'!! {name}: no genes match var_names; skipping')\n"
        "        continue\n"
        "    sc.tl.score_genes(adata, present, score_name=f'{name}_score',\n"
        "                      use_raw=False)\n"
        "    obs[f'{name}_score'] = adata.obs[f'{name}_score'].values\n"
        "    print(f'  {name}: scored {len(present)} / {len(genes)} genes')\n"
    ),
    md("### Mixed model per signature"),
    code(
        "def fit_to_df(fit, signature):\n"
        "    ci = fit.conf_int()\n"
        "    ci.columns = ['ci_low', 'ci_high']\n"
        "    out = pd.DataFrame({\n"
        "        'coef'   : fit.params,\n"
        "        'se'     : fit.bse,\n"
        "        'z'      : fit.tvalues,\n"
        "        'pvalue' : fit.pvalues,\n"
        "    }).join(ci, how='left')\n"
        "    out.index.name = 'term'\n"
        "    out = out.reset_index()\n"
        "    out.insert(0, 'signature', signature)\n"
        "    return out\n"
        "\n"
        "results = []\n"
        "for sig in ['PIG_score', 'microglia_score', 'endocytosis_score',\n"
        "            'DAM_score', 'ARM_score']:\n"
        "    if sig not in obs.columns:\n"
        "        continue\n"
        "    df = obs[[sig, 'treatment', 'region', 'mouse']].dropna()\n"
        "    df = df.rename(columns={sig: 'score'})\n"
        "    md_ = smf.mixedlm('score ~ C(treatment, Treatment(\"WT\"))'\n"
        "                      ' * C(region)',\n"
        "                      data=df, groups=df['mouse'])\n"
        "    fit = md_.fit(method='lbfgs', reml=False, maxiter=200)\n"
        "    results.append(fit_to_df(fit, sig))\n"
        "\n"
        "if results:\n"
        "    mixed_df = pd.concat(results, ignore_index=True)\n"
        "    mixed_df.to_csv(TBL / 'signature_mixed_model.tsv',\n"
        "                    sep='\\t', index=False)\n"
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
        "from utils.gene_sets import PIG_CHEN2020, filter_to_var\n"
        "\n"
        "adata = sc.read_h5ad(H5AD)\n"
        "if 'PIG_score' not in adata.obs.columns:\n"
        "    present = filter_to_var(PIG_CHEN2020, adata.var_names)\n"
        "    print(f'PIG: scoring {len(present)} / {len(PIG_CHEN2020)} genes')\n"
        "    sc.tl.score_genes(adata, present, score_name='PIG_score',\n"
        "                      use_raw=False)\n"
    ),
    md("### z-shift relative to regional WT"),
    md(
        "**WT n=1 caveat:** the regional WT distribution comes from a "
        "single mouse, so SD reflects within-mouse spot variance only. "
        "Treat z-shift maps as descriptive. The treatment-blind "
        "fallback below uses the PBS median as a reference instead "
        "and avoids the WT bottleneck."
    ),
    code(
        "z = responder_zshift(adata,\n"
        "                     score_key='PIG_score',\n"
        "                     treatment_key=TREATMENT_KEY,\n"
        "                     region_key=REGION_KEY,\n"
        "                     sample_key=SAMPLE_KEY,\n"
        "                     ref='WT')\n"
        "adata.obs['PIG_zshift_vsWT'] = z\n"
        "\n"
        "# Treatment-blind fallback: z relative to PBS median per region.\n"
        "z_pbs = responder_zshift(adata,\n"
        "                          score_key='PIG_score',\n"
        "                          treatment_key=TREATMENT_KEY,\n"
        "                          region_key=REGION_KEY,\n"
        "                          sample_key=SAMPLE_KEY,\n"
        "                          ref='PBS')\n"
        "adata.obs['PIG_zshift_vsPBS'] = z_pbs\n"
        "adata.obs['responder'] = adata.obs['PIG_zshift_vsPBS'] < -1.0\n"
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
        "        c = sub_l.obs['PIG_zshift_vsPBS'].clip(-3, 3).values\n"
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
        "# 03/05 — Permutation null for the rescue slope\n"
        "\n"
        "Permute the BRI/PBS labels at the *mouse* level (WT held out — "
        "n=1 is fixed). Recompute the rescue slope per region. With "
        "3 BRICHOS vs 7 PBS mice we have C(10,3) = 120 unique "
        "permutations, enough resolution for an empirical p down to "
        "~1/120.\n"
        "\n"
        "Output: `results/tables/attenuation/rescue_permutation_null.tsv`."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from itertools import combinations\n"
        "from utils.attenuation import pseudobulk_lfc, rescue_slope\n"
        "\n"
        "counts = pd.read_csv(TBL / 'pseudobulk_counts.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "meta   = pd.read_csv(TBL / 'pseudobulk_meta.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "disease_df = pd.read_csv(TBL / 'disease_signature.tsv',\n"
        "                          sep='\\t', index_col=0)\n"
        "obs_stats = pd.read_csv(TBL / 'regional_rescue_stats.tsv',\n"
        "                         sep='\\t', index_col=0)\n"
    ),
    code(
        "# Build the permutation set: keep WT fixed, swap BRI/PBS labels\n"
        "# at the mouse level. Each unique partition of non-WT mice into\n"
        "# (BRI, PBS) groups of original sizes is one permutation.\n"
        "treat_by_sample = (meta.groupby('sample')['treatment']\n"
        "                   .first())\n"
        "non_wt = treat_by_sample[treat_by_sample != 'WT'].index.tolist()\n"
        "n_bri = (treat_by_sample == 'BRICHOS').sum()\n"
        "perms = list(combinations(non_wt, n_bri))\n"
        "print(f'{len(perms)} unique BRI-vs-PBS reassignments')\n"
    ),
    code(
        "padj_thresh, lfc_thresh = 0.1, 0.5\n"
        "disease_up = disease_df.index[(disease_df.lfc >  lfc_thresh) &\n"
        "                               (disease_df.padj < padj_thresh)]\n"
        "sig_mask = disease_df.index.isin(disease_up)\n"
        "sig_series = pd.Series(sig_mask, index=disease_df.index)\n"
        "\n"
        "def slope_for(meta_local, region):\n"
        "    try:\n"
        "        res = pseudobulk_lfc(counts, meta_local,\n"
        "                             group_a='BRICHOS', group_b='PBS',\n"
        "                             region=region, min_n=2)\n"
        "    except ValueError:\n"
        "        return np.nan\n"
        "    return rescue_slope(disease_df.lfc, res.lfc,\n"
        "                        sig_mask=sig_series, n_boot=0)['slope']\n"
    ),
    code(
        "regions = obs_stats.index.tolist()\n"
        "rows = []\n"
        "for region in regions:\n"
        "    obs_slope = obs_stats.loc[region, 'slope']\n"
        "    null_slopes = []\n"
        "    for bri_set in perms:\n"
        "        local = meta.copy()\n"
        "        new_treat = local['treatment'].copy()\n"
        "        non_wt_mask = new_treat != 'WT'\n"
        "        new_treat[non_wt_mask] = np.where(\n"
        "            local.loc[non_wt_mask, 'sample'].isin(bri_set),\n"
        "            'BRICHOS', 'PBS')\n"
        "        local['treatment'] = new_treat\n"
        "        s = slope_for(local, region)\n"
        "        if np.isfinite(s):\n"
        "            null_slopes.append(s)\n"
        "    null_slopes = np.array(null_slopes)\n"
        "    if len(null_slopes) == 0 or not np.isfinite(obs_slope):\n"
        "        rows.append(dict(region=region, obs_slope=obs_slope,\n"
        "                         null_median=np.nan, p_one_sided=np.nan,\n"
        "                         n_perms=0))\n"
        "        continue\n"
        "    # rescue test: smaller (more negative) slope is rare under null\n"
        "    p = (null_slopes <= obs_slope).mean()\n"
        "    rows.append(dict(region=region, obs_slope=obs_slope,\n"
        "                     null_median=float(np.nanmedian(null_slopes)),\n"
        "                     p_one_sided=float(p),\n"
        "                     n_perms=int(len(null_slopes))))\n"
        "perm_df = pd.DataFrame(rows).set_index('region')\n"
        "perm_df.to_csv(TBL / 'rescue_permutation_null.tsv', sep='\\t')\n"
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
# 07 — broader signatures (non-PIG biology)
# -----------------------------------------------------------------------------

write("07-broader-signatures.ipynb", [
    md(
        "# 03/07 — Broader signatures, beyond PIG\n"
        "\n"
        "PIG is microglia-centric. BRICHOS may also modulate astrocyte "
        "reactivity, oligodendrocyte stress, barrier integrity, neuronal "
        "/ synaptic biology, and chaperone networks. This notebook:\n"
        "\n"
        "1. **Multi-signature heatmap.** For every signature in "
        "`utils.gene_sets.SIGNATURE_REGISTRY`, compute per-region "
        "rescue concordance + median signed LFC. Heatmap of regions "
        "x signatures shows where BRICHOS hits which biology.\n"
        "2. **Discovery.** Per region, top BRI-vs-PBS DEGs *outside* "
        "PIG. Annotates each with literature-set membership so novel "
        "modulators stand out.\n"
        "3. **Pathway enrichment** (optional). GSEA via `gseapy` if "
        "installed."
    ),
    code(COMMON_HEADER_CODE),
    code(
        "from utils.attenuation import pseudobulk_lfc, signed_rescue\n"
        "from utils.gene_sets import (\n"
        "    SIGNATURE_REGISTRY, PIG_CHEN2020, filter_to_var,\n"
        ")\n"
        "\n"
        "counts = pd.read_csv(TBL / 'pseudobulk_counts.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "meta   = pd.read_csv(TBL / 'pseudobulk_meta.tsv',\n"
        "                     sep='\\t', index_col=0)\n"
        "regions = sorted(meta['region'].dropna().unique())\n"
        "print(f'{len(regions)} regions, {len(SIGNATURE_REGISTRY)} signatures')\n"
    ),
    md(
        "### Per-region BRI-vs-PBS LFC (cache once)\n"
        "\n"
        "Same call as nb 02; recompute here so this notebook can run "
        "standalone."
    ),
    code(
        "rescue_lfc = {}\n"
        "for region in regions:\n"
        "    try:\n"
        "        res = pseudobulk_lfc(counts, meta,\n"
        "                             group_a='BRICHOS', group_b='PBS',\n"
        "                             region=region, min_n=2)\n"
        "    except ValueError as e:\n"
        "        print(f'  {region}: skip ({e})')\n"
        "        continue\n"
        "    rescue_lfc[region] = res\n"
        "print(f'usable regions: {len(rescue_lfc)}')\n"
    ),
    md(
        "### 1. Multi-signature heatmap\n"
        "\n"
        "For every (region, signature) pair, signed median LFC of "
        "signature genes under BRICHOS, with disease direction taken "
        "from the global PBS-vs-WT signature (`disease_signature.tsv`)."
    ),
    code(
        "from utils.gene_sets import SIGNATURE_DIRECTION\n"
        "\n"
        "disease_df = pd.read_csv(TBL / 'disease_signature.tsv',\n"
        "                          sep='\\t', index_col=0)\n"
        "lfc_thresh, padj_thresh = 0.5, 0.1\n"
        "disease_up = set(disease_df.index[(disease_df.lfc >  lfc_thresh) &\n"
        "                                   (disease_df.padj < padj_thresh)])\n"
        "disease_dn = set(disease_df.index[(disease_df.lfc < -lfc_thresh) &\n"
        "                                   (disease_df.padj < padj_thresh)])\n"
        "\n"
        "def rescue_for_signature(lfc, sig_genes, prior_direction):\n"
        "    \"\"\"Score rescue (positive=BRICHOS counters disease).\n"
        "\n"
        "    Each gene's expected sign comes from the data-driven disease\n"
        "    set first; if the gene is not in it, the signature-level\n"
        "    `prior_direction` is used. Genes with no information are\n"
        "    skipped.\n"
        "    \"\"\"\n"
        "    sig = pd.Index(sig_genes).intersection(lfc.index)\n"
        "    if len(sig) == 0:\n"
        "        return dict(median_signed=np.nan, frac_concordant=np.nan,\n"
        "                    n=0)\n"
        "    s = lfc.loc[sig].copy()\n"
        "    sign = pd.Series(np.nan, index=sig)\n"
        "    # data-driven (per-gene) sign overrides prior\n"
        "    sign.loc[sig.intersection(disease_up)] = -1.0\n"
        "    sign.loc[sig.intersection(disease_dn)] = +1.0\n"
        "    # signature-level prior fills the rest\n"
        "    if prior_direction == 'up':\n"
        "        sign = sign.fillna(-1.0)\n"
        "    elif prior_direction == 'down':\n"
        "        sign = sign.fillna(+1.0)\n"
        "    s = s.loc[sign.dropna().index]\n"
        "    sign = sign.dropna()\n"
        "    if len(s) == 0:\n"
        "        return dict(median_signed=np.nan, frac_concordant=np.nan,\n"
        "                    n=0)\n"
        "    signed = s * sign  # positive = rescue direction\n"
        "    return dict(median_signed=float(signed.median()),\n"
        "                frac_concordant=float((signed > 0).mean()),\n"
        "                n=int(len(signed)))\n"
        "\n"
        "rows = []\n"
        "for region, res in rescue_lfc.items():\n"
        "    for sig_name, sig_genes in SIGNATURE_REGISTRY.items():\n"
        "        prior = SIGNATURE_DIRECTION.get(sig_name)\n"
        "        out = rescue_for_signature(res.lfc, sig_genes, prior)\n"
        "        rows.append(dict(region=region, signature=sig_name,\n"
        "                         prior=prior, **out))\n"
        "broad_df = pd.DataFrame(rows)\n"
        "broad_df.to_csv(TBL / 'broader_signature_rescue.tsv',\n"
        "                sep='\\t', index=False)\n"
        "broad_df.head()\n"
    ),
    code(
        "# Order signatures so disease-up programs sit on top, homeostatic\n"
        "# (disease-down) below; label arrows indicate the prior used.\n"
        "order_up = [s for s, d in SIGNATURE_DIRECTION.items()\n"
        "            if d == 'up' and s in SIGNATURE_REGISTRY]\n"
        "order_dn = [s for s, d in SIGNATURE_DIRECTION.items()\n"
        "            if d == 'down' and s in SIGNATURE_REGISTRY]\n"
        "ordered = order_up + order_dn\n"
        "\n"
        "pivot = (broad_df.pivot(index='signature', columns='region',\n"
        "                         values='median_signed')\n"
        "         .reindex(ordered).fillna(0))\n"
        "ylabels = [\n"
        "    ('\\u2191 ' if SIGNATURE_DIRECTION.get(s) == 'up' else '\\u2193 ') + s\n"
        "    for s in pivot.index\n"
        "]\n"
        "\n"
        "fig, ax = plt.subplots(figsize=(0.55 * pivot.shape[1] + 2.8,\n"
        "                                 0.32 * pivot.shape[0] + 1.0))\n"
        "vmax = max(0.3, np.nanpercentile(np.abs(pivot.values), 95))\n"
        "im = ax.imshow(pivot.values, cmap='RdBu_r',\n"
        "               vmin=-vmax, vmax=vmax, aspect='auto')\n"
        "ax.set_xticks(range(pivot.shape[1]))\n"
        "ax.set_xticklabels(pivot.columns, rotation=40, ha='right',\n"
        "                   fontsize=7)\n"
        "ax.set_yticks(range(pivot.shape[0]))\n"
        "ax.set_yticklabels(ylabels, fontsize=7)\n"
        "# divider between disease-up and disease-down blocks\n"
        "ax.axhline(len(order_up) - 0.5, color='#222', linewidth=0.6)\n"
        "ax.text(-0.5, len(order_up) / 2 - 0.5, 'disease-up',\n"
        "        rotation=90, ha='right', va='center', fontsize=6,\n"
        "        color='#666')\n"
        "ax.text(-0.5, len(order_up) + len(order_dn) / 2 - 0.5,\n"
        "        'homeostatic / disease-down',\n"
        "        rotation=90, ha='right', va='center', fontsize=6,\n"
        "        color='#666')\n"
        "for s in ax.spines.values():\n"
        "    s.set_visible(False)\n"
        "cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)\n"
        "cb.set_label('median signed LFC (rescue+)', fontsize=7)\n"
        "ax.set_title('Where BRICHOS hits which biology', loc='left',\n"
        "             fontsize=8, pad=4)\n"
        "fig.tight_layout()\n"
        "fig.savefig(FIG / 'broader_signatures_heatmap.svg',\n"
        "            bbox_inches='tight')\n"
        "fig.savefig(FIG / 'broader_signatures_heatmap.png',\n"
        "            bbox_inches='tight', dpi=200)\n"
        "plt.show()\n"
    ),
    md(
        "### 2. Per-region top BRI-vs-PBS DEGs OUTSIDE PIG\n"
        "\n"
        "Discovery mode: which novel genes (not in PIG) are strongly "
        "modulated by BRICHOS? Each region gets a top-N by |LFC| with "
        "literature-set membership annotated."
    ),
    code(
        "TOPN = 30\n"
        "pig_set = set(PIG_CHEN2020)\n"
        "\n"
        "def annotate_membership(gene):\n"
        "    return ','.join(name for name, gs in SIGNATURE_REGISTRY.items()\n"
        "                    if gene in gs)\n"
        "\n"
        "discovery_tables = {}\n"
        "for region, res in rescue_lfc.items():\n"
        "    df = pd.DataFrame({'lfc': res.lfc, 'pvalue': res.pvalue,\n"
        "                       'padj': res.padj}).copy()\n"
        "    df = df[df.padj < 0.1]                  # adjust as desired\n"
        "    df = df[~df.index.isin(pig_set)]        # PIG-orthogonal\n"
        "    df['memberships'] = [annotate_membership(g) for g in df.index]\n"
        "    df = df.reindex(df.lfc.abs().sort_values(ascending=False).index)\n"
        "    discovery_tables[region] = df.head(TOPN)\n"
        "\n"
        "# write per-region tables to a single multi-sheet TSV-like file\n"
        "with open(TBL / 'discovery_top_DEGs_per_region.tsv', 'w') as f:\n"
        "    for region, sub in discovery_tables.items():\n"
        "        f.write(f'## {region}\\n')\n"
        "        sub.to_csv(f, sep='\\t')\n"
        "        f.write('\\n')\n"
        "print(f'wrote discovery_top_DEGs_per_region.tsv with '\n"
        "      f'{len(discovery_tables)} regions')\n"
    ),
    code(
        "# Display top-15 for each region inline\n"
        "for region, sub in discovery_tables.items():\n"
        "    print(f'\\n=== {region} ===')\n"
        "    print(sub.head(15).round(3).to_string())\n"
    ),
    md(
        "### 3. Pathway enrichment (optional)\n"
        "\n"
        "Requires `gseapy` (`pip install gseapy`). For each rescue region "
        "ranks genes by signed LFC and runs GSEA prerank against MSigDB "
        "Hallmark + GO_BP. Skip if package missing."
    ),
    code(
        "try:\n"
        "    import gseapy\n"
        "    have_gseapy = True\n"
        "except ImportError:\n"
        "    print('gseapy not installed — skipping enrichment')\n"
        "    have_gseapy = False\n"
        "\n"
        "if have_gseapy:\n"
        "    enrich_results = {}\n"
        "    for region, res in rescue_lfc.items():\n"
        "        rnk = res.lfc.dropna().sort_values(ascending=False)\n"
        "        try:\n"
        "            pre = gseapy.prerank(\n"
        "                rnk=rnk.reset_index().rename(\n"
        "                    columns={'index':'gene', 'lfc':'lfc'}),\n"
        "                gene_sets='GO_Biological_Process_2021',\n"
        "                organism='mouse', no_plot=True,\n"
        "                outdir=None, threads=4)\n"
        "            enrich_results[region] = pre.res2d.head(20)\n"
        "        except Exception as e:\n"
        "            print(f'  {region}: {e}')\n"
        "    for region, df in enrich_results.items():\n"
        "        df.to_csv(TBL / f'gsea_{region.replace(\" \", \"_\")}.tsv',\n"
        "                  sep='\\t', index=False)\n"
        "        print(f'\\n=== {region} top 5 ===')\n"
        "        print(df[['Term', 'NES', 'NOM p-val']].head(5).to_string())\n"
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
    "1. `01-pseudobulk.ipynb` — collapse spot counts to (sample, region) "
    "(min_spots=10 to keep small regions).\n"
    "2. `02-attenuation-stats.ipynb` — global PBS-vs-WT disease "
    "signature (n=1 WT pooled-variance fallback), then per-region "
    "BRI-vs-PBS rescue tests: sign concordance + rescue slope with "
    "bootstrap CI. Outputs `disease_signature.tsv`, "
    "`regional_rescue_stats.tsv`.\n"
    "3. `03-signature-mixed-model.ipynb` — mixed-effects on PIG / "
    "microglia / endocytosis scores at spot level (mouse random "
    "effect, robust to n=1 WT).\n"
    "4. `04-responder-spots.ipynb` — per-spot z-shift vs PBS median "
    "(treatment-blind, avoids n=1 WT bottleneck), responder maps, "
    "Mann-Whitney on per-mouse responder fraction.\n"
    "5. `05-permutation-null.ipynb` — mouse-level permutation of "
    "BRI/PBS labels (WT held out). C(10,3)=120 unique permutations, "
    "empirical p per region for the rescue slope.\n"
    "6. `06-confounder-checks.ipynb` — QC filter, slide-batch, "
    "composition adjustment, meninges-rim sensitivity.\n"
    "7. `07-broader-signatures.ipynb` — beyond PIG: rescue heatmap "
    "across ~20 cell-type and pathway signatures (astrocyte, oligo, "
    "barrier, synapse, heat-shock, etc.); discovery of top non-PIG "
    "DEGs per region; optional GSEA via gseapy.\n"
    "\n"
    "Design note — the WT cohort is a single mouse, which blocks "
    "per-region WT-anchored tests. The pipeline anchors disease once "
    "globally and then tests rescue per region using BRI vs PBS, "
    "where replication is fine.\n"
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
