#!/usr/bin/env python
"""Second-pass cell2location deconvolution with a PRUNED reference.

Pass 1 (notebook 03/09) over-assigned the immune compartment (~24 %/spot)
because cell2location spread mass across 12 similar rare immune signatures.
Here we keep only the biologically plausible CNS-resident *innate/myeloid*
identities and drop the adaptive lymphoid states (CD4 Th17, CD8 TRM, MAIT,
NK, Treg, TEx, Plasmablasts) that behaved as signature bleed.

Reuses the cached combined reference so we skip the Linnarsson+atlas rebuild.
Writes to *_pruned paths so the first-pass results stay intact.

Run headless on the peer:
    python utils/run_deconv_pruned.py
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore', category=FutureWarning)

def log(*a):
    print(*a, flush=True)

# --- project root ------------------------------------------------------
ROOT = Path.cwd().resolve()
while not (ROOT / 'utils').exists():
    if ROOT.parent == ROOT:
        raise RuntimeError('could not locate project root')
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))
from utils import deconvolution as dv  # noqa: E402

ACCELERATOR = 'cpu'

# --- paths -------------------------------------------------------------
BASEDIR = Path('/Volumes/processing2/ST_BRICHOS/data')
H5AD_ORIENTED = BASEDIR / 'ST_BRICHOS_region_subcluster_oriented.h5ad'
H5AD_BASE     = BASEDIR / 'ST_BRICHOS_region_subcluster.h5ad'
H5AD = H5AD_ORIENTED if H5AD_ORIENTED.exists() else H5AD_BASE
COUNT_LAYER = 'counts'

REF_DIR_P = ROOT / 'data' / 'external' / 'deconv_reference_pruned'
REF_DIR_P.mkdir(parents=True, exist_ok=True)
COMBINED_H5AD = ROOT / 'data' / 'external' / 'deconv_reference' / 'combined_reference.h5ad'
C2L_REF_DIR = REF_DIR_P / 'c2l_regression'
C2L_VIS_DIR = REF_DIR_P / 'c2l_spatial'

SAMPLE_KEY, REGION_KEY, TREATMENT_KEY = 'sample_id', 're_annotation_regions', 'treatment'

TBL = ROOT / 'results' / 'tables' / 'deconvolution_pruned'; TBL.mkdir(parents=True, exist_ok=True)
FIG = ROOT / 'results' / 'figures' / 'manuscript'; FIG.mkdir(parents=True, exist_ok=True)
FIG_PREFIX = 'deconv_pruned_'

# --- which cell types to keep -----------------------------------------
KEEP_MYELOID = {'MYE:Microglia', 'MYE:CAM', 'MYE:MDM', 'MYE:cDC1', 'MYE:Granulocyte'}
# brain structural types kept automatically (everything without MYE: prefix)

# ======================================================================
# 1 · load cached combined reference and prune
# ======================================================================
log('loading combined reference:', COMBINED_H5AD)
ref = sc.read_h5ad(COMBINED_H5AD)
log('  full reference:', ref.shape, '|', ref.obs.cell_type.nunique(), 'types')

ct = ref.obs['cell_type'].astype(str)
keep_mask = (~ct.str.startswith('MYE:')) | ct.isin(KEEP_MYELOID)
dropped = sorted(ct[~keep_mask].unique().tolist())
ref = ref[keep_mask.values].copy()
ref.obs['cell_type'] = ref.obs['cell_type'].astype(str).astype('category')
log('  dropped immune states:', dropped)
log('  pruned reference:', ref.shape, '|', ref.obs.cell_type.nunique(), 'types')
log(ref.obs['cell_type'].value_counts().to_string())

# ======================================================================
# 2 · reference regression -> per-gene signatures
# ======================================================================
import cell2location  # noqa: E402
from cell2location.models import RegressionModel  # noqa: E402
from cell2location.utils.filtering import filter_genes  # noqa: E402

ref.X = ref.layers['counts'].copy()
selected = filter_genes(ref, cell_count_cutoff=5,
                        cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
ref = ref[:, selected].copy()
log('genes after filter:', ref.n_vars)

RegressionModel.setup_anndata(ref, batch_key='ref_source', labels_key='cell_type')
reg = RegressionModel(ref)
log('training RegressionModel (250 epochs)...')
reg.train(max_epochs=250, accelerator=ACCELERATOR)
ref = reg.export_posterior(
    ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator': ACCELERATOR})
reg.save(str(C2L_REF_DIR), overwrite=True)

if 'means_per_cluster_mu_fg' in ref.varm:
    inf_aver = ref.varm['means_per_cluster_mu_fg'].copy()
else:
    inf_aver = ref.var[[c for c in ref.var.columns if 'means_per_cluster_mu_fg' in c]].copy()
inf_aver.columns = [c.replace('means_per_cluster_mu_fg_', '') for c in inf_aver.columns]
inf_aver.to_csv(REF_DIR_P / 'inf_aver_signatures.csv')
log('signatures:', inf_aver.shape)

# ======================================================================
# 3 · spatial mapping
# ======================================================================
from cell2location.models import Cell2location  # noqa: E402

assert H5AD.exists(), f'mount the processing volume: {H5AD}'
log('loading visium:', H5AD)
vis = sc.read_h5ad(H5AD)
if COUNT_LAYER in vis.layers:
    vis.X = vis.layers[COUNT_LAYER].copy()
vis.var_names_make_unique()

shared = [g for g in inf_aver.index if g in vis.var_names]
log('shared genes vis n signatures:', len(shared))
vis = vis[:, shared].copy()
sig = inf_aver.loc[shared]

Cell2location.setup_anndata(vis, batch_key=SAMPLE_KEY)
mod = Cell2location(vis, cell_state_df=sig, N_cells_per_location=30, detection_alpha=20)
log('training Cell2location (5000 epochs, CPU)...')
mod.train(max_epochs=5000, batch_size=None, train_size=1, accelerator=ACCELERATOR)
vis = mod.export_posterior(
    vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator': ACCELERATOR})
mod.save(str(C2L_VIS_DIR), overwrite=True)
log('spatial mapping done; obsm:', list(vis.obsm.keys()))

# ======================================================================
# 4 · post-process -> fractions, tables, figures
# ======================================================================
frac = dv.abundance_to_fractions(vis)
dv.add_fractions_to_obs(vis, frac)
frac.to_csv(TBL / 'per_spot_fractions.csv')

mye_cols = [c for c in frac.columns if c.startswith('MYE:')]
vis.obs['myeloid_frac'] = frac[mye_cols].sum(axis=1).values
vis.obs['dominant_type'] = frac.idxmax(axis=1).values
log('\nmean composition:')
log(frac.mean().sort_values(ascending=False).round(4).to_string())
log('\ntotal immune fraction (mean/spot):', round(frac[mye_cols].sum(axis=1).mean(), 4))
vis.write_h5ad(REF_DIR_P / 'ST_BRICHOS_deconvolved.h5ad')

# region x treatment composition
comp = frac.copy()
comp[REGION_KEY] = vis.obs[REGION_KEY].values
comp[TREATMENT_KEY] = vis.obs[TREATMENT_KEY].values
comp = comp.groupby([REGION_KEY, TREATMENT_KEY], observed=True).mean().reset_index()
comp.to_csv(TBL / 'region_treatment_composition.csv', index=False)

# spatial maps: aggregate myeloid + top 3 myeloid states
to_map = ['myeloid_frac'] + list(frac[mye_cols].mean().sort_values(ascending=False).head(3).index)
for feat in to_map:
    col = feat if feat == 'myeloid_frac' else f'prop_{feat}'
    vmax = float(np.nanpercentile(vis.obs[col], 99)) or 1e-6
    for treat in ['WT', 'PBS', 'BRICHOS']:
        sub = vis[vis.obs[TREATMENT_KEY] == treat]
        if sub.n_obs == 0:
            continue
        libs = sub.obs[SAMPLE_KEY].unique()
        fig, axes = plt.subplots(1, len(libs), figsize=(3 * len(libs), 3.2), squeeze=False)
        for ax, lib in zip(axes.flat, libs):
            sl = sub[sub.obs[SAMPLE_KEY] == lib]
            sp = vis.uns['spatial'][lib]
            sf = sp['scalefactors']['tissue_hires_scalef']
            ax.imshow(sp['images']['hires'], origin='upper')
            xy = sl.obsm['spatial'] * sf
            ax.scatter(xy[:, 0], xy[:, 1], c=sl.obs[col].values, s=0.6,
                       cmap='magma', vmin=0, vmax=vmax, edgecolors='none')
            ax.set_title(f'{lib} [{treat}]', fontsize=8, loc='left')
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
        fig.suptitle(f'{feat}: proportion (pruned ref)', fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG / f'{FIG_PREFIX}{feat.replace(":", "_")}_{treat}.png',
                    dpi=200, bbox_inches='tight')
        plt.close(fig)

# region heatmap (PBS)
comp_pbs = comp[comp[TREATMENT_KEY] == 'PBS'].set_index(REGION_KEY)[mye_cols]
fig, ax = plt.subplots(figsize=(max(6, len(mye_cols) * 1.1), 6))
im = ax.imshow(comp_pbs.values, aspect='auto', cmap='viridis')
ax.set_xticks(range(len(mye_cols))); ax.set_xticklabels(mye_cols, rotation=90)
ax.set_yticks(range(len(comp_pbs.index))); ax.set_yticklabels(comp_pbs.index)
ax.set_title('Myeloid-state mean fraction by region (PBS, pruned ref)')
fig.colorbar(im, ax=ax, shrink=0.7)
fig.tight_layout(); fig.savefig(FIG / f'{FIG_PREFIX}myeloid_region_heatmap.png', dpi=200,
                                bbox_inches='tight')
plt.close(fig)

# BRI vs PBS shift, sample-level
from scipy.stats import mannwhitneyu  # noqa: E402
df = frac.copy()
df[[SAMPLE_KEY, REGION_KEY, TREATMENT_KEY]] = vis.obs[[SAMPLE_KEY, REGION_KEY, TREATMENT_KEY]].values
samp = df.groupby([REGION_KEY, TREATMENT_KEY, SAMPLE_KEY], observed=True)[mye_cols].mean().reset_index()
rows = []
for region in samp[REGION_KEY].unique():
    r = samp[samp[REGION_KEY] == region]
    bri = r[r[TREATMENT_KEY] == 'BRICHOS']; pbs_ = r[r[TREATMENT_KEY] == 'PBS']
    if len(bri) < 2 or len(pbs_) < 2:
        continue
    for st in mye_cols:
        try:
            _, p = mannwhitneyu(bri[st], pbs_[st])
        except ValueError:
            p = np.nan
        rows.append(dict(region=region, state=st, mean_pbs=pbs_[st].mean(),
                         mean_bri=bri[st].mean(), delta=bri[st].mean() - pbs_[st].mean(), mwu_p=p))
shift = pd.DataFrame(rows).sort_values('mwu_p')
shift.to_csv(TBL / 'myeloid_state_BRIvPBS_shift.tsv', sep='\t', index=False)
log('\nwrote tables to', TBL)
log('DONE_PRUNED')
