#!/usr/bin/env python
"""Orthogonal deconvolution with stereoscope (scvi-tools) — microglia-only ref.

Cross-method check of cell2location run 3 (`run_deconv_micro.py`): same
combined reference pruned to the 6 brain classes + `MYE:Microglia`, same gene
set (the c2l-filtered signature genes, so the two runs see identical input),
but a different model (Andersson et al. 2020; scvi.external implementation,
not the archived almaan/stereoscope CLI).

Stereoscope has no per-spot cell-density prior and no batch term; instead it
fits a per-gene noise "cell type" that soaks up ambient signal. We save its
weight separately (`noise_fraction.csv`) and report proportions renormalized
over the real cell types so they're directly comparable to c2l fractions.

Writes to *_stereoscope paths; c2l results are untouched. Compare in
notebooks/03/09d-stereoscope-vs-cell2location.ipynb.

Run headless on the peer (conda env `c2l`):
    python utils/run_deconv_stereoscope.py
Override epochs/device via env vars, e.g.
    STEREO_SC_EPOCHS=400 STEREO_ST_EPOCHS=2000 STEREO_BATCH=1024 STEREO_ACCEL=cpu python utils/run_deconv_stereoscope.py
"""
from __future__ import annotations
import inspect, os, sys, warnings
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

ACCELERATOR = os.environ.get('STEREO_ACCEL', 'cpu')
SC_EPOCHS = int(os.environ.get('STEREO_SC_EPOCHS', 400))
ST_EPOCHS = int(os.environ.get('STEREO_ST_EPOCHS', 2000))
BATCH = int(os.environ.get('STEREO_BATCH', 1024))   # scvi default 128 is overhead-bound on CPU
N_HVG_FALLBACK = 5000

# --- paths -------------------------------------------------------------
BASEDIR = Path('/Volumes/processing2/ST_BRICHOS/data')
H5AD_ORIENTED = BASEDIR / 'ST_BRICHOS_region_subcluster_oriented.h5ad'
H5AD_BASE     = BASEDIR / 'ST_BRICHOS_region_subcluster.h5ad'
H5AD = H5AD_ORIENTED if H5AD_ORIENTED.exists() else H5AD_BASE
COUNT_LAYER = 'counts'

COMBINED_H5AD = ROOT / 'data' / 'external' / 'deconv_reference' / 'combined_reference.h5ad'
C2L_MICRO_SIG = ROOT / 'data' / 'external' / 'deconv_reference_micro' / 'inf_aver_signatures.csv'
REF_DIR_S = ROOT / 'data' / 'external' / 'deconv_reference_stereoscope'
REF_DIR_S.mkdir(parents=True, exist_ok=True)

SAMPLE_KEY, REGION_KEY, TREATMENT_KEY = 'sample_id', 're_annotation_regions', 'treatment'

TBL = ROOT / 'results' / 'tables' / 'deconvolution_stereoscope'; TBL.mkdir(parents=True, exist_ok=True)
FIG = ROOT / 'results' / 'figures' / 'manuscript'; FIG.mkdir(parents=True, exist_ok=True)
FIG_PREFIX = 'deconv_stereo_'

KEEP_MYELOID = {'MYE:Microglia'}   # identical to c2l run 3


def train_kwargs(fn) -> dict:
    """scvi-tools >=1.0 takes `accelerator`; older builds take `use_gpu`."""
    params = inspect.signature(fn).parameters
    if 'accelerator' in params:
        return {'accelerator': ACCELERATOR}
    return {'use_gpu': ACCELERATOR not in ('cpu', None)}


# ======================================================================
# 1 · reference: prune to brain classes + microglia
# ======================================================================
log('loading combined reference:', COMBINED_H5AD)
ref = sc.read_h5ad(COMBINED_H5AD)
ct = ref.obs['cell_type'].astype(str)
keep_mask = (~ct.str.startswith('MYE:')) | ct.isin(KEEP_MYELOID)
ref = ref[keep_mask.values].copy()
ref.obs['cell_type'] = ref.obs['cell_type'].astype(str).astype('category')
log('  pruned reference:', ref.shape, '|', ref.obs.cell_type.nunique(), 'types')
log(ref.obs['cell_type'].value_counts().to_string())

# ======================================================================
# 2 · visium + shared gene set
# ======================================================================
assert H5AD.exists(), f'mount the processing volume: {H5AD}'
log('loading visium:', H5AD)
vis = sc.read_h5ad(H5AD)
if COUNT_LAYER in vis.layers:
    vis.X = vis.layers[COUNT_LAYER].copy()
vis.var_names_make_unique()
vis.layers['counts'] = vis.X.copy()

if C2L_MICRO_SIG.exists():
    genes = pd.read_csv(C2L_MICRO_SIG, index_col=0, usecols=[0]).index
    log('gene set: c2l run-3 signature genes', len(genes))
else:
    log(f'no {C2L_MICRO_SIG.name}; falling back to {N_HVG_FALLBACK} seurat_v3 HVGs on the reference')
    sc.pp.highly_variable_genes(ref, n_top_genes=N_HVG_FALLBACK, flavor='seurat_v3',
                                layer='counts', batch_key='ref_source')
    genes = ref.var_names[ref.var['highly_variable']]
genes = [g for g in genes if g in ref.var_names and g in vis.var_names]
log('shared genes used:', len(genes))
ref = ref[:, genes].copy()
vis = vis[:, genes].copy()

# ======================================================================
# 3 · stereoscope
# ======================================================================
from scvi.external import RNAStereoscope, SpatialStereoscope  # noqa: E402

RNAStereoscope.setup_anndata(ref, layer='counts', labels_key='cell_type')
sc_model = RNAStereoscope(ref)
log(f'training RNAStereoscope ({SC_EPOCHS} epochs, batch {BATCH}, {ACCELERATOR})...')
sc_model.train(max_epochs=SC_EPOCHS, batch_size=BATCH, **train_kwargs(sc_model.train))
sc_model.save(str(REF_DIR_S / 'rna_model'), overwrite=True)

SpatialStereoscope.setup_anndata(vis, layer='counts')
st_model = SpatialStereoscope.from_rna_model(vis, sc_model)
log(f'training SpatialStereoscope ({ST_EPOCHS} epochs, batch {BATCH}, {ACCELERATOR})...')
st_model.train(max_epochs=ST_EPOCHS, batch_size=BATCH, **train_kwargs(st_model.train))
st_model.save(str(REF_DIR_S / 'spatial_model'), overwrite=True)

for name, m in (('rna', sc_model), ('spatial', st_model)):
    hist = m.history.get('elbo_train', m.history.get('train_loss_epoch'))
    if hist is not None:
        hist.to_csv(TBL / f'loss_{name}.csv')
        log(f'  {name} loss first/last:', float(hist.iloc[0, 0]), '->', float(hist.iloc[-1, 0]))

# ======================================================================
# 4 · proportions -> tables, figures
# ======================================================================
raw = st_model.get_proportions(keep_noise=True)
noise_cols = [c for c in raw.columns if c not in ref.obs['cell_type'].cat.categories]
noise = raw[noise_cols].sum(axis=1)
frac = raw.drop(columns=noise_cols)
frac = frac.div(frac.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
frac.index = vis.obs_names
noise.index = vis.obs_names

frac.to_csv(TBL / 'per_spot_fractions.csv')
noise.rename('noise_fraction').to_csv(TBL / 'noise_fraction.csv')
log('\nmean composition:')
log(frac.mean().sort_values(ascending=False).round(4).to_string())
log('noise term (mean/spot):', round(float(noise.mean()), 4))

comp = frac.copy()
comp[REGION_KEY] = vis.obs[REGION_KEY].values
comp[TREATMENT_KEY] = vis.obs[TREATMENT_KEY].values
comp = comp.groupby([REGION_KEY, TREATMENT_KEY], observed=True).mean().reset_index()
comp.to_csv(TBL / 'region_treatment_composition.csv', index=False)

# microglia spatial map per treatment
col = 'MYE:Microglia'
vmax = float(np.nanpercentile(frac[col], 99)) or 1e-6
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
        ax.scatter(xy[:, 0], xy[:, 1], c=frac.loc[sl.obs_names, col].values, s=0.6,
                   cmap='magma', vmin=0, vmax=vmax, edgecolors='none')
        ax.set_title(f'{lib} [{treat}]', fontsize=8, loc='left')
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle(f'{col}: proportion (stereoscope)', fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / f'{FIG_PREFIX}MYE_Microglia_{treat}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

log('\nwrote tables to', TBL)
log('DONE_STEREOSCOPE')
