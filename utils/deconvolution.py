"""Helpers for cell2location deconvolution of the BRICHOS Visium data.

Scope
-----
The reference is a *myeloid-only* atlas (Chhatbar et al. 2026, Nat Immunol;
Zenodo 16938034). A myeloid-only reference cannot deconvolve a whole 55 um
Visium spot on its own, so we build a *combined* reference:

    Allen cortex+hippocampus (major non-myeloid types: neurons, astro,
    oligo, OPC, endo, VLMC, ...)          <-- structural background
  + Chhatbar mouse myeloid atlas states    <-- fine microglia/myeloid states
  = combined reference with one unified `cell_type` column

We strip Allen's own microglia/myeloid before grafting so the myeloid
compartment is described only by the fine atlas states.

This module is deliberately dependency-light (scanpy/anndata/scipy) so the
reference wrangling runs without cell2location/torch installed. The
cell2location model calls themselves live in the notebook.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad


# --------------------------------------------------------------------------
# 1 · load the Matrix-Market bundle written by utils/export_seurat_atlas.R
# --------------------------------------------------------------------------
def load_mtx_bundle(bundle_dir: str | Path) -> ad.AnnData:
    """Read matrix.mtx.gz + features/barcodes/metadata into AnnData.

    The R exporter writes counts as genes x cells (10x convention); we
    transpose to cells x genes. Raw counts land in both ``.X`` and
    ``.layers['counts']`` (cell2location's reference model wants raw counts).
    """
    d = Path(bundle_dir)
    with gzip.open(d / "matrix.mtx.gz", "rb") as fh:
        m = sio.mmread(fh).tocsr()            # genes x cells
    genes = _read_lines(d / "features.tsv.gz")
    bcs = _read_lines(d / "barcodes.tsv.gz")
    meta = pd.read_csv(d / "metadata.csv.gz")

    X = m.T.tocsr()                            # cells x genes
    if X.shape[0] != len(bcs) or X.shape[1] != len(genes):
        raise ValueError(
            f"shape mismatch: matrix {X.shape} vs {len(bcs)} cells / "
            f"{len(genes)} genes")

    obs = meta.copy()
    if "barcode" in obs.columns:
        obs.index = obs["barcode"].astype(str).values
    else:
        obs.index = np.asarray(bcs, dtype=str)
    var = pd.DataFrame(index=pd.Index(genes, name=None))

    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    return a


def _read_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [ln.rstrip("\n").strip('"') for ln in fh]


# --------------------------------------------------------------------------
# 2 · build the combined reference (Allen background + myeloid states)
# --------------------------------------------------------------------------
# Substrings (lower-cased) used to detect Allen's own myeloid/immune cells so
# we can strip them before grafting the atlas myeloid states. Tune to the
# label vocabulary of whichever Allen build you load.
DEFAULT_MYELOID_PATTERNS = (
    "microglia", "micro-", "pvm", "perivascular macrophage", "macrophage",
    "myeloid", "monocyte", "dendritic", "border-associated", "bam",
    "immune",
)


def strip_myeloid(
    brain_ref: ad.AnnData,
    label_col: str,
    patterns: tuple[str, ...] = DEFAULT_MYELOID_PATTERNS,
) -> tuple[ad.AnnData, list[str]]:
    """Drop cells whose ``label_col`` matches any myeloid pattern.

    Returns (filtered_ref, dropped_labels) so the notebook can print exactly
    which Allen labels were removed (sanity-check before grafting).
    """
    lab = brain_ref.obs[label_col].astype(str)
    low = lab.str.lower()
    hit = np.zeros(len(lab), dtype=bool)
    for p in patterns:
        hit |= low.str.contains(p, regex=False)
    dropped = sorted(lab[hit].unique().tolist())
    return brain_ref[~hit].copy(), dropped


def build_combined_reference(
    brain_ref: ad.AnnData,
    myeloid_ref: ad.AnnData,
    *,
    brain_label_col: str,
    myeloid_label_col: str,
    myeloid_prefix: str = "MYE:",
    counts_layer: str = "counts",
) -> ad.AnnData:
    """Intersect genes, concatenate, and write a unified ``cell_type`` col.

    Myeloid state labels are prefixed (default ``MYE:``) so they never
    collide with Allen major-type names and are easy to pull back out of the
    deconvolution result. Both inputs must carry raw counts in
    ``.layers[counts_layer]``.
    """
    for name, a in (("brain_ref", brain_ref), ("myeloid_ref", myeloid_ref)):
        if counts_layer not in a.layers:
            raise KeyError(f"{name} missing raw counts in .layers['{counts_layer}']")

    shared = brain_ref.var_names.intersection(myeloid_ref.var_names)
    if len(shared) < 2000:
        raise ValueError(
            f"only {len(shared)} shared genes — check gene symbol casing "
            "(both must be mouse MGI symbols)")
    b = brain_ref[:, shared].copy()
    m = myeloid_ref[:, shared].copy()

    b.obs["cell_type"] = b.obs[brain_label_col].astype(str).values
    m.obs["cell_type"] = myeloid_prefix + m.obs[myeloid_label_col].astype(str).values
    b.obs["ref_source"] = "brain"
    m.obs["ref_source"] = "myeloid"

    # keep only what the model needs; counts in .X
    b.X = b.layers[counts_layer].copy()
    m.X = m.layers[counts_layer].copy()
    keep = ["cell_type", "ref_source"]
    b.obs = b.obs[keep]
    m.obs = m.obs[keep]

    combined = ad.concat([b, m], join="outer", label=None, index_unique="-")
    combined.layers["counts"] = combined.X.copy()
    combined.obs["cell_type"] = combined.obs["cell_type"].astype("category")
    combined.obs["ref_source"] = combined.obs["ref_source"].astype("category")
    return combined


def subsample_per_type(
    a: ad.AnnData, label_col: str, n_per: int = 500, seed: int = 0
) -> ad.AnnData:
    """Cap cells per label — keeps the reference regression fast on CPU."""
    rng = np.random.default_rng(seed)
    idx: list[int] = []
    codes = a.obs[label_col].astype(str).values
    for lab in pd.unique(codes):
        pos = np.where(codes == lab)[0]
        if len(pos) > n_per:
            pos = rng.choice(pos, n_per, replace=False)
        idx.extend(pos.tolist())
    idx.sort()
    return a[idx].copy()


# --------------------------------------------------------------------------
# 3 · post-process cell2location output
# --------------------------------------------------------------------------
def abundance_to_fractions(
    adata_vis: ad.AnnData, q05_key: str = "q05_cell_abundance_w_sf"
) -> pd.DataFrame:
    """Turn cell2location q05 abundances into per-spot fractions (rows sum≈1).

    cell2location stores the posterior 5% quantile of absolute abundance in
    ``adata_vis.obsm[q05_key]`` with columns named
    ``q05cell_abundance_w_sf_<cell_type>``. We normalize each spot to a
    composition. Absolute abundances remain available in ``.obsm``.
    """
    if q05_key not in adata_vis.obsm:
        raise KeyError(
            f"{q05_key} not in .obsm — run cell2location export_posterior first")
    ab = adata_vis.obsm[q05_key].copy()
    ab.columns = [c.split("w_sf_")[-1] for c in ab.columns]
    frac = ab.div(ab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return frac


def add_fractions_to_obs(
    adata_vis: ad.AnnData, frac: pd.DataFrame, prefix: str = "prop_"
) -> None:
    """Write per-cell-type fractions into ``.obs`` as ``prop_<cell_type>``."""
    for ct in frac.columns:
        adata_vis.obs[f"{prefix}{ct}"] = frac[ct].reindex(adata_vis.obs_names).values


def region_composition(
    adata_vis: ad.AnnData,
    frac: pd.DataFrame,
    region_key: str,
    treatment_key: str,
) -> pd.DataFrame:
    """Mean per-spot fraction by region x treatment (long form)."""
    df = frac.copy()
    df[region_key] = adata_vis.obs[region_key].values
    df[treatment_key] = adata_vis.obs[treatment_key].values
    long = df.groupby([region_key, treatment_key], observed=True).mean()
    return long
