"""Statistical helpers for the BRICHOS attenuation analysis.

Functions here build on pseudobulk counts (one row per mouse x region)
and provide:

* `make_pseudobulk`     — collapse spot-level counts to mouse x region.
* `pseudobulk_lfc`      — log2 fold-change with limma-voom-ish weights.
* `sign_concordance`    — fraction of disease-DEGs reversing under treatment.
* `attenuation_slope`   — Pearson slope/intercept of LFC_residual on
                          LFC_disease, with mouse-level bootstrap CI.
* `responder_zshift`    — per-spot z-shift toward WT for a signature score.

Conventions
-----------
LFC_disease  = mean(PBS) - mean(WT)        (what disease does)
LFC_residual = mean(BRI) - mean(WT)        (what is left after BRICHOS)
LFC_rescue   = mean(BRI) - mean(PBS)       (the BRI-vs-PBS contrast)

For attenuation, slope b of LFC_residual on LFC_disease (over genes
significant in PBS vs WT) summarises rescue:
    b ~ 1   no rescue
    b < 1   attenuation
    b ~ 0   full rescue
    b < 0   overshoot
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats


# ---------- pseudobulk -------------------------------------------------------

def make_pseudobulk(adata, sample_key: str, region_key: str,
                     treatment_key: str = "treatment",
                     min_spots: int = 20,
                     layer: str | None = None):
    """Collapse spot counts to (sample, region) pseudobulk.

    Parameters
    ----------
    adata : AnnData with raw counts in .X (or .layers[layer]).
    sample_key, region_key, treatment_key : obs columns.
    min_spots : drop (sample, region) cells with fewer spots than this.
    layer : count layer to sum. None = .X (must be raw integer counts).

    Returns
    -------
    counts : DataFrame (n_groups, n_genes) — summed counts.
    meta   : DataFrame (n_groups, [sample, region, treatment, n_spots]).
    """
    X = adata.X if layer is None else adata.layers[layer]
    if sp.issparse(X):
        X_csr = X.tocsr()
    else:
        X_csr = sp.csr_matrix(X)

    obs = adata.obs[[sample_key, region_key, treatment_key]].copy()
    obs["__row__"] = np.arange(adata.n_obs)
    grouped = obs.groupby([sample_key, region_key, treatment_key],
                          observed=True)

    rows, meta_records = [], []
    for (samp, reg, treat), idx in grouped.indices.items():
        if len(idx) < min_spots:
            continue
        rows.append(np.asarray(X_csr[idx].sum(axis=0)).ravel())
        meta_records.append(dict(sample=samp, region=reg, treatment=treat,
                                 n_spots=len(idx)))

    counts = pd.DataFrame(np.vstack(rows),
                          columns=adata.var_names,
                          index=[f"{m['sample']}__{m['region']}" for m in meta_records])
    meta = pd.DataFrame(meta_records, index=counts.index)
    return counts, meta


# ---------- per-gene log-fold change ----------------------------------------

@dataclass
class LFCResult:
    lfc: pd.Series
    se: pd.Series
    pvalue: pd.Series
    padj: pd.Series
    n_a: int
    n_b: int

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({"lfc": self.lfc, "se": self.se,
                             "pvalue": self.pvalue, "padj": self.padj})


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def pseudobulk_lfc(counts: pd.DataFrame, meta: pd.DataFrame,
                    group_a: str, group_b: str,
                    region: str | None = None,
                    pseudocount: float = 1.0,
                    min_count: int = 10,
                    min_n: int = 2) -> LFCResult:
    """log2 fold-change A vs B on log-CPM pseudobulk.

    Welch t-test on log2(CPM + pseudocount) when both groups have
    `>= 2` replicates. With `min_n=1` and one group having a single
    replicate, falls back to pooled-variance (assume equal variance,
    borrow from the multi-replicate group). Inflates significance —
    use only for descriptive ranking when WT n=1, never to claim FDR
    on individual genes.

    For publication-grade per-gene stats, swap in pyDESeq2.
    """
    m = meta
    if region is not None:
        m = m[m["region"] == region]
    a = m.index[m["treatment"] == group_a]
    b = m.index[m["treatment"] == group_b]
    n_a, n_b = len(a), len(b)
    if n_a < min_n or n_b < min_n:
        raise ValueError(f"need >={min_n} replicates in each group; got "
                         f"{group_a}: {n_a}, {group_b}: {n_b}")
    sub = counts.loc[list(a) + list(b)]
    keep = (sub.sum(axis=0) >= min_count).values
    sub = sub.loc[:, keep]

    lib = sub.sum(axis=1).values[:, None]
    cpm = sub.values / np.where(lib == 0, 1, lib) * 1e6
    log2cpm = np.log2(cpm + pseudocount)

    A = log2cpm[: n_a]
    B = log2cpm[n_a:]
    mean_a = A.mean(axis=0)
    mean_b = B.mean(axis=0)

    if n_a >= 2 and n_b >= 2:
        var_a = A.var(axis=0, ddof=1)
        var_b = B.var(axis=0, ddof=1)
        se = np.sqrt(var_a / n_a + var_b / n_b)
        se = np.where(se == 0, np.nan, se)
        t = (mean_a - mean_b) / se
        df = (var_a / n_a + var_b / n_b) ** 2 / (
            (var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1)
        )
    else:
        # single replicate on one side: pool variance from the multi side.
        # df = n_multi - 1; SE assumes equal variance.
        if n_a == 1:
            var_pool = B.var(axis=0, ddof=1)
            df = n_b - 1
        else:
            var_pool = A.var(axis=0, ddof=1)
            df = n_a - 1
        se = np.sqrt(var_pool * (1.0 / n_a + 1.0 / n_b))
        se = np.where(se == 0, np.nan, se)
        t = (mean_a - mean_b) / se

    p = 2 * stats.t.sf(np.abs(t), df=df)
    padj = _bh(np.nan_to_num(p, nan=1.0))

    genes = sub.columns
    return LFCResult(
        lfc=pd.Series(mean_a - mean_b, index=genes, name="lfc"),
        se=pd.Series(se, index=genes, name="se"),
        pvalue=pd.Series(p, index=genes, name="pvalue"),
        padj=pd.Series(padj, index=genes, name="padj"),
        n_a=n_a, n_b=n_b,
    )


# ---------- rescue helpers (BRI vs PBS evaluated against a disease list) ----

def signed_rescue(lfc_rescue: pd.Series,
                   disease_up: Iterable[str],
                   disease_dn: Iterable[str] | None = None) -> dict:
    """Test whether disease-up genes go DOWN under treatment (and vice-versa).

    Inputs
    ------
    lfc_rescue : per-gene LFC of treatment vs disease-control
                 (positive => treatment increases the gene).
    disease_up : genes elevated in disease (PBS vs WT, or literature PIG).
    disease_dn : optional, genes lowered in disease.

    Returns dict with sign-concordance (binomial vs 0.5) and a one-sample
    Wilcoxon test on signed LFC pooled across the signature.
    """
    lfc_rescue = lfc_rescue.dropna()

    up = pd.Index(disease_up).intersection(lfc_rescue.index)
    dn = (pd.Index(disease_dn).intersection(lfc_rescue.index)
          if disease_dn is not None else pd.Index([]))

    # signed expectation: rescue moves up-genes down (-) and dn-genes up (+).
    pooled = pd.concat([
        -lfc_rescue.loc[up],
        lfc_rescue.loc[dn],
    ])
    pooled = pooled.replace([np.inf, -np.inf], np.nan).dropna()
    if len(pooled) == 0:
        return dict(n=0, frac_concordant=np.nan,
                    binom_p=np.nan, wilcoxon_p=np.nan,
                    median_signed_lfc=np.nan)

    concordant = (pooled > 0).sum()
    n = int(len(pooled))
    binom_p = stats.binomtest(int(concordant), n, p=0.5,
                              alternative="greater").pvalue
    try:
        w_p = stats.wilcoxon(pooled, alternative="greater").pvalue
    except ValueError:
        w_p = np.nan
    return dict(n=n, n_up=int(len(up)), n_dn=int(len(dn)),
                frac_concordant=float(concordant) / n,
                binom_p=float(binom_p),
                wilcoxon_p=float(w_p),
                median_signed_lfc=float(pooled.median()))


def rescue_slope(lfc_disease: pd.Series, lfc_rescue: pd.Series,
                  sig_mask=None,
                  n_boot: int = 1000, seed: int = 0) -> dict:
    """Slope of LFC_rescue (BRI-vs-PBS) on LFC_disease (PBS-vs-WT).

    Negative slope ~ rescue (treatment opposes disease). Slope magnitude
    near -1 = full rescue, 0 = no effect, >0 = exacerbation.

    `sig_mask` may be:
      * pandas Series of bool indexed like `lfc_disease`,
      * boolean ndarray aligned with `lfc_disease.index`,
      * iterable of gene names to keep,
      * None (use all common genes).
    """
    common = lfc_disease.index.intersection(lfc_rescue.index)
    d = lfc_disease.loc[common].values
    r = lfc_rescue.loc[common].values
    if sig_mask is not None:
        if isinstance(sig_mask, pd.Series):
            sig = sig_mask.reindex(common).fillna(False).astype(bool).values
        elif isinstance(sig_mask, np.ndarray) and sig_mask.dtype == bool:
            sig = (pd.Series(sig_mask, index=lfc_disease.index)
                   .reindex(common).fillna(False).astype(bool).values)
        else:
            sig = common.isin(list(sig_mask))
        d, r = d[sig], r[sig]
    if len(d) < 3:
        return dict(slope=np.nan, intercept=np.nan,
                    ci_low=np.nan, ci_high=np.nan, n=len(d))
    rng = np.random.default_rng(seed)
    slopes = np.empty(n_boot)
    n = len(d)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        slopes[i] = np.polyfit(d[idx], r[idx], 1)[0]
    fit = np.polyfit(d, r, 1)
    return dict(slope=float(fit[0]), intercept=float(fit[1]),
                ci_low=float(np.nanpercentile(slopes, 2.5)),
                ci_high=float(np.nanpercentile(slopes, 97.5)),
                n=int(n))


# ---------- sign concordance ------------------------------------------------

def sign_concordance(lfc_disease: pd.Series, lfc_rescue: pd.Series,
                     sig_mask: pd.Series) -> dict:
    """Fraction of disease-DEGs whose treatment LFC reverses sign.

    A gene "reverses" when sign(lfc_disease) != sign(lfc_rescue) AND
    both are non-zero. Tests that fraction against 0.5 (binomial).
    """
    common = lfc_disease.index.intersection(lfc_rescue.index)\
                              .intersection(sig_mask.index)
    d = lfc_disease.loc[common]
    r = lfc_rescue.loc[common]
    sig = sig_mask.loc[common].astype(bool)
    if sig.sum() == 0:
        return dict(n_sig=0, n_reverse=0, fraction=np.nan, p=np.nan)
    d, r = d[sig], r[sig]
    nz = (d != 0) & (r != 0)
    d, r = d[nz], r[nz]
    reverses = (np.sign(d) != np.sign(r))
    k = int(reverses.sum())
    n = int(len(reverses))
    p = stats.binomtest(k, n, p=0.5, alternative="greater").pvalue
    return dict(n_sig=n, n_reverse=k, fraction=k / n if n else np.nan,
                p=p)


# ---------- attenuation slope ----------------------------------------------

def _slope(x: np.ndarray, y: np.ndarray) -> float:
    """Slope of y on x through origin? No — standard OLS with intercept."""
    if len(x) < 3:
        return np.nan
    s = np.polyfit(x, y, 1)
    return float(s[0])


def attenuation_slope(lfc_disease: pd.Series, lfc_residual: pd.Series,
                       sig_mask: pd.Series, n_boot: int = 1000,
                       seed: int = 0) -> dict:
    """Slope of LFC_residual on LFC_disease over disease-DEGs, with CI.

    The CI here is gene-level bootstrap (resample genes with replacement).
    For mouse-level CI, bootstrap mice in the upstream pseudobulk and
    pass the resulting slopes through `bootstrap_mice` instead.
    """
    common = lfc_disease.index.intersection(lfc_residual.index)\
                              .intersection(sig_mask.index)
    d = lfc_disease.loc[common].values
    r = lfc_residual.loc[common].values
    sig = sig_mask.loc[common].astype(bool).values
    d, r = d[sig], r[sig]
    if len(d) < 3:
        return dict(slope=np.nan, intercept=np.nan, ci_low=np.nan,
                    ci_high=np.nan, n=len(d))

    rng = np.random.default_rng(seed)
    slopes = np.empty(n_boot)
    n = len(d)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        slopes[i] = _slope(d[idx], r[idx])
    slope_obs = _slope(d, r)
    intercept = float(np.polyfit(d, r, 1)[1])
    return dict(slope=slope_obs, intercept=intercept,
                ci_low=float(np.nanpercentile(slopes, 2.5)),
                ci_high=float(np.nanpercentile(slopes, 97.5)),
                n=int(n))


# ---------- responder spot scoring ------------------------------------------

def responder_zshift(adata, score_key: str, treatment_key: str,
                      region_key: str, sample_key: str,
                      ref: str = "WT") -> pd.Series:
    """Per-spot z-shift of a signature score relative to regional WT mean.

    For each region, compute mean and SD of `score_key` across `ref`
    spots. Return (score - ref_mean) / ref_sd for every spot. Negative
    z in BRI/PBS spots = closer to WT than that region's WT itself
    (rare); positive z = pathological.
    """
    obs = adata.obs
    out = pd.Series(np.nan, index=obs.index, name=f"{score_key}_zshift")
    for region in obs[region_key].dropna().unique():
        m_region = obs[region_key] == region
        ref_mask = m_region & (obs[treatment_key] == ref)
        if ref_mask.sum() < 5:
            continue
        ref_mean = obs.loc[ref_mask, score_key].mean()
        ref_sd = obs.loc[ref_mask, score_key].std(ddof=1)
        if not np.isfinite(ref_sd) or ref_sd == 0:
            continue
        out.loc[m_region] = (obs.loc[m_region, score_key] - ref_mean) / ref_sd
    return out


# ---------- mouse-level bootstrap helper -----------------------------------

def bootstrap_mice(counts: pd.DataFrame, meta: pd.DataFrame,
                    region: str, statistic, n_boot: int = 200,
                    seed: int = 0) -> np.ndarray:
    """Resample mice with replacement within each treatment, recompute
    `statistic(counts_boot, meta_boot, region)`, return array of values.

    `statistic` must accept (counts, meta, region) and return a scalar
    (e.g., the attenuation slope for the region).
    """
    rng = np.random.default_rng(seed)
    sub = meta[meta["region"] == region]
    by_treat = {t: sub[sub["treatment"] == t]["sample"].unique()
                for t in sub["treatment"].unique()}
    out = np.empty(n_boot)
    for i in range(n_boot):
        picked = []
        for t, mice in by_treat.items():
            picked.extend(rng.choice(mice, size=len(mice), replace=True))
        # rebuild meta/counts using picked mice; duplicates allowed
        rows = [meta.index[(meta["sample"] == m) & (meta["region"] == region)]
                for m in picked]
        rows = [r for sub_rows in rows for r in sub_rows]
        if not rows:
            out[i] = np.nan
            continue
        out[i] = statistic(counts.loc[rows], meta.loc[rows], region)
    return out
