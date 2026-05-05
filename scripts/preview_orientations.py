"""Render every Visium section in a grid for orientation eyeballing.

Reads the combined h5ad (path configurable via CLI), draws each library's
hires image with spots overlaid, labels with sample id + treatment +
current angle/flip from data/rotations.yaml. Output:

    results/figures/manuscript/orientation_preview.png
    results/figures/manuscript/orientation_preview.svg

Usage
-----
    python scripts/preview_orientations.py
    python scripts/preview_orientations.py --h5ad /path/to/combined.h5ad
    python scripts/preview_orientations.py --apply-rotations  # preview AFTER applying yaml
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from utils.canonical_orientation import apply_orientation_table  # noqa: E402

BASEDIR = Path("/Volumes/processing2/ST_BRICHOS/data")
DEFAULT_H5AD = BASEDIR / "ST_BRICHOS_region_subcluster.h5ad"
ROT_YAML = ROOT / "data" / "rotations.yaml"
OUT_DIR = ROOT / "results" / "figures" / "manuscript"

TREAT_COLOR = {"WT": "#4a4a4a", "PBS": "#1f77b4", "BRICHOS": "#d62728"}


def load_rotations() -> dict:
    if not ROT_YAML.exists():
        return {}
    with open(ROT_YAML) as f:
        return yaml.safe_load(f) or {}


def draw_sample(ax, adata, lib_id: str, sample_key: str,
                treatment_key: str, rot_spec: dict):
    sp = adata.uns["spatial"][lib_id]
    img = sp["images"]["hires"]
    sf = sp["scalefactors"]["tissue_hires_scalef"]
    ax.imshow(img, origin="upper")

    mask = (adata.obs[sample_key] == lib_id).to_numpy()
    xy = adata.obsm["spatial"][mask] * sf
    treat = (adata.obs.loc[mask, treatment_key].iloc[0]
             if treatment_key in adata.obs.columns and mask.any()
             else "?")
    color = TREAT_COLOR.get(treat, "#888")
    ax.scatter(xy[:, 0], xy[:, 1], s=0.4, c=color, alpha=0.55,
               edgecolors="none")

    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    angle = rot_spec.get("angle_deg", 0)
    flip = rot_spec.get("flip_h", False)
    suffix = f" / rot {angle:+g}" + ("/flip" if flip else "")
    ax.set_title(f"{lib_id}  [{treat}]{suffix}", fontsize=8,
                 color=color, loc="left", pad=2)

    # crosshair through image centre to help judging axis alignment
    H, W = img.shape[:2]
    ax.axhline(H / 2, color="#888", linewidth=0.4, alpha=0.5)
    ax.axvline(W / 2, color="#888", linewidth=0.4, alpha=0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--sample-key", default="library_id",
                    help="obs column matching adata.uns['spatial'] keys "
                         "(default: library_id; try sample_id if missing)")
    ap.add_argument("--treatment-key", default="treatment")
    ap.add_argument("--apply-rotations", action="store_true",
                    help="apply rotations.yaml before previewing")
    ap.add_argument("--out", default=str(OUT_DIR / "orientation_preview"))
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"reading {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    if args.sample_key not in adata.obs.columns:
        candidates = [c for c in adata.obs.columns
                      if "sample" in c.lower() or "library" in c.lower()]
        raise SystemExit(
            f"obs column {args.sample_key!r} not found. "
            f"candidates: {candidates}"
        )

    rotations = load_rotations()
    if args.apply_rotations:
        print("applying rotations.yaml")
        apply_orientation_table(adata, rotations,
                                sample_key=args.sample_key, verbose=True)

    libs = list(adata.uns.get("spatial", {}).keys())
    libs.sort()
    n = len(libs)
    if n == 0:
        raise SystemExit("no spatial libraries found in adata.uns['spatial']")

    ncols = min(4, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.0 * ncols, 3.2 * nrows),
                             squeeze=False)
    for ax, lib in zip(axes.flat, libs):
        draw_sample(ax, adata, lib, args.sample_key, args.treatment_key,
                    rotations.get(lib, {}))
    for ax in axes.flat[n:]:
        ax.axis("off")

    fig.suptitle(
        ("Orientation preview — AFTER rotations.yaml"
         if args.apply_rotations else
         "Orientation preview — current (no rotation applied)"),
        fontsize=10, y=0.995)
    fig.tight_layout()
    out_png = Path(args.out + ".png")
    out_svg = Path(args.out + ".svg")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_svg}")


if __name__ == "__main__":
    main()
