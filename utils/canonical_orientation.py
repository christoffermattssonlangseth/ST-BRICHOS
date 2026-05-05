"""Canonical orientation for Visium sections.

Rotates each sample's spot coordinates and tissue images so that all
sections face the same direction in figure panels. Use together with
`data/rotations.yaml` (per-sample angle + flip).

Conventions
-----------
* `angle_deg` is counter-clockwise in the displayed image (positive
  angle rotates the section CCW). Sign matches `scipy.ndimage.rotate`.
* `flip_h` mirrors left/right *after* rotation. Use for sections that
  were captured on the wrong side of the slide.
* `obsm['spatial']` is treated as full-resolution pixel coordinates,
  per scanpy/squidpy convention. Both `hires` and `lowres` images are
  rotated.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np
from scipy.ndimage import rotate as nd_rotate


def _rotation_matrix(angle_deg: float) -> np.ndarray:
    """2D rotation matrix matching scipy.ndimage.rotate convention.

    scipy rotates the array CCW by +angle (in image (row, col) space).
    To move the *points* the same way, we apply the same CCW rotation
    in (x, y) space. Image y-axis points down, so the matrix is:

        [[cos a, -sin a],
         [sin a,  cos a]]

    applied to (x, y) coordinates relative to the rotation centre.
    """
    a = np.deg2rad(angle_deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s], [s, c]])


def rotate_visium_sample(adata, library_id: str, angle_deg: float,
                          flip_h: bool = False, sample_key: str = "library_id",
                          verbose: bool = False) -> None:
    """Rotate one Visium sample in-place.

    Rotates `adata.obsm['spatial']` rows belonging to `library_id` and
    the `hires`/`lowres` images stored under
    `adata.uns['spatial'][library_id]['images']`. Centre of rotation is
    the centre of the hires image (in fullres pixel units), so spots
    and image stay aligned.
    """
    if library_id not in adata.uns.get("spatial", {}):
        raise KeyError(f"library_id {library_id!r} not in adata.uns['spatial']")

    sp = adata.uns["spatial"][library_id]
    sf_h = sp["scalefactors"]["tissue_hires_scalef"]

    img_h = sp["images"].get("hires")
    img_l = sp["images"].get("lowres")
    if img_h is None:
        raise KeyError(f"hires image missing for {library_id}")

    H_h, W_h = img_h.shape[:2]
    cx_full = (W_h / sf_h) / 2.0
    cy_full = (H_h / sf_h) / 2.0

    if angle_deg != 0:
        img_h = nd_rotate(img_h, angle_deg, reshape=False, order=1,
                          mode="constant",
                          cval=1.0 if img_h.dtype.kind == "f" else 255)
        if img_l is not None:
            img_l = nd_rotate(img_l, angle_deg, reshape=False, order=1,
                              mode="constant",
                              cval=1.0 if img_l.dtype.kind == "f" else 255)

    if flip_h:
        img_h = img_h[:, ::-1].copy()
        if img_l is not None:
            img_l = img_l[:, ::-1].copy()

    sp["images"]["hires"] = img_h
    if img_l is not None:
        sp["images"]["lowres"] = img_l

    if sample_key not in adata.obs.columns:
        raise KeyError(f"obs column {sample_key!r} missing — cannot select spots")
    mask = (adata.obs[sample_key] == library_id).to_numpy()
    if not mask.any():
        if verbose:
            print(f"  no spots for {library_id}, image rotated only")
        return

    xy = adata.obsm["spatial"][mask].astype(float).copy()
    xy[:, 0] -= cx_full
    xy[:, 1] -= cy_full
    if angle_deg != 0:
        xy = xy @ _rotation_matrix(angle_deg).T
    if flip_h:
        xy[:, 0] *= -1
    xy[:, 0] += cx_full
    xy[:, 1] += cy_full

    new_xy = adata.obsm["spatial"].astype(float).copy()
    new_xy[mask] = xy
    adata.obsm["spatial"] = new_xy

    if verbose:
        print(f"  {library_id}: rotated {angle_deg:+.1f} deg"
              f"{' + flip' if flip_h else ''}, {mask.sum()} spots")


def apply_orientation_table(adata, table: dict,
                             sample_key: str = "library_id",
                             verbose: bool = True) -> None:
    """Apply a {library_id: {angle_deg, flip_h}} table to adata in-place.

    Missing samples are left as-is.
    """
    libs: Iterable[str] = list(adata.uns.get("spatial", {}).keys())
    for lib in libs:
        spec = table.get(lib, {})
        angle = float(spec.get("angle_deg", 0.0))
        flip = bool(spec.get("flip_h", False))
        if angle == 0 and not flip:
            if verbose:
                print(f"  {lib}: skip (no rotation/flip)")
            continue
        rotate_visium_sample(adata, lib, angle_deg=angle,
                             flip_h=flip, sample_key=sample_key,
                             verbose=verbose)
