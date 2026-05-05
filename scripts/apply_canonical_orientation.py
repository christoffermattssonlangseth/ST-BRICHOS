"""Read rotations.yaml, apply to combined h5ad, write oriented copy.

Usage
-----
    python scripts/apply_canonical_orientation.py
    python scripts/apply_canonical_orientation.py \
        --in  /Volumes/processing2/ST_BRICHOS/data/ST_BRICHOS_new.h5ad \
        --out /Volumes/processing2/ST_BRICHOS/data/ST_BRICHOS_oriented.h5ad

The default I/O path uses notebooks/data (a symlink to the processing
volume). Downstream notebooks should switch their `read_h5ad` calls to
the oriented file once you are happy with the preview.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import scanpy as sc
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from utils.canonical_orientation import apply_orientation_table  # noqa: E402

BASEDIR = Path("/Volumes/processing2/ST_BRICHOS/data")
DEFAULT_IN = BASEDIR / "ST_BRICHOS_region_subcluster.h5ad"
DEFAULT_OUT = BASEDIR / "ST_BRICHOS_region_subcluster_oriented.h5ad"
ROT_YAML = ROOT / "data" / "rotations.yaml"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", default=str(DEFAULT_IN))
    ap.add_argument("--out", dest="out_path", default=str(DEFAULT_OUT))
    ap.add_argument("--rotations", default=str(ROT_YAML))
    ap.add_argument("--sample-key", default="library_id")
    args = ap.parse_args()

    print(f"reading {args.in_path}")
    adata = sc.read_h5ad(args.in_path)

    with open(args.rotations) as f:
        table = yaml.safe_load(f) or {}
    print(f"loaded {len(table)} rotation entries from {args.rotations}")

    apply_orientation_table(adata, table,
                            sample_key=args.sample_key, verbose=True)

    print(f"writing {args.out_path}")
    adata.write_h5ad(args.out_path, compression="gzip")
    print("done.")


if __name__ == "__main__":
    main()
