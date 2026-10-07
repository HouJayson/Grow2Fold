#!/usr/bin/env python3
# =============================================================================
# Copyright (c) 2026 Jixin Hou et al.
# All rights reserved.
#
# This code is provided as part of the research software accompanying:
# [Grow2Fold: mapping heterogeneous developmental growth to human brain folding]
#
# Use, modification, and redistribution are permitted under the terms of the
# license provided in the LICENSE file of this repository.
#
# Repository: [https://github.com/BioDMX-UGA/Grow2Fold]
# =============================================================================
import argparse
from pathlib import Path
import numpy as np
import pyvista as pv
import nibabel as nib


def write_surf_gii(points: np.ndarray, tris: np.ndarray, out_path: Path) -> None:
    coordsys = nib.gifti.GiftiCoordSystem(
        dataspace=nib.nifti1.xform_codes["NIFTI_XFORM_UNKNOWN"],
        xformspace=nib.nifti1.xform_codes["NIFTI_XFORM_UNKNOWN"],
        xform=np.eye(4, dtype=np.float64),
    )

    da_pts = nib.gifti.GiftiDataArray(
        points.astype(np.float32),
        intent=nib.nifti1.intent_codes["NIFTI_INTENT_POINTSET"],
        datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_FLOAT32"],
        coordsys=coordsys,
    )
    da_tri = nib.gifti.GiftiDataArray(
        tris.astype(np.int32),
        intent=nib.nifti1.intent_codes["NIFTI_INTENT_TRIANGLE"],
        datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_INT32"],
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.gifti.GiftiImage(darrays=[da_pts, da_tri]), str(out_path))


def write_metric_gii(values: np.ndarray, out_path: Path) -> None:
    v = np.asarray(values).reshape(-1).astype(np.float32)
    da = nib.gifti.GiftiDataArray(
        v,
        intent=nib.nifti1.intent_codes["NIFTI_INTENT_SHAPE"],
        datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_FLOAT32"],
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.gifti.GiftiImage(darrays=[da]), str(out_path))


def vtk_to_gii(vtk_path: Path, out_dir: Path, name: str, metrics: list[str], allow_missing=False) -> None:
    mesh = pv.read(str(vtk_path))
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()
    # Do not clean or merge points: vertex order must remain identical to the
    # reference VTK so target-space parcellations can later be copied safely.
    mesh = mesh.triangulate()
    pts = np.asarray(mesh.points, dtype=np.float32)

    faces = np.asarray(mesh.faces, dtype=np.int64).reshape(-1, 4)
    tris = faces[:, 1:4].astype(np.int32)

    surf_out = out_dir / f"{name}.surf.gii"
    write_surf_gii(pts, tris, surf_out)

    n = pts.shape[0]
    available = set(mesh.point_data.keys())

    for m in metrics:
        if m not in available:
            if allow_missing:
                arr = np.zeros(n, dtype=np.float32)
            else:
                raise RuntimeError(f"Missing point-data array '{m}' in {vtk_path}")
        else:
            arr = np.asarray(mesh.point_data[m]).squeeze()
            if arr.ndim != 1:
                raise RuntimeError(
                    f"{m}: expected a scalar point array, got shape {arr.shape}"
                )
            if arr.shape[0] != n:
                raise RuntimeError(f"{m}: length mismatch ({arr.shape[0]} vs {n})")

        write_metric_gii(arr, out_dir / f"{name}.{m}.func.gii")

    print(f"[OK] {vtk_path.name} -> {surf_out.name} (+{len(metrics)} metrics)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vtk", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--name", required=True, type=str)
    ap.add_argument("--metrics", nargs="*", default=[])
    ap.add_argument("--allow-missing", action="store_true")
    args = ap.parse_args()

    vtk_to_gii(args.vtk, args.out, args.name, args.metrics, allow_missing=args.allow_missing)


if __name__ == "__main__":
    main()
