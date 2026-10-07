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

"""Write legacy ASCII VTK from a surface GIFTI and scalar data.

Continuous arrays are supplied as metric GIFTIs with repeated ``--scalar``
arguments.  Discrete arrays such as parcellations can be copied directly from
a reference VTK with ``--reference-vtk`` and repeated ``--copy-scalar``
arguments.  When reference arrays are copied, both vertex count and triangle
topology are verified before writing.
"""

import argparse
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pyvista as pv


def load_surf_gii(path: Path):
    image = nib.load(str(path))
    points = None
    triangles = None

    for data_array in image.darrays:
        array = np.asarray(data_array.data)
        if data_array.intent == nib.nifti1.intent_codes["NIFTI_INTENT_POINTSET"]:
            points = array.astype(np.float64)
        elif data_array.intent == nib.nifti1.intent_codes["NIFTI_INTENT_TRIANGLE"]:
            triangles = array.astype(np.int64)

    if points is None or triangles is None:
        raise ValueError(f"Surface GIFTI lacks pointset or triangles: {path}")
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"Surface points must have shape N x 3: {path}")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError(f"Surface triangles must have shape M x 3: {path}")
    if triangles.size and (
        triangles.min() < 0 or triangles.max() >= points.shape[0]
    ):
        raise ValueError(f"Surface contains invalid triangle indices: {path}")
    return points, triangles


def load_metric_gii(path: Path, expected_length: int) -> np.ndarray:
    image = nib.load(str(path))
    if not image.darrays:
        raise ValueError(f"Metric GIFTI contains no data arrays: {path}")

    values = np.asarray(image.darrays[0].data).squeeze()
    if values.ndim != 1:
        raise ValueError(f"Metric must be one-dimensional: {path}; got {values.shape}")
    if values.shape[0] != expected_length:
        raise ValueError(
            f"Metric length mismatch for {path}: "
            f"{values.shape[0]} versus {expected_length}"
        )
    return values.astype(np.float64)


def parse_scalar_item(item: str):
    if "=" not in item:
        raise ValueError(
            f"Invalid --scalar '{item}'; expected name=/path/to/metric.func.gii"
        )
    name, path = (part.strip() for part in item.split("=", 1))
    if not name or not path:
        raise ValueError(f"Invalid --scalar '{item}'")
    return name, Path(path)


def load_reference_arrays(
    reference_path: Path,
    names: list[str],
    expected_points: int,
    expected_triangles: np.ndarray,
):
    reference = pv.read(str(reference_path))
    if not isinstance(reference, pv.PolyData):
        reference = reference.extract_surface()
    reference = reference.triangulate()

    if reference.n_points != expected_points:
        raise ValueError(
            f"Reference VTK has {reference.n_points} vertices, but the registered "
            f"surface has {expected_points}"
        )

    reference_triangles = (
        np.asarray(reference.faces, dtype=np.int64).reshape(-1, 4)[:, 1:4]
    )
    if (
        reference_triangles.shape != expected_triangles.shape
        or not np.array_equal(reference_triangles, expected_triangles)
    ):
        raise ValueError(
            "Reference and registered surfaces do not have identical triangle "
            "topology; copying labels by vertex index would be unsafe"
        )

    arrays = {}
    available = list(reference.point_data.keys())
    for name in names:
        if name not in reference.point_data:
            raise KeyError(
                f"Reference VTK lacks point scalar '{name}'. Available: {available}"
            )
        values = np.asarray(reference.point_data[name]).squeeze()
        if values.ndim != 1 or values.shape[0] != expected_points:
            raise ValueError(f"Reference scalar '{name}' is not a length-N array")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Reference label scalar '{name}' contains NaN or Inf")

        # Parcellations are categorical even when the source VTK stored them as
        # floating-point values.  Fail rather than silently changing noninteger
        # labels.
        rounded = np.rint(values)
        if not np.allclose(values, rounded, rtol=0.0, atol=1.0e-6):
            raise ValueError(f"Reference label scalar '{name}' is not integer-valued")
        arrays[name] = rounded.astype(np.int32)
    return arrays


def write_legacy_vtk_polydata(
    output_path: Path,
    points: np.ndarray,
    triangles: np.ndarray,
    continuous_scalars: dict[str, np.ndarray],
    integer_scalars: dict[str, np.ndarray],
    title: str,
):
    number_of_points = points.shape[0]
    number_of_triangles = triangles.shape[0]

    with output_path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("# vtk DataFile Version 2.0\n")
        stream.write(f"{title.replace(chr(10), ' ')}\n")
        stream.write("ASCII\n")
        stream.write("DATASET POLYDATA\n")

        stream.write(f"POINTS {number_of_points} float\n")
        for point in points:
            stream.write(f"{point[0]:.9g} {point[1]:.9g} {point[2]:.9g}\n")

        stream.write(
            f"POLYGONS {number_of_triangles} {number_of_triangles * 4}\n"
        )
        for triangle in triangles:
            stream.write(
                f"3 {int(triangle[0])} {int(triangle[1])} {int(triangle[2])}\n"
            )

        if continuous_scalars or integer_scalars:
            stream.write(f"POINT_DATA {number_of_points}\n")

        for name, values in continuous_scalars.items():
            stream.write(f"SCALARS {name} float 1\n")
            stream.write("LOOKUP_TABLE default\n")
            for value in values:
                if np.isfinite(value):
                    stream.write(f"{float(value):.9g}\n")
                else:
                    stream.write("nan\n")

        for name, values in integer_scalars.items():
            stream.write(f"SCALARS {name} int 1\n")
            stream.write("LOOKUP_TABLE default\n")
            for value in values:
                stream.write(f"{int(value)}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Write legacy ASCII VTK from GIFTI geometry and scalar arrays."
    )
    parser.add_argument("--surf", required=True, type=Path)
    parser.add_argument(
        "--scalar",
        action="append",
        default=[],
        help="Continuous scalar: name=/path/to/metric.func.gii",
    )
    parser.add_argument(
        "--reference-vtk",
        type=Path,
        help="Reference VTK providing categorical point scalars",
    )
    parser.add_argument(
        "--copy-scalar",
        action="append",
        default=[],
        help="Point scalar to copy from --reference-vtk; repeat as needed",
    )
    parser.add_argument("--title", default="surface_with_scalars")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    if not args.surf.is_file():
        raise FileNotFoundError(f"Surface file not found: {args.surf}")
    if args.copy_scalar and args.reference_vtk is None:
        raise ValueError("--reference-vtk is required with --copy-scalar")
    if args.reference_vtk is not None and not args.reference_vtk.is_file():
        raise FileNotFoundError(f"Reference VTK not found: {args.reference_vtk}")

    points, triangles = load_surf_gii(args.surf)
    continuous = {}
    for item in args.scalar:
        name, metric_path = parse_scalar_item(item)
        if name in continuous:
            raise ValueError(f"Duplicate scalar name: {name}")
        if not metric_path.is_file():
            raise FileNotFoundError(f"Metric file not found: {metric_path}")
        continuous[name] = load_metric_gii(metric_path, points.shape[0])

    integer = {}
    if args.copy_scalar:
        duplicates = set(continuous).intersection(args.copy_scalar)
        if duplicates:
            raise ValueError(f"Scalar names supplied twice: {sorted(duplicates)}")
        integer = load_reference_arrays(
            args.reference_vtk,
            args.copy_scalar,
            points.shape[0],
            triangles,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_legacy_vtk_polydata(
        args.out,
        points,
        triangles,
        continuous,
        integer,
        args.title,
    )
    print(f"[OK] Wrote {args.out}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
