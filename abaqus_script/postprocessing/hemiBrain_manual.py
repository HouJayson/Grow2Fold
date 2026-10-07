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

"""Split whole-brain cortical surfaces into left and right hemispheres.

The midsagittal plane is estimated from a reference slice VTK.  Each frame is
clipped by that infinite plane and only the newly created planar cut contours
are capped.  Each closed hemisphere is then conservatively remeshed over its
entire surface and reprojected onto the pre-remeshing geometry.

Requirements
------------
    numpy, pyvista, vtk, pymeshlab

Optional (only when SAVE_GIFTI=True)
-------------------------------------
    nibabel
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import re

import numpy as np
import pyvista as pv
import vtk
import pymeshlab


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
CASE_DIR = Path( r"path-to-vtkfile")
WHOLE_DIR = CASE_DIR / "whole"
TEMP_DIR = CASE_DIR / "temp"
SLICE_FILE = WHOLE_DIR / "slice_file.vtk"

# x is the left-right direction for this dataset.  If the coordinates use a
# different convention, change this to "y" or "z".
LEFT_RIGHT_AXIS = "x"

EXPECTED_FRAME_COUNT = 26
CLEAN_TOLERANCE = 1.0e-6
SAVE_GIFTI = True

# Each worker processes one complete frame.  Start with two because PyMeshLab
# remeshing is memory-intensive.  Set to 1 to disable multiprocessing.
MAX_WORKERS = 4

# Conservative whole-hemisphere isotropic remeshing.  A target scale of 1.0
# matches the target edge length to the original median edge length.
REMESH_ENTIRE_HEMISPHERE = True
REMESH_TARGET_EDGE_SCALE = 1.0
REMESH_ITERATIONS = 10

# With a target equal to the original median edge length, the main need is to
# split the very long cap edges.  Disabling collapse makes ten iterations much
# faster and better preserves the existing cortical sampling density.
REMESH_COLLAPSE_EDGES = True

# Maximum permitted departure from the original clipped surface, expressed as
# a fraction of the target edge length.  Reprojection is also enabled below.
# Increase cautiously to 0.10 if the remesher cannot regularize the cap.
REMESH_MAX_DEVIATION_SCALE = 0.05

# Edges with a dihedral angle above this threshold are treated as features.
# The relatively low value helps preserve the cortex-cap junction.
REMESH_FEATURE_ANGLE_DEG = 15.0

# Computing point-to-original-surface distances is useful but expensive.
# Only the listed zero-based frames receive the remesh_distance VTK scalar.
# Use set() for maximum speed or set(range(26)) to check every frame.
REMESH_DEVIATION_QC_FRAMES = {0, 12, 25}


# ---------------------------------------------------------------------------
# Mesh and plane utilities
# ---------------------------------------------------------------------------
def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_polydata(path: Path, clean_tol: float = CLEAN_TOLERANCE) -> pv.PolyData:
    """Read a surface as cleaned, triangular PolyData."""
    mesh = pv.read(str(path))
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()
    if mesh.n_points == 0 or mesh.n_cells == 0:
        raise ValueError(f"Mesh is empty: {path}")
    return mesh.triangulate().clean(tolerance=clean_tol)


def fit_plane_from_slice(slice_path: Path):
    """Fit a best-fit plane to all points in the reference slice VTK.

    Returns
    -------
    origin : ndarray, shape (3,)
        Centroid of the slice points.
    normal : ndarray, shape (3,)
        Unit normal of the least-squares plane.
    diagnostics : dict
        RMS and maximum point-to-plane distances and slice dimensions.
    """
    slice_mesh = pv.read(str(slice_path))
    points = np.asarray(slice_mesh.points, dtype=np.float64)
    if points.shape[0] < 3:
        raise ValueError(f"At least three slice points are required: {slice_path}")

    origin = points.mean(axis=0)
    centered = points - origin
    _, singular_values, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    normal /= np.linalg.norm(normal)

    # Give the normal a reproducible sign, approximately toward the positive
    # left-right coordinate.  Clipping itself does not depend on this choice.
    axis_idx = {"x": 0, "y": 1, "z": 2}[LEFT_RIGHT_AXIS.lower()]
    if normal[axis_idx] < 0:
        normal *= -1.0

    distances = centered @ normal
    diagnostics = {
        "n_points": int(points.shape[0]),
        "rms_distance": float(np.sqrt(np.mean(distances**2))),
        "max_distance": float(np.max(np.abs(distances))),
        "singular_values": singular_values,
    }
    return origin, normal, diagnostics


def clip_and_cap_one_side(
    mesh: pv.PolyData,
    origin: np.ndarray,
    normal: np.ndarray,
) -> pv.PolyData:
    """Keep one plane half-space and cap only contours made by this clip.

    vtkClipClosedSurface generates faces solely on the supplied clipping plane;
    it does not invoke a general hole-filling operation on the rest of the mesh.
    """
    plane = vtk.vtkPlane()
    plane.SetOrigin(*(float(v) for v in origin))
    plane.SetNormal(*(float(v) for v in normal))

    planes = vtk.vtkPlaneCollection()
    planes.AddItem(plane)

    clipper = vtk.vtkClipClosedSurface()
    clipper.SetInputData(mesh)
    clipper.SetClippingPlanes(planes)
    clipper.SetGenerateFaces(1)
    clipper.SetTriangulationErrorDisplay(1)
    clipper.Update()

    output = pv.wrap(clipper.GetOutput()).triangulate()
    if output.n_points == 0 or output.n_cells == 0:
        raise RuntimeError("Clipping returned an empty surface. Check the slice plane.")
    return output.clean(tolerance=CLEAN_TOLERANCE)


def split_hemispheres(mesh: pv.PolyData, origin, normal):
    """Clip both sides and label them from their centroid coordinates."""
    side_a = clip_and_cap_one_side(mesh, origin, normal)
    side_b = clip_and_cap_one_side(mesh, origin, -normal)

    axis_idx = {"x": 0, "y": 1, "z": 2}[LEFT_RIGHT_AXIS.lower()]
    centroid_a = side_a.center[axis_idx]
    centroid_b = side_b.center[axis_idx]

    # Conventional Cartesian labeling: lower LR coordinate = left.  This is
    # easy to reverse here if the source files use the opposite convention.
    if centroid_a < centroid_b:
        left, right = side_a, side_b
    else:
        left, right = side_b, side_a
    return left, right


def unique_edge_lengths(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Return lengths of all unique triangle edges."""
    edges = np.vstack(
        (
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        )
    )
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    lengths = np.linalg.norm(points[edges[:, 0]] - points[edges[:, 1]], axis=1)
    return lengths[np.isfinite(lengths) & (lengths > 0.0)]


def remesh_entire_hemisphere(
    mesh: pv.PolyData,
    origin: np.ndarray,
    normal: np.ndarray,
    target_edge_scale: float = REMESH_TARGET_EDGE_SCALE,
    iterations: int = REMESH_ITERATIONS,
    max_deviation_scale: float = REMESH_MAX_DEVIATION_SCALE,
    feature_angle_deg: float = REMESH_FEATURE_ANGLE_DEG,
    collapse_edges: bool = REMESH_COLLAPSE_EDGES,
    measure_deviation: bool = False,
) -> pv.PolyData:
    """Conservatively remesh the complete capped hemisphere.

    The target edge length is based on the median original edge length.  Every
    remeshing operation is distance-limited and new vertices are reprojected
    onto the original clipped surface.  The input should not yet contain any
    metrics that need to be retained because topology-changing remeshing cannot
    preserve point-wise arrays without a separate interpolation step.
    """
    surface = mesh.triangulate().clean(tolerance=CLEAN_TOLERANCE)
    points = np.asarray(surface.points, dtype=np.float64)
    faces = np.asarray(surface.faces).reshape(-1, 4)[:, 1:4].astype(np.int32)

    edge_lengths = unique_edge_lengths(points, faces)
    if edge_lengths.size == 0:
        raise RuntimeError("Cannot determine a target edge length for remeshing.")

    median_edge = float(np.median(edge_lengths))
    target_edge = float(target_edge_scale) * median_edge
    max_deviation = float(max_deviation_scale) * target_edge

    bounds_size = points.max(axis=0) - points.min(axis=0)
    bbox_diagonal = float(np.linalg.norm(bounds_size))
    if bbox_diagonal <= 0.0:
        raise RuntimeError("Hemisphere bounding box has zero size.")

    # PercentageValue is measured relative to the bounding-box diagonal.
    target_percent = 100.0 * target_edge / bbox_diagonal
    deviation_percent = 100.0 * max_deviation / bbox_diagonal

    print(
        "     remesh settings: "
        f"median edge={median_edge:.6g}, target={target_edge:.6g}, "
        f"max deviation={max_deviation:.6g}"
    )

    mesh_set = pymeshlab.MeshSet()
    mesh_set.add_mesh(
        pymeshlab.Mesh(vertex_matrix=points, face_matrix=faces),
        "capped_hemisphere",
    )

    mesh_set.apply_filter(
        "meshing_isotropic_explicit_remeshing",
        iterations=int(iterations),
        adaptive=False,
        selectedonly=False,
        targetlen=pymeshlab.PercentageValue(target_percent),
        featuredeg=float(feature_angle_deg),
        checksurfdist=True,
        maxsurfdist=pymeshlab.PercentageValue(deviation_percent),
        splitflag=True,
        collapseflag=bool(collapse_edges),
        swapflag=True,
        smoothflag=True,
        reprojectflag=True,
    )

    # Remove only artifacts that can be introduced by topology operations.  Do
    # not invoke MeshLab's general close-holes filter here.
    for filter_name in (
        "meshing_remove_duplicate_faces",
        "meshing_remove_duplicate_vertices",
        "meshing_remove_null_faces",
        "meshing_remove_unreferenced_vertices",
    ):
        try:
            mesh_set.apply_filter(filter_name)
        except Exception:
            pass

    result = mesh_set.current_mesh()
    new_points = np.asarray(result.vertex_matrix(), dtype=np.float64)
    new_faces = np.asarray(result.face_matrix(), dtype=np.int64)
    if new_points.size == 0 or new_faces.size == 0:
        raise RuntimeError("Whole-hemisphere remeshing produced an empty mesh.")

    vtk_faces = np.column_stack(
        (np.full(new_faces.shape[0], 3, dtype=np.int64), new_faces)
    ).ravel()
    output = pv.PolyData(new_points, vtk_faces).clean(tolerance=CLEAN_TOLERANCE)

    if measure_deviation:
        # Quantify geometric change only for requested QC frames.  The scalar
        # is retained in VTK for inspection in ParaView.
        distance_mesh = output.compute_implicit_distance(surface, inplace=False)
        remesh_distance = np.abs(
            np.asarray(distance_mesh.point_data["implicit_distance"], dtype=float)
        )
        output.point_data["remesh_distance"] = remesh_distance
        print(
            "     measured surface deviation: "
            f"median={np.median(remesh_distance):.6g}, "
            f"95th={np.percentile(remesh_distance, 95):.6g}, "
            f"max={np.max(remesh_distance):.6g}"
        )
        if np.max(remesh_distance) > 1.25 * max_deviation:
            print(
                "     WARNING: measured maximum deviation exceeds the requested "
                "distance limit by more than 25%."
            )

    # Recreate an explicit artificial-cap mask for later curvature/MSM masking.
    # Reprojection should keep the cap planar; the scale-aware tolerance also
    # accommodates small numerical deviations introduced by MeshLab.
    normal = np.asarray(normal, dtype=np.float64)
    normal /= np.linalg.norm(normal)
    origin = np.asarray(origin, dtype=np.float64)
    cap_tolerance = max(10.0 * CLEAN_TOLERANCE, 1.0e-6 * bbox_diagonal)

    output_faces = np.asarray(output.faces).reshape(-1, 4)[:, 1:4]
    plane_distance = np.abs((np.asarray(output.points) - origin) @ normal)
    cap_faces = np.all(plane_distance[output_faces] <= cap_tolerance, axis=1)
    output.cell_data["cap_face"] = cap_faces.astype(np.int32)

    medial_wall = np.zeros(output.n_points, dtype=np.int32)
    if np.any(cap_faces):
        medial_wall[np.unique(output_faces[cap_faces])] = 1
    else:
        print("     WARNING: the planar cap could not be reidentified after remeshing.")
    output.point_data["medial_wall"] = medial_wall

    return output


def count_boundary_edges_near_plane(
    mesh: pv.PolyData,
    origin: np.ndarray,
    normal: np.ndarray,
    plane_tolerance: float,
) -> int:
    """Count uncapped boundary line cells whose vertices lie on the cut plane."""
    edges = mesh.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    if edges.n_cells == 0:
        return 0

    distances = np.abs((np.asarray(edges.points) - origin) @ normal)
    count = 0
    offset = 0
    lines = np.asarray(edges.lines)
    while offset < lines.size:
        n_vertices = int(lines[offset])
        ids = lines[offset + 1 : offset + 1 + n_vertices]
        if ids.size and np.all(distances[ids] <= plane_tolerance):
            count += 1
        offset += n_vertices + 1
    return count


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------
def save_legacy_vtk(mesh: pv.PolyData, output_path: Path):
    ensure_dir(output_path.parent)
    mesh.save(str(output_path), binary=False)


def save_gifti(mesh: pv.PolyData, output_path: Path):
    try:
        import nibabel as nib
    except ImportError as exc:
        raise ImportError("Install nibabel or set SAVE_GIFTI=False.") from exc

    surface = mesh.triangulate()
    points = np.asarray(surface.points, dtype=np.float32)
    faces = np.asarray(surface.faces).reshape(-1, 4)[:, 1:4].astype(np.int32)
    image = nib.gifti.GiftiImage(
        darrays=[
            nib.gifti.GiftiDataArray(
                points,
                intent=nib.nifti1.intent_codes["NIFTI_INTENT_POINTSET"],
                datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_FLOAT32"],
            ),
            nib.gifti.GiftiDataArray(
                faces,
                intent=nib.nifti1.intent_codes["NIFTI_INTENT_TRIANGLE"],
                datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_INT32"],
            ),
        ]
    )
    ensure_dir(output_path.parent)
    nib.save(image, str(output_path))


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------
def natural_sort_key(path: Path):
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r"(\d+)", path.name)]


def find_frame_files(folder: Path, slice_path: Path):
    """Return naturally sorted frame VTKs, explicitly excluding the slice."""
    slice_resolved = slice_path.resolve()
    files = [
        path for path in folder.glob("*.vtk")
        if path.resolve() != slice_resolved and "slice" not in path.stem.lower()
    ]
    return sorted(files, key=natural_sort_key)


def process_frame(
    frame_index: int,
    input_path: Path,
    left_dir: Path,
    right_dir: Path,
    gifti_dir: Path,
    origin: np.ndarray,
    normal: np.ndarray,
    plane_tolerance: float,
    measure_deviation: bool,
):
    whole = load_polydata(input_path)
    left, right = split_hemispheres(whole, origin, normal)

    if REMESH_ENTIRE_HEMISPHERE:
        print("     remeshing left hemisphere")
        left = remesh_entire_hemisphere(
            left,
            origin,
            normal,
            measure_deviation=measure_deviation,
        )
        print("     remeshing right hemisphere")
        right = remesh_entire_hemisphere(
            right,
            origin,
            normal,
            measure_deviation=measure_deviation,
        )

    left_vtk = left_dir / f"{input_path.stem}_L.vtk"
    right_vtk = right_dir / f"{input_path.stem}_R.vtk"
    save_legacy_vtk(left, left_vtk)
    save_legacy_vtk(right, right_vtk)

    if SAVE_GIFTI:
        save_gifti(left, gifti_dir / "L" / f"{input_path.stem}_L.surf.gii")
        save_gifti(right, gifti_dir / "R" / f"{input_path.stem}_R.surf.gii")

    # This checks only for uncapped edges at the newly introduced cut.  Any
    # boundary edges away from the plane are intentionally left unchanged.
    left_open = count_boundary_edges_near_plane(
        left, origin, normal, plane_tolerance
    )
    right_open = count_boundary_edges_near_plane(
        right, origin, normal, plane_tolerance
    )
    print(
        f"[OK] {input_path.name}: "
        f"L={left.n_points} points, R={right.n_points} points, "
        f"open cut edges L/R={left_open}/{right_open}"
    )
    if left_open or right_open:
        print("     WARNING: at least one interhemispheric cut contour is not capped.")
    return frame_index, input_path.name


def main():
    if not WHOLE_DIR.is_dir():
        raise FileNotFoundError(f"Whole-brain folder not found: {WHOLE_DIR}")
    if not SLICE_FILE.is_file():
        raise FileNotFoundError(f"Slice VTK not found: {SLICE_FILE}")

    origin, normal, plane_info = fit_plane_from_slice(SLICE_FILE)
    print(f"Plane origin: {np.array2string(origin, precision=6)}")
    print(f"Plane normal: {np.array2string(normal, precision=6)}")
    print(
        "Slice planarity: "
        f"RMS={plane_info['rms_distance']:.6g}, "
        f"max={plane_info['max_distance']:.6g}"
    )

    frames = find_frame_files(WHOLE_DIR, SLICE_FILE)
    if not frames:
        raise FileNotFoundError(f"No frame VTK files found in: {WHOLE_DIR}")
    if len(frames) != EXPECTED_FRAME_COUNT:
        print(
            f"WARNING: found {len(frames)} frame VTKs; "
            f"expected {EXPECTED_FRAME_COUNT}. Processing all files found."
        )

    output_root = ensure_dir(TEMP_DIR / "01_cut_hemispheres")
    left_dir = ensure_dir(output_root / "L")
    right_dir = ensure_dir(output_root / "R")
    gifti_dir = ensure_dir(output_root / "gifti")

    # Scale the cap check tolerance to the reference slice size.
    slice_mesh = pv.read(str(SLICE_FILE))
    scale = max(float(slice_mesh.length), 1.0)
    plane_tolerance = max(10.0 * CLEAN_TOLERANCE, 1.0e-7 * scale)

    worker_count = max(1, min(int(MAX_WORKERS), len(frames)))
    print(f"Processing {len(frames)} frames with {worker_count} worker(s).")

    if worker_count == 1:
        for index, frame in enumerate(frames):
            print(f"[{index + 1:02d}/{len(frames):02d}] Processing {frame.name}")
            process_frame(
                index,
                frame,
                left_dir,
                right_dir,
                gifti_dir,
                origin,
                normal,
                plane_tolerance,
                index in REMESH_DEVIATION_QC_FRAMES,
            )
    else:
        # Separate processes are used instead of threads because every frame
        # owns independent VTK/MeshLab objects and writes unique output files.
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            future_to_frame = {
                executor.submit(
                    process_frame,
                    index,
                    frame,
                    left_dir,
                    right_dir,
                    gifti_dir,
                    origin,
                    normal,
                    plane_tolerance,
                    index in REMESH_DEVIATION_QC_FRAMES,
                ): frame
                for index, frame in enumerate(frames)
            }

            completed = 0
            for future in as_completed(future_to_frame):
                frame = future_to_frame[future]
                # result() propagates worker errors instead of silently
                # continuing with a missing hemisphere.
                _, frame_name = future.result()
                completed += 1
                print(
                    f"[completed {completed:02d}/{len(frames):02d}] "
                    f"{frame_name}"
                )

    print(f"Finished. Outputs: {output_root}")


if __name__ == "__main__":
    # Required for safe process spawning on Windows.
    multiprocessing.freeze_support()
    main()
