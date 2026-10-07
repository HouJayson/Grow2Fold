#!/usr/bin/env bash
set -euo pipefail

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


# Sequential MSM registration for four independent surface series:
#   L/gray  -> lh.gray,   L/white -> lh.white
#   R/gray  -> rh.gray,   R/white -> rh.white
# For every series:
#   frame 0: curvature -> matching real-reference curvature
#   frame n: sulc(frame n) -> sulc(frame n-1), n >= 1
# Incremental spheres are composed into that series' real-reference sphere space.

log() { echo "[$(date +'%F %T')] $*" >&2; }
die() { log "[ERROR] $*"; exit 1; }

# =============================================================================
# USER SETTINGS
# =============================================================================
CASE_ID="CC00976XX20_unif"

ROOT="/mnt/e/PhD/PINN/WholeBrain_Pipeline/postprocessing/dataSummary/${CASE_ID}"

REF_DIR="${ROOT}/ref"
CUT_ROOT="${ROOT}/temp/01_cut_hemispheres"
SIM_GII_ROOT="${CUT_ROOT}/gifti"
REG_TEMP="${ROOT}/temp/02_sequential_registration"
OUT_ROOT="${ROOT}/hemi"

# Keep the four scripts together in this directory.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
HEMI_BRAIN="${SCRIPT_DIR}/hemiBrain_fast_local_cut.py"
VTK_TO_GII="${SCRIPT_DIR}/vtk_to_gii_metrics.py"
WRITE_VTK="${SCRIPT_DIR}/write_metrics_vtk.py"
CONFIG_DIR="${SCRIPT_DIR}/configs"

MSM_CONFIG="${CONFIG_DIR}/msm_config_init"
NEWMSM_CONFIG="${CONFIG_DIR}/newmsm_config_init"

: "${MSM_BIN:=/usr/local/bin/msm}"
: "${NEWMSM_BIN:=/home/jh30459/msm-env/bin/newmsm}"
: "${WB_BIN:=/usr/local/workbench-1.5.0/bin_linux64/wb_command}"
: "${PYTHON_BIN:=/home/jh30459/miniconda3/envs/surfReg/bin/python}"

USE_NEWMSM=1
EXPECTED_FRAME_COUNT=26
HEMIS=(L R)
SURFACE_TYPES=(gray white)

# Set to 1 to run cutting/capping/remeshing before registration, or 0 to reuse
# existing temp/01_cut_hemispheres outputs.
PERFORM_CUT=0
CUT_OVERWRITE=0
CUT_WORKERS=4

# Independent preparation of spheres, curv, and sulc. MSM remains sequential
# inside each anatomical series, but the four L/R x gray/white series can run
# at the same time. Total peak preparation jobs are approximately
# PREP_JOBS_PER_SERIES * 4 when PARALLEL_SERIES=1.
PREP_JOBS_PER_SERIES=8
PARALLEL_SERIES=2

FEATURE_SMOOTH_SIGMA=2.0

# Retained from the attached script.  Set to 0 if the complete cap/medial wall
# should contribute to the MSM cost without suppression.
MASK_MEDIAL_WALL=0

# When 0, an existing nonempty final VTK is treated as complete. A fully
# completed L/R x gray/white series is skipped before any preparation or MSM.
# Set to 1 only when every final VTK must be regenerated.
OVERWRITE_FINAL_VTK=0

mkdir -p "${REG_TEMP}" "${OUT_ROOT}"

# =============================================================================
# VALIDATION AND STRUCTURE METADATA
# =============================================================================
hemi_fs() { [[ "$1" == "L" ]] && echo "lh" || echo "rh"; }
hemi_struct() { [[ "$1" == "L" ]] && echo "CORTEX_LEFT" || echo "CORTEX_RIGHT"; }

require_executable() { [[ -x "$1" ]] || die "Executable not found: $1"; }

ensure_environment() {
  require_executable "${WB_BIN}"
  require_executable "${PYTHON_BIN}"
  [[ -f "${VTK_TO_GII}" ]] || die "Missing helper: ${VTK_TO_GII}"
  [[ -f "${WRITE_VTK}" ]] || die "Missing helper: ${WRITE_VTK}"
  if [[ "${PERFORM_CUT}" -eq 1 ]]; then
    [[ -f "${HEMI_BRAIN}" ]] || die "Missing cutter: ${HEMI_BRAIN}"
  fi

  command -v mris_convert >/dev/null 2>&1 || die "mris_convert not found"
  command -v mris_smooth  >/dev/null 2>&1 || die "mris_smooth not found"
  command -v mris_sphere  >/dev/null 2>&1 || die "mris_sphere not found"
  command -v mris_inflate >/dev/null 2>&1 || die "mris_inflate not found"

  if [[ "${USE_NEWMSM}" -eq 1 ]]; then
    require_executable "${NEWMSM_BIN}"
    [[ -f "${NEWMSM_CONFIG}" ]] || die "Missing config: ${NEWMSM_CONFIG}"
  else
    require_executable "${MSM_BIN}"
    [[ -f "${MSM_CONFIG}" ]] || die "Missing config: ${MSM_CONFIG}"
  fi

  "${PYTHON_BIN}" - <<PY >/dev/null
import nibabel, numpy, pyvista
if ${PERFORM_CUT}:
    import pymeshlab, vtk
PY
}

set_struct_sphere() {
  "${WB_BIN}" -set-structure "$2" "$(hemi_struct "$1")" \
    -surface-type SPHERICAL >/dev/null 2>&1 || true
}

set_struct_surface() {
  "${WB_BIN}" -set-structure "$2" "$(hemi_struct "$1")" \
    >/dev/null 2>&1 || true
}

set_struct_metric() {
  "${WB_BIN}" -set-structure "$2" "$(hemi_struct "$1")" \
    >/dev/null 2>&1 || true
}

# =============================================================================
# CUTTING
# =============================================================================
run_cutting_if_requested() {
  if [[ "${PERFORM_CUT}" -ne 1 ]]; then
    log "PERFORM_CUT=0: reuse existing hemisphere surfaces"
    return 0
  fi

  log "Run cutting/capping/remeshing for ${CASE_ID}"
  local -a arguments=(
    --case-id "${CASE_ID}"
    --case-root "${ROOT}"
    --workers "${CUT_WORKERS}"
    --expected-frames "$((EXPECTED_FRAME_COUNT * ${#SURFACE_TYPES[@]}))"
  )
  [[ "${CUT_OVERWRITE}" -eq 1 ]] && arguments+=(--overwrite)
  "${PYTHON_BIN}" "${HEMI_BRAIN}" "${arguments[@]}"
}

# =============================================================================
# SPHERE AND FEATURE GENERATION
# =============================================================================
make_sphere() {
  local surface="$1" H="$2" work_dir="$3" tag="$4"
  local HFS; HFS="$(hemi_fs "${H}")"
  local output="${work_dir}/${tag}.sphere.surf.gii"
  mkdir -p "${work_dir}"

  if [[ ! -f "${output}" ]]; then
    log "Generate sphere: ${tag}"
    mris_convert "${surface}" "${work_dir}/${HFS}.white" >/dev/null 2>&1
    mris_smooth -n 1 "${work_dir}/${HFS}.white" \
      "${work_dir}/${HFS}.smoothwm" >/dev/null 2>&1
    mris_sphere "${work_dir}/${HFS}.smoothwm" \
      "${work_dir}/${HFS}.sphere" >/dev/null 2>&1
    mris_convert "${work_dir}/${HFS}.sphere" "${output}" >/dev/null 2>&1
  fi
  set_struct_sphere "${H}" "${output}"
}

curv_metric() {
  local H="$1" surface="$2" output="$3"
  if [[ ! -f "${output}" ]]; then
    mkdir -p "$(dirname "${output}")"
    "${WB_BIN}" -surface-curvature "${surface}" -mean "${output}" \
      >/dev/null 2>&1
  fi
  set_struct_metric "${H}" "${output}"
}

sulc_metric() {
  local H="$1" surface="$2" work_dir="$3" output="$4"
  local HFS; HFS="$(hemi_fs "${H}")"
  if [[ -f "${output}" ]]; then
    set_struct_metric "${H}" "${output}"
    return 0
  fi

  mkdir -p "${work_dir}" "$(dirname "${output}")"
  mris_convert "${surface}" "${work_dir}/${HFS}.white" >/dev/null 2>&1
  mris_inflate "${work_dir}/${HFS}.white" \
    "${work_dir}/${HFS}.inflated" >/dev/null 2>&1

  local fs_sulc="${work_dir}/${HFS}.sulc"
  [[ -f "${fs_sulc}" ]] || die "FreeSurfer did not create ${fs_sulc}"

  "${PYTHON_BIN}" - "${fs_sulc}" "${output}" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np
from nibabel.freesurfer.io import read_morph_data

source, output = sys.argv[1:3]
values = read_morph_data(source).astype(np.float32)
array = nib.gifti.GiftiDataArray(
    values,
    intent=nib.nifti1.intent_codes["NIFTI_INTENT_SHAPE"],
    datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_FLOAT32"],
)
Path(output).parent.mkdir(parents=True, exist_ok=True)
nib.save(nib.gifti.GiftiImage(darrays=[array]), output)
PY
  set_struct_metric "${H}" "${output}"
}

smooth_metric() {
  local H="$1" sphere="$2" input="$3" output="$4"
  if [[ ! -f "${output}" ]]; then
    "${WB_BIN}" -metric-smoothing "${sphere}" "${input}" \
      "${FEATURE_SMOOTH_SIGMA}" "${output}" >/dev/null 2>&1
  fi
  set_struct_metric "${H}" "${output}"
}

mask_feature() {
  local H="$1" input="$2" source_vtk="$3" mode="$4" output="$5"
  if [[ "${MASK_MEDIAL_WALL}" -ne 1 ]]; then
    [[ -f "${output}" ]] || cp -f "${input}" "${output}"
    set_struct_metric "${H}" "${output}"
    return 0
  fi
  if [[ -f "${output}" ]]; then
    set_struct_metric "${H}" "${output}"
    return 0
  fi

  "${PYTHON_BIN}" - "${input}" "${source_vtk}" "${mode}" "${output}" <<'PY'
import sys
from pathlib import Path
import nibabel as nib
import numpy as np
import pyvista as pv

metric_path, vtk_path, mode, output_path = sys.argv[1:5]
values = np.asarray(nib.load(metric_path).darrays[0].data, dtype=np.float32).copy()
mesh = pv.read(vtk_path)
if values.size != mesh.n_points:
    raise ValueError("Feature and VTK vertex counts differ")

if mode == "simulation":
    if "medial_wall" not in mesh.point_data:
        raise KeyError(f"medial_wall not found in {vtk_path}")
    invalid = np.asarray(mesh.point_data["medial_wall"]).squeeze() != 0
elif mode == "reference":
    if "par_huang" not in mesh.point_data:
        raise KeyError(f"par_huang not found in {vtk_path}")
    invalid = np.asarray(mesh.point_data["par_huang"]).squeeze() == 0
else:
    raise ValueError(mode)

values[invalid] = 0.0
array = nib.gifti.GiftiDataArray(
    values,
    intent=nib.nifti1.intent_codes["NIFTI_INTENT_SHAPE"],
    datatype=nib.nifti1.data_type_codes["NIFTI_TYPE_FLOAT32"],
)
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
nib.save(nib.gifti.GiftiImage(darrays=[array]), output_path)
PY
  set_struct_metric "${H}" "${output}"
}

# =============================================================================
# MSM AND RESAMPLING
# =============================================================================
run_msm() {
  local in_sphere="$1" ref_sphere="$2" in_feature="$3" ref_feature="$4"
  local prefix="$5" feature_name="$6"
  local binary config
  if [[ "${USE_NEWMSM}" -eq 1 ]]; then
    binary="${NEWMSM_BIN}"; config="${NEWMSM_CONFIG}"
  else
    binary="${MSM_BIN}"; config="${MSM_CONFIG}"
  fi

  mkdir -p "$(dirname "${prefix}")"
  local output="${prefix}sphere.reg.surf.gii"
  [[ -f "${output}" ]] && return 0

  log "Run MSM feature=${feature_name}: $(basename "${prefix}")"
  if ! "${binary}" --inmesh="${in_sphere}" --refmesh="${ref_sphere}" \
      --indata="${in_feature}" --refdata="${ref_feature}" \
      --out="${prefix}" --conf="${config}" >"${prefix}msm.log" 2>&1; then
    tail -n 50 "${prefix}msm.log" >&2 || true
    return 1
  fi
  [[ -f "${output}" ]] || die "Missing MSM sphere: ${output}"
}

compose_spheres() {
  local H="$1" incremental="$2" previous_native="$3"
  local previous_cumulative="$4" output="$5"
  if [[ ! -f "${output}" ]]; then
    "${WB_BIN}" -surface-sphere-project-unproject "${incremental}" \
      "${previous_native}" "${previous_cumulative}" "${output}" \
      >/dev/null 2>&1
  fi
  set_struct_sphere "${H}" "${output}"
}

resample_surface() {
  [[ -f "$4" ]] || "${WB_BIN}" -surface-resample "$1" "$2" "$3" \
    BARYCENTRIC "$4" >/dev/null 2>&1
}

resample_metric() {
  local H="$1"
  [[ -f "$5" ]] || "${WB_BIN}" -metric-resample "$2" "$3" "$4" \
    BARYCENTRIC "$5" >/dev/null 2>&1
  set_struct_metric "${H}" "$5"
}

# =============================================================================
# PARALLEL FRAME PREPARATION
# =============================================================================
prepare_frame() {
  local H="$1" frame_surface="$2" sim_vtk_dir="$3" work="$4" index="$5"
  local filename base tag frame_work frame_vtk native_sphere
  filename="$(basename "${frame_surface}")"
  base="${filename%.surf.gii}"
  tag="$(printf 'frame%02d' "${index}")"
  frame_work="${work}/frames/${tag}"
  frame_vtk="${sim_vtk_dir}/${base}.vtk"
  native_sphere="${frame_work}/sphere/${base}.sphere.surf.gii"

  [[ -f "${frame_vtk}" ]] || die "Missing simulation VTK: ${frame_vtk}"
  mkdir -p "${frame_work}/features"
  set_struct_surface "${H}" "${frame_surface}"
  make_sphere "${frame_surface}" "${H}" "${frame_work}/sphere" "${base}"

  local curv="${frame_work}/features/native.curv.func.gii"
  local sulc="${frame_work}/features/native.sulc.func.gii"
  local curv_sm="${frame_work}/features/native.curv.sm.func.gii"
  local sulc_sm="${frame_work}/features/native.sulc.sm.func.gii"
  local mask_tag="full"; [[ "${MASK_MEDIAL_WALL}" -eq 1 ]] && mask_tag="masked"
  local curv_feature="${frame_work}/features/native.curv.sm.${mask_tag}.func.gii"
  local sulc_feature="${frame_work}/features/native.sulc.sm.${mask_tag}.func.gii"

  curv_metric "${H}" "${frame_surface}" "${curv}"
  sulc_metric "${H}" "${frame_surface}" "${frame_work}/sulc_work" "${sulc}"
  smooth_metric "${H}" "${native_sphere}" "${curv}" "${curv_sm}"
  smooth_metric "${H}" "${native_sphere}" "${sulc}" "${sulc_sm}"
  mask_feature "${H}" "${curv_sm}" "${frame_vtk}" simulation "${curv_feature}"
  mask_feature "${H}" "${sulc_sm}" "${frame_vtk}" simulation "${sulc_feature}"
}

prepare_all_frames() {
  local H="$1" sim_vtk_dir="$2" work="$3"; shift 3
  local -a surfaces=("$@")
  local running=0 failed=0 index
  for index in "${!surfaces[@]}"; do
    prepare_frame "${H}" "${surfaces[index]}" "${sim_vtk_dir}" \
      "${work}" "${index}" &
    running=$((running + 1))
    if [[ "${running}" -ge "${PREP_JOBS_PER_SERIES}" ]]; then
      if ! wait -n; then failed=1; fi
      running=$((running - 1))
    fi
  done
  while [[ "${running}" -gt 0 ]]; do
    if ! wait -n; then failed=1; fi
    running=$((running - 1))
  done
  [[ "${failed}" -eq 0 ]] || return 1
}

# =============================================================================
# ONE ANATOMICAL SERIES (hemisphere x gray/white)
# =============================================================================
process_series() {
  local H="$1" surface_type="$2" HFS
  HFS="$(hemi_fs "${H}")"
  [[ "${surface_type}" == "gray" || "${surface_type}" == "white" ]] || \
    die "Unsupported surface type: ${surface_type}"

  local ref_vtk="${REF_DIR}/${CASE_ID}.${HFS}.${surface_type}.vtk"
  local sim_gii_dir="${SIM_GII_ROOT}/${H}"
  local sim_vtk_dir="${CUT_ROOT}/${H}"
  local work="${REG_TEMP}/${H}/${surface_type}"
  local output_dir="${OUT_ROOT}/${H}"
  mkdir -p "${work}" "${output_dir}"

  [[ -f "${ref_vtk}" ]] || die "Missing reference VTK: ${ref_vtk}"
  [[ -d "${sim_gii_dir}" ]] || die "Missing GIFTI directory: ${sim_gii_dir}"
  [[ -d "${sim_vtk_dir}" ]] || die "Missing VTK directory: ${sim_vtk_dir}"

  # Discover frames before doing any expensive reference preparation. This
  # allows a completed anatomical series to return immediately.
  local -a frames=() all_hemi_frames=()
  mapfile -t all_hemi_frames < <(find "${sim_gii_dir}" -maxdepth 1 -type f \
    -name "*_${H}.surf.gii" -print | sort -V)
  local candidate candidate_name
  for candidate in "${all_hemi_frames[@]}"; do
    candidate_name="$(basename "${candidate}")"
    candidate_name="${candidate_name,,}"
    if [[ "${candidate_name}" == *"${surface_type}"* ]]; then
      frames+=("${candidate}")
    fi
  done
  [[ "${#frames[@]}" -gt 0 ]] || \
    die "No ${surface_type} frames found for ${H} in ${sim_gii_dir}"
  [[ "${#frames[@]}" -eq "${EXPECTED_FRAME_COUNT}" ]] || \
    die "Found ${#frames[@]} ${H}/${surface_type} frames; expected ${EXPECTED_FRAME_COUNT}"

  # Identify missing final outputs. Use -s rather than -f so a zero-byte file
  # left by an interrupted writer is regenerated instead of skipped.
  local -a missing_indices=()
  local index filename base expected_vtk
  for index in "${!frames[@]}"; do
    filename="$(basename "${frames[index]}")"
    base="${filename%.surf.gii}"
    expected_vtk="${output_dir}/${base}.vtk"
    if [[ "${OVERWRITE_FINAL_VTK}" -eq 1 || ! -s "${expected_vtk}" ]]; then
      missing_indices+=("${index}")
    fi
  done

  if [[ "${#missing_indices[@]}" -eq 0 ]]; then
    log "Skip completed ${H}/${surface_type}: all ${#frames[@]} final VTKs exist"
    return 0
  fi

  # Registration is sequential. To produce the last missing frame, only its
  # preceding chain is needed; completed frames after it can be ignored.
  local last_required_index
  last_required_index="${missing_indices[$((${#missing_indices[@]} - 1))]}"
  log "${H}/${surface_type}: ${#missing_indices[@]} final VTK(s) missing; process through frame $((last_required_index + 1))"

  local ref_base="${CASE_ID}.${HFS}.${surface_type}"
  local ref_surf="${work}/reference/${ref_base}.surf.gii"
  if [[ ! -f "${ref_surf}" ]]; then
    "${PYTHON_BIN}" "${VTK_TO_GII}" --vtk "${ref_vtk}" \
      --out "${work}/reference" --name "${ref_base}"
  fi
  set_struct_surface "${H}" "${ref_surf}"

  local ref_sphere="${work}/reference/sphere/${ref_base}.sphere.surf.gii"
  local ref_curv="${work}/reference/features/ref.curv.func.gii"
  make_sphere "${ref_surf}" "${H}" "${work}/reference/sphere" "${ref_base}"
  curv_metric "${H}" "${ref_surf}" "${ref_curv}"

  local ref_curv_sm="${work}/reference/features/ref.curv.sm.func.gii"
  smooth_metric "${H}" "${ref_sphere}" "${ref_curv}" "${ref_curv_sm}"
  local mask_tag="full"; [[ "${MASK_MEDIAL_WALL}" -eq 1 ]] && mask_tag="masked"
  local ref_curv_feature="${work}/reference/features/ref.curv.sm.${mask_tag}.func.gii"
  mask_feature "${H}" "${ref_curv_sm}" "${ref_vtk}" reference "${ref_curv_feature}"

  local -a required_frames=(
    "${frames[@]:0:$((last_required_index + 1))}"
  )
  log "Prepare ${#required_frames[@]} required ${H}/${surface_type} frames with ${PREP_JOBS_PER_SERIES} jobs"
  prepare_all_frames "${H}" "${sim_vtk_dir}" "${work}" "${required_frames[@]}"

  local previous_native="${ref_sphere}"
  local previous_cumulative="${ref_sphere}"
  local previous_sulc_feature=""
  for ((index = 0; index <= last_required_index; index++)); do
    local frame_surface="${frames[index]}" filename base tag frame_work
    filename="$(basename "${frame_surface}")"; base="${filename%.surf.gii}"
    tag="$(printf 'frame%02d' "${index}")"; frame_work="${work}/frames/${tag}"

    local native_sphere="${frame_work}/sphere/${base}.sphere.surf.gii"
    local curv="${frame_work}/features/native.curv.func.gii"
    local sulc="${frame_work}/features/native.sulc.func.gii"
    local curv_feature="${frame_work}/features/native.curv.sm.${mask_tag}.func.gii"
    local sulc_feature="${frame_work}/features/native.sulc.sm.${mask_tag}.func.gii"

    local in_feature ref_feature feature_name
    if [[ "${index}" -eq 0 ]]; then
      in_feature="${curv_feature}"
      ref_feature="${ref_curv_feature}"
      feature_name="curv"
    else
      in_feature="${sulc_feature}"
      ref_feature="${previous_sulc_feature}"
      feature_name="sulc"
    fi

    log "[${H}/${surface_type} $((index + 1))/${#frames[@]}] ${base}: MSM ${feature_name}"
    local prefix="${frame_work}/msm/${base}_to_previous.${feature_name}.${mask_tag}."
    run_msm "${native_sphere}" "${previous_native}" "${in_feature}" \
      "${ref_feature}" "${prefix}" "${feature_name}"

    local incremental="${prefix}sphere.reg.surf.gii" cumulative
    set_struct_sphere "${H}" "${incremental}"
    if [[ "${index}" -eq 0 ]]; then
      cumulative="${incremental}"
    else
      cumulative="${frame_work}/msm/${base}.cumulative_to_ref.${mask_tag}.sphere.surf.gii"
      compose_spheres "${H}" "${incremental}" "${previous_native}" \
        "${previous_cumulative}" "${cumulative}"
    fi

    local out_work="${work}/registered/${tag}"
    local reg_surface="${out_work}/${base}.registered.surf.gii"
    local reg_curv="${out_work}/${base}.curv.func.gii"
    local reg_sulc="${out_work}/${base}.sulc.func.gii"
    local out_vtk="${output_dir}/${base}.vtk"
    mkdir -p "${out_work}"

    if [[ "${OVERWRITE_FINAL_VTK}" -eq 0 && -s "${out_vtk}" ]]; then
      log "Skip existing final VTK: ${out_vtk}"
    else
      resample_surface "${frame_surface}" "${cumulative}" "${ref_sphere}" "${reg_surface}"
      set_struct_surface "${H}" "${reg_surface}"
      resample_metric "${H}" "${curv}" "${cumulative}" "${ref_sphere}" "${reg_curv}"
      resample_metric "${H}" "${sulc}" "${cumulative}" "${ref_sphere}" "${reg_sulc}"

      "${PYTHON_BIN}" "${WRITE_VTK}" --surf "${reg_surface}" \
        --scalar "Sulc=${reg_sulc}" --scalar "Curv=${reg_curv}" \
        --reference-vtk "${ref_vtk}" \
        --copy-scalar par_Clustering --copy-scalar par_huang \
        --copy-scalar par_FS2009 --out "${out_vtk}"
    fi

    previous_native="${native_sphere}"
    previous_cumulative="${cumulative}"
    previous_sulc_feature="${sulc_feature}"
  done
  log "Completed ${H}/${surface_type}: ${output_dir}"
}

# =============================================================================
# MAIN
# =============================================================================
ensure_environment
run_cutting_if_requested

if [[ "${PARALLEL_SERIES}" -eq 1 ]]; then
  pids=()
  series_names=()
  for H in "${HEMIS[@]}"; do
    for surface_type in "${SURFACE_TYPES[@]}"; do
      process_series "${H}" "${surface_type}" &
      pids+=("$!")
      series_names+=("${H}/${surface_type}")
    done
  done
  failed=0
  for index in "${!pids[@]}"; do
    if ! wait "${pids[index]}"; then
      log "[ERROR] Series failed: ${series_names[index]}"
      failed=1
    fi
  done
  [[ "${failed}" -eq 0 ]] || die "At least one anatomical series failed"
else
  for H in "${HEMIS[@]}"; do
    for surface_type in "${SURFACE_TYPES[@]}"; do
      process_series "${H}" "${surface_type}"
    done
  done
fi

log "DONE. Final legacy VTKs: ${OUT_ROOT}"
