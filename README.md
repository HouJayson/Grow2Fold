# Grow2Fold

## Mapping heterogeneous developmental growth to human brain folding

Research code accompanying the manuscript **Grow2Fold: mapping heterogeneous developmental growth to human brain folding**.

**Authors:** Jixin Hou, Akbar Solhtalab, Kun Jiang, Beikang Gu, Tiantain Li, Xianyan Chen, Taotao Wu, Kenan Song, Yang Liu, Wenzhan Song, Dajiang Zhu, Tianming Liu, Ellen Kuhl, Gang Li, Mir Jalil Razavi, and Xianqiao Wang.

**Corresponding author:** Xianqiao Wang.

## Overview

Grow2Fold investigates how spatially heterogeneous developmental growth contributes to prenatal human cortical folding. The framework combines anatomically realistic brain geometry, cortical regionalization, and finite-element growth simulations in Abaqus/Explicit.

This repository contains the cortical regionalization script, tools for preparing labeled hexahedral finite-element models, a Fortran subroutine for region-specific orthotropic growth, and scripts for extracting and processing simulated cortical surfaces.

This is the code release prepared for manuscript submission. Subject-specific inputs and several supporting files are not included in this snapshot, so reproducing the complete workflow requires additional materials and configuration described below. Symbolic-regression and manuscript figure-generation scripts are not included in this directory.

## Repository structure

```text
.
├── README.md
├── Cortical parcellation atlas/
│   └── Cortical_Regionalization.m
└── abaqus_script/
    ├── meshing/
    │   ├── SurfToHex_whole.m
    │   └── Functions/                       # Mesh preparation and export utilities
    ├── modeling/
    │   └── vuexpan_orth_wholebrain.f
    └── postprocessing/
        ├── ExtractPointsFromABAQUS.py
        ├── surfReconstruction.m
        ├── hemiBrain_manual.py
        ├── Registration_hemi.sh
        ├── vtk_to_gii_metrics.py
        └── write_metrics_vtk.py
```

## Software requirements

Requirements depend on the stage being run:

| Stage | Required software |
| --- | --- |
| Cortical regionalization and model preparation | MATLAB; Statistics and Machine Learning Toolbox for functions such as `kmeans`, `nnmf`, and `knnsearch`; Parallel Computing Toolbox for the supplied `parfor` loops |
| Finite-element simulation | Licensed Abaqus/Explicit and a Fortran compiler supported by the installed Abaqus release for compiling `VUEXPAN` |
| Abaqus output extraction | Abaqus Python with `odbAccess`, NumPy, and SciPy available in that environment |
| Surface reconstruction and folding metrics | MATLAB, the relevant toolboxes above, and the additional MATLAB helpers listed below |
| Hemisphere separation and format conversion | Python 3 with NumPy, PyVista, VTK, PyMeshLab, and NiBabel |
| Surface registration | Bash on Linux or WSL, FreeSurfer, Connectome Workbench (`wb_command`), MSM or newMSM, and the Python surface-processing packages |

Exact tested software versions are not recorded in this snapshot. Use an Abaqus-compatible compiler and Python environment; the standalone Python environment does not provide `odbAccess`.

For the standalone Python surface-processing scripts, install the imported packages in a dedicated environment:

```bash
python -m pip install numpy scipy pyvista vtk pymeshlab nibabel
```

This command does not install MATLAB, Abaqus, FreeSurfer, Connectome Workbench, or MSM/newMSM.

## Workflow

The commands below are usage templates. Configure the paths, subject identifiers, inputs, and stage-specific settings before running them.

### 1. Cortical regionalization

`Cortical parcellation atlas/Cortical_Regionalization.m` constructs cortical partitions from vertex-area trajectories across gestational ages. It includes spectral clustering, NMF, orthogonal NMF, and graph-regularized orthogonal NMF backends, with configurable model-selection and stability analyses.

- Supply left and right hemisphere VTK surfaces with consistent vertex correspondence across scans and the required surface label arrays.
- Configure `cfg.data_root`, `cfg.pattern`, `cfg.inflated`, and the output paths, including `cfg.spec_save_similarity_dir`.
- Set the age filters, clustering backend, candidate region counts (`cfg.k_range`), and selection settings for the intended analysis. The distributed defaults are not a complete manuscript reproduction configuration.
- Add `abaqus_script/meshing/Functions` to the MATLAB path to make `mvtk_read` available, then run `Cortical_Regionalization` from its directory.

Outputs include cortical label maps and model-selection results in the configured output directory. Regionalization is a separate analysis stage; mapping its labels to the subject-specific simulation surfaces must be prepared before model generation.

### 2. Prepare the finite-element model

`abaqus_script/meshing/SurfToHex_whole.m` reads an existing hexahedral brain mesh and reference surfaces, assigns regional labels and material orientations, constructs boundary sets, and exports Abaqus input and surface-connectivity files. An initial hexahedral mesh must be supplied; this script does not generate that mesh from MRI alone.

Run MATLAB from `abaqus_script/meshing` so that the relative `Functions` path resolves. Configure `case_idx`, `sessionId_min`, `ga_min`, and the input/output directories. Replace the literal `*` session placeholder with the actual identifier.

The script expects inputs under `inputFiles/longitudinal/GW21/<case_idx>/` by default:

```text
sub-<case_idx>_ses-<sessionId_min>.GA<ga_min>.W.gray_hull_label.vtk
sub-<case_idx>_ses-<sessionId_min>.GA<ga_min>.W.white_hull_reg_label.vtk
sub-<case_idx>_ses-<sessionId_min>.GA<ga_min>.slice.vtk
sub-<case_idx>_ses-<sessionId_min>.GA<ga_min>.inner_skull.vtk
Brain_<case_idx>_Hex.inp
```

The reference surfaces must provide the label arrays accessed by the script, including `par_Clustering`, `par_huang`, `par_FS2009`, and `par_MMP` where used. Check mesh element-set names and subject-specific smoothing/growth-indicator settings against the supplied inputs.

```matlab
SurfToHex_whole
```

Outputs include a configured `.inp` model and `Sconn_<case_idx>.xlsx` in `dataStorage/<case_idx>/`, plus a VTK mesh for inspection in `temp/<case_idx>/`.

### 3. Run the growth simulation

`abaqus_script/modeling/vuexpan_orth_wholebrain.f` implements region-specific orthotropic growth through Abaqus's `VUEXPAN` interface. Regional material names in the input model must match those dispatched by the subroutine.

From a working directory containing the configured input model and a copy of the subroutine, a typical submission command is:

```bash
abaqus job=brain_growth input=model.inp user=vuexpan_orth_wholebrain.f interactive
```

Replace `model.inp` with the generated input filename. Configure solver resources and the compiler environment for your system. Simulation output is written to an Abaqus `.odb` database; runtime and memory requirements depend on mesh size, solver settings, and hardware.

### 4. Extract simulated surface coordinates

Edit the user settings in `abaqus_script/postprocessing/ExtractPointsFromABAQUS.py`, including `odb_path`, `file_path`, `instance_name`, `user_id`, step names, and frame counts. Despite its placeholder name (`path-to-excel`), `file_path` is the destination directory for exported MATLAB `.mat` files. Configure the additional Python package path only if needed by your Abaqus environment.

Run through Abaqus Python:

```bash
abaqus python ExtractPointsFromABAQUS.py
```

The script exports surface coordinates and displacement magnitudes, node labels, and step/frame metadata. Its default two-step frame plan exports 26 frames and skips the duplicate initial frame of the second step. Adapt this plan to the actual database.

### 5. Reconstruct surfaces and compute folding metrics

`abaqus_script/postprocessing/surfReconstruction.m` combines extracted coordinates with the connectivity workbook and reference surfaces. It reconstructs gray and white matter surfaces and computes cortical thickness, curvature, gyrification, and sulcal-depth measures.

Configure the reference, connectivity, input, and output paths, case identifiers, filenames, and frame counts before running it. Ensure that extracted node ordering matches the connectivity files.

**Additional preparation is required in this snapshot:**

- `caseId_ref` is referenced but not defined; define the corresponding reference-case identifiers.
- The script points to an external `functionsMat` directory. Helpers such as `quadToTriConnectivity`, `corticalThickness`, `dataSmoothing`, `Curvatures`, `GyrificationIndex`, `SulcalDepth`, and the alpha-shape variants are not included here and must be supplied on the MATLAB path.
- The supplied meshing utilities provide `mvtk_read`, `mvtk_write`, and `laplacianSmoothing`; add their directory to the MATLAB path as needed.

### 6. Separate hemispheres and register surfaces

`hemiBrain_manual.py` splits whole-brain surfaces using a plane estimated from a reference slice, caps the cut boundaries, and remeshes the hemispheres. Configure `CASE_DIR`, `SLICE_FILE`, axis conventions, frame counts, and output settings, then run:

```bash
python hemiBrain_manual.py
```

`vtk_to_gii_metrics.py` converts VTK surfaces and selected metrics to GIFTI. `write_metrics_vtk.py` writes GIFTI surface geometry and metrics back to VTK. Inspect their command-line options with:

```bash
python vtk_to_gii_metrics.py --help
python write_metrics_vtk.py --help
```

`Registration_hemi.sh` performs sequential registration of left/right gray and white surface series. It aligns the initial frame using curvature, then registers later frames to the preceding frame using sulcal depth and composes the mappings into reference sphere space.

Before running registration, configure the case paths, reference surfaces, executable locations, and frame count. Supply the selected configuration file (`configs/msm_config_init` or `configs/newmsm_config_init`), which is absent from this snapshot. If `PERFORM_CUT=1`, the script also expects `hemiBrain_fast_local_cut.py`, which is not included. The included `hemiBrain_manual.py` has its own settings-based interface; it is not a direct replacement for that command-line helper. With `PERFORM_CUT=0`, prepare hemisphere outputs in the paths expected by the registration script first.

Once these prerequisites are satisfied:

```bash
bash Registration_hemi.sh
```

## Data and reproducibility

This snapshot contains source code. It does not include prenatal MRI data, subject-specific reference surfaces, initial hexahedral meshes, complete simulation inputs or outputs, registration reference/configuration files, or all MATLAB metric helpers. Input-data access and redistribution remain subject to the original providers' terms.

The study uses prenatal MRI data from the Developing Human Connectome Project (dHCP), as described in the manuscript. Consult the manuscript's data-availability statement for the dataset, access procedure, and any supplementary materials. No supplementary-data download link or archived release DOI is provided in this snapshot.

To reproduce a particular experiment, use the matching subject inputs, regional labels, growth laws, solver settings, frame plan, and analysis configuration. A complete end-to-end run has not been verified for this release.

## Citation

Please cite the accompanying manuscript when using this work:

> Hou, J., et al. *Grow2Fold: mapping heterogeneous developmental growth to human brain folding.* Manuscript prepared for submission.

Publication details and a DOI can be added when available.

## License

The source headers refer to a repository `LICENSE` file, but no such file is included in this snapshot. A license must be selected and added by the authors to specify terms for use, modification, and redistribution. Retain applicable notices and license terms for third-party utilities.
