# AIDAmri v3.0 Release Notes

## Acknowledgements

We would like to thank Paul Camacho from the Biomedical Imaging Center, University of Illinois Urbana-Champaign, for his valuable support and contributions to this release.


This release summarizes the major changes introduced in the Version 3.

## Highlights

- Updated the container workflow for Python 3.10, DSI Studio 2025, FSL
- Added improved support for animal brain extraction, bias field correction, DWI preprocessing, DSI Studio connectivity analysis, and incidence map generation.
- Introduced ANTs N4 bias correction, DIPY Patch2self, fslpy, and bet4animal.
- Improvement of batchProc.py by adding multiple new command-line arguments
- Improved NIfTI handling across the pipeline, including orientation checks

## Docker and Dependencies

- Base image changed from `ubuntu:18.04` to `ubuntu:22.04`.
- Added platform-aware Docker builds
- Updated Python setup to use the system Python 3 on Ubuntu 22.04.
- Updated `nipype` from `1.1.2` to `1.7.0`.
- Added ANTs 2.6.2 for `N4BiasFieldCorrection`.
- Updated DSI Studio to the Ubuntu 22.04 build from 2025.04.16.
- Added installation support for `immv` and `bet4animal`.

## Preprocessing

- Added stricter LIP orientation handling and header checks for T2 and DWI preprocessing.
- Added ANTs N4 bias field correction as an alternative to MICO.
- Added options to skip bias correction or BET where supported.
- Added bet4animal`support for animal brain extraction.
- Added user-adjustable BET parameters: horizontal gradient and center coordinates.

## DWI and DSI Studio

- Added b0 averaging support.
- Added Patch2Self denoising support.
- Added automatic bval/bvec repetition adjustment after conversion.
- Added DSI Studio reconstruction options for `dti` and `gqi`.
- Added `in_vivo` and `ex_vivo` parameter presets.
- Added isotropic resampling options
- Added tracking presets and support for custom tracking parameters.
- Added thread-count control for tractography.
- Added support for skipping slice-wise motion correction.
- Added legacy file type support for older DSI Studio output formats.
- Improved connectivity seed/ROI file detection and connectivity output naming.
- Added process logging improvements for direct and batch execution.

## Batch Processing

- Enhanced `batchProc.py` argument handling.
- Added grouped command-line options for batch control, CPU usage, T2 preprocessing, DWI preprocessing, animal BET settings, and DSI Studio processing.

- Added batch-level forwarding of new T2, DWI, BET, and DSI Studio options.

## Incidence and Stroke Mask Analysis

- Improved incidence map generation with session-specific input selection.
- Added validation for output names and session names.
- Incidence maps now use session/input-based output names.
- Heatmaps are now exported as PNG, PDF, and SVG.
- Improved affected-region and stroke-volume output handling.
- Added CSV output for affected region summaries.
- Added support for creating stroke mask label files from NIfTI images.
- Added `annotation_label_IDs.csv` for atlas structure identification.

## T2, T2 Map, fMRI, and ROI Processing

- Improved T2 value extraction output naming.
- Improved error handling when atlas or acronym files are missing.
- Renamed T2 value CSV headers to clarify that values are mean T2 values.
- Improved NIfTI compatibility in ROI, fMRI activity, and DTI data extraction scripts.

## Conversion and Helper Tools

- Updated `conv2Nifti_auto.py` so the default `proc_data` directory is created next to the raw data folder.
- Added automatic bval/bvec adjustment during conversion.
- Added QC image generation and HTML reports for converted NIfTI files.
- Added `plot_sourcedata_niftis.py` for source data QC mosaics.
- Added `crop_T2.py` for cropping T2-weighted images before preprocessing.
- Added `fieldmap_json_edit.py` for editing BIDS `IntendedFor` fields.
- Added a helper script for cleaning subject list text files by removing spaces and caret characters.
- Added an `immv` installation helper so `bet4animal` can use the FSL image-moving utility inside the Docker container.

## Documentation and Figures

- Added a new Markdown manual for AIDAmri v3.0.
- Updated Docker, preprocessing, registration, T2, DTI, fMRI, and ROI analysis documentation.
- Added SVG versions of documentation figures.
- Removed the committed `manual.pdf`.
- Updated README wording and Docker/reorientation notes.

## Data and Atlas Files

- Updated several atlas/template files in `lib/`.
- Added `lib/annotation_label_IDs.csv`.
- Removed deprecated `lib/DTI_Jones30.txt`.
- Removed deprecated `lib/MPI_maskBIG_for_incidence2.nii`.