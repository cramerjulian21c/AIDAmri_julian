# AIDAmri Helper Tools

This folder contains utility scripts for data preparation, quality control, file cleanup, and batch support tasks used in AIDAmri workflows.

## Overview

| Script | Purpose |
| --- | --- |
| `DistributeStrokeMasks.py` | Resample and propagate existing stroke masks across timepoints with NiftyReg `reg_resample`. |
| `adjustbvecRep.py` | Adjust repeated DWI `.bval` and `.bvec` sidecars after conversion. |
| `MRI_files_summarizer.py` | Create a CSV inventory of NIfTI files under `**/brkraw/*.nii.gz`. |
| `ReorientBatch.py` | Reorient NIfTI files to a target orientation, defaulting to AIDAmri's LIP orientation. |
| `crop_T2.py` | Crop T2-weighted images to a defined field of view before preprocessing. |
| `fieldmap_json_edit.py` | Add or update BIDS `IntendedFor` fields in fieldmap JSON files. |
| `getAtlasRegionSize_BIDS.py` | Compute region-wise atlas volumes for BIDS-style folder structures. |
| `getAtlasRegionSize_noBIDS.py` | Compute region-wise atlas volumes for non-BIDS folder structures. |
| `plot_sourcedata_niftis.py` | Generate QC mosaic images and an HTML report for source NIfTI files. |
| `remove_carets_spaces.sh` | Clean plain-text `subject` files by removing spaces and caret (`^`) characters. |
| `reset_naming.py` | Clean Bruker `subject` files before PV-to-NIfTI conversion. |

## Atlas Region Size

### `getAtlasRegionSize_BIDS.py`

Computes region-wise volumes from registered annotation NIfTI files in BIDS-style datasets. The script searches recursively under the input folder and writes one `.txt` and one `.mat` output per annotation file.

Supported annotation patterns:

- `**/*_Anno_parental.nii.gz`
- `**/*_AnnoSplit_parental.nii.gz`
- `**/*_AnnorsfMRI.nii.gz`
- `**/*_Anno.nii.gz`
- `**/*_AnnoSplit.nii.gz`

Outputs are written to:

```text
<inputFolder>/output_region_size/
```

Usage:

```bash
python getAtlasRegionSize_BIDS.py -i /path/to/input_folder
```

### `getAtlasRegionSize_noBIDS.py`

Computes region-wise volumes for non-BIDS workflows. It expects an input folder containing registered annotation files and paired masks.

Usage:

```bash
python getAtlasRegionSize_noBIDS.py -i /path/to/T2w_folder
```

Both atlas region scripts report voxel counts and physical volumes in `mm^3`. The volume calculation uses the configured voxel dimensions written in the output headers.

## Bruker Subject File Cleanup

### `reset_naming.py`

Prepares Bruker raw data before PV-to-NIfTI conversion. The script scans recursively for files named `subject` and updates them in place.

Main operations:

- Removes the first underscore in `SUBJECT_id` and `SUBJECT_study_name` lines.
- Replaces `baseline` with `PT0` in study names.

Usage:

```bash
python reset_naming.py -i /path/to/raw_data
```

### `remove_carets_spaces.sh`

Removes spaces and caret (`^`) characters from each line of a plain-text input file, typically a Bruker `subject` file.

Usage:

```bash
bash remove_carets_spaces.sh /path/to/subject
```

Outputs:

- `subject_clean.txt`: cleaned lines with spaces and caret characters removed.
- `subject_orig.txt`: copy of the original input file.

Note: despite the historical script comments, the current implementation does not remove underscores.

## Reorientation

### `ReorientBatch.py`

Batch-reorients `.nii` and `.nii.gz` files under an input root while mirroring the folder structure to an output root.

Key behavior:

- Reorients NIfTI images to a target orientation, defaulting to `LIP`.
- Copies non-NIfTI files unchanged.
- Copies `.bval` sidecars.
- Reorients `.bvec` sidecars consistently with the image orientation change.

Usage:

```bash
python ReorientBatch.py -i /path/to/input_root -o /path/to/output_root
```

Common options:

- `-t`, `--target`: target orientation, for example `LIP`.
- `-n`, `--non_interactive`: run without prompts.
- `-l`, `--logfile`: log file name written into the output root.

## Quality Control and Data Summaries

### `adjustbvecRep.py`

Adjusts `.bval` and `.bvec` files for DWI acquisitions where the gradient table needs to be repeated to match the number of image volumes. Files with unusable volume counts are moved into hidden subfolders such as `.single_volume` or `.low_volume_count`.

Usage:

```bash
python adjustbvecRep.py /path/to/sub-001/ses-001
```

This helper is called automatically by `conv2Nifti_auto.py` during conversion when DWI data are present.

### `MRI_files_summarizer.py`

Creates a CSV overview of NIfTI files under `**/brkraw/*.nii.gz`. Filename tokens are used to extract subject, session, run, and modality information.

Usage:

```bash
python MRI_files_summarizer.py -i /path/to/project -o /path/to/output_folder
```

Output:

```text
MRI_files_overview.csv
```

### `plot_sourcedata_niftis.py`

Creates QC mosaic PNG files and an HTML report for NIfTI files in a subject/session folder.

Usage:

```bash
python plot_sourcedata_niftis.py /path/to/sub-001/ses-001 /path/to/qc_output
```

Optional argument:

- `--n_slices`: number of slices per orientation, default `10`.

## T2 Cropping

### `crop_T2.py`

Crops a T2-weighted image to a defined field of view. The script uses FSL `fslroi` through Nipype and writes quick-look PNG images with FSL `slicer` when available.

Usage:

```bash
python crop_T2.py \
  -i /path/to/input_T2w.nii.gz \
  -x_max 120 \
  -y_max 120 \
  -z_max 40 \
  -o /path/to/cropped_T2w.nii.gz
```

Optional lower bounds:

- `-x_min`, default `0`
- `-y_min`, default `0`
- `-z_min`, default `0`

## Fieldmap JSON Editing

### `fieldmap_json_edit.py`

Updates BIDS fieldmap JSON files with `IntendedFor` entries for matching DWI and functional files.

Usage:

```bash
python fieldmap_json_edit.py \
  --bids_root /path/to/proc_data \
  --participant 001 \
  --session PT0
```

Optional argument:

- `--overwrite`: replace existing `IntendedFor` fields.

## Stroke Mask Distribution

### `DistributeStrokeMasks.py`

Finds existing stroke masks and propagates them to other timepoints by resampling with NiftyReg.

Search pattern:

```text
**/anat/*Stroke_mask.nii.gz
```

Workflow:

1. Resample each stroke mask into incidence space using `*IncidenceData.nii.gz` and `*MatrixInv.txt`.
2. Resample from incidence space into other timepoints using `*MatrixBspline.nii` and a `*BiasBet.nii.gz` reference.
3. Write `missing_files_log.txt` when required files are missing.

Usage:

```bash
python DistributeStrokeMasks.py -i /path/to/subject_root_or_dataset_root
```

Outputs:

- New `*Stroke_mask.nii.gz` files in timepoints that do not already have a stroke mask.
- `missing_files_log.txt` in the input root.

Requirement:

- `reg_resample` must be available on `PATH`.
