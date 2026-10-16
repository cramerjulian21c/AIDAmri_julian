"""Run DTIdata_extract.py across proc_folder/subject/session/dwi folders.

Match the ROI name passed with -r immediately before .nii or .nii.gz,
ignoring case (e.g. annosplit does not match AnnoSplit_parental).
For each session, run the extractor sequentially on the selected maps
in DSI_studio. Select FA, MD, RD or AD with -m/--metric. The default
selection, all, processes all four metrics and legacy .fa0 maps.
Forward an optional -t translation table without choosing one automatically.

The extractor creates individual CSV files beside the parameter maps.
Report missing inputs, ambiguous ROIs and colliding output names, and
continue after individual failures. Exit with status 1 if any extraction
or session fails, or no extraction succeeds; otherwise exit with status 0.

Usage:
    python batch_DTIdata_extract.py -i PROC_FOLDER -r annosplit_parental [-t TXT] [-m fa]

Created on 08.04.2019 by Niklas Pallast and Markus Aswendt.
"""

import argparse
from pathlib import Path
import subprocess
import sys


METRICS = ("fa", "md", "rd", "ad")


def _nifti_stem(path):
    for suffix in (".nii.gz", ".nii"):
        if path.name.lower().endswith(suffix):
            return path.name[:-len(suffix)]
    return None


def findData(path):
    """Return sorted session DWI directories without assuming naming prefixes."""
    return sorted(folder for folder in Path(path).glob("*/*/dwi") if folder.is_dir())


def _find_rois(dwi, roi_name):
    return sorted(
        path for path in dwi.iterdir()
        if path.is_file() and (stem := _nifti_stem(path)) is not None
        and stem.lower().endswith(roi_name.lower())
    )


def _find_metrics(folder, selected_metric="all"):
    selected = METRICS if selected_metric == "all" else (selected_metric,)
    metrics = {metric: [] for metric in selected}
    for path in sorted(folder.iterdir()):
        stem = _nifti_stem(path)
        if not path.is_file() or stem is None:
            continue
        metric = stem.rsplit(".", 1)[-1].lower()
        if metric == "fa0" and selected_metric == "all":
            metric = "fa"
        if metric in metrics:
            metrics[metric].append(path)
    return metrics


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("-i", "--input", required=True,
                        help="Cohort root containing subject/session/dwi folders")
    parser.add_argument("-r", "--roi_name", required=True,
                        help="ROI filename suffix before .nii[.gz], e.g. annosplit or annosplit_parental")
    parser.add_argument("-t", "--translatorTXT", help="Optional translation TXT passed to DTIdata_extract.py")
    parser.add_argument("-m", "--metric", choices=(*METRICS, "all"), default="all", type=str.lower,
                        help="Metric suffix before .nii[.gz] (default: all; includes legacy .fa0 maps)")
    args = parser.parse_args(argv)

    root = Path(args.input).expanduser().resolve()
    translator = Path(args.translatorTXT).expanduser().resolve() if args.translatorTXT is not None else None
    extractor = Path(__file__).resolve().with_name("DTIdata_extract.py")
    roi_name = args.roi_name.strip()
    if not root.is_dir():
        print(f"Error: Proc folder does not exist or is not a directory: '{root}'.", file=sys.stderr)
        return 1
    if not roi_name:
        print("Error: ROI name must not be empty.", file=sys.stderr)
        return 1
    for path in (extractor, translator):
        if path is not None and not path.is_file():
            print(f"Error: Required file does not exist: '{path}'.", file=sys.stderr)
            return 1

    succeeded = skipped = failed = 0
    try:
        sessions = findData(root)
    except OSError as error:
        print(f"Error: Cannot search '{root}': {error}", file=sys.stderr)
        return 1
    if not sessions:
        print(f"Warning: No subject/session/dwi folders found in '{root}'.", file=sys.stderr)

    for dwi in sessions:
        context = str(dwi.relative_to(root))
        try:
            rois = _find_rois(dwi, roi_name)
            if not rois:
                print(f"Warning: [{context}] No ROI matching '*{roi_name}.nii[.gz]'. Skipping session.",
                      file=sys.stderr)
                skipped += 1
                continue
            if len(rois) > 1:
                print(f"Error: [{context}] Multiple matching ROIs: {', '.join(p.name for p in rois)}. "
                      "Skipping session.", file=sys.stderr)
                failed += 1
                continue
            metrics_folder = dwi / "DSI_studio"
            if not metrics_folder.is_dir():
                print(f"Warning: [{context}] DSI_studio folder missing. Skipping session.", file=sys.stderr)
                skipped += 1
                continue
            metrics = _find_metrics(metrics_folder, args.metric)
        except OSError as error:
            print(f"Error: [{context}] Cannot read session inputs: {error}", file=sys.stderr)
            failed += 1
            continue

        # Different NIfTI encodings of one image would target the same CSV.
        # This should be avoided, so check for collisions and skip all colliding maps.
        by_stem = {}
        for images in metrics.values():
            for image in images:
                by_stem.setdefault(_nifti_stem(image), []).append(image)
        collisions = set()
        for stem, images in by_stem.items():
            if len(images) > 1:
                print(f"Error: [{context}] Colliding CSV outputs for '{stem}': "
                      f"{', '.join(p.name for p in images)}. Skipping these maps.", file=sys.stderr)
                collisions.update(images)
                failed += len(images)

        for metric, images in metrics.items():
            if not images:
                print(f"Warning: [{context}] No {metric.upper()} map found. Skipping parameter.", file=sys.stderr)
                skipped += 1
            for image in images:
                if image in collisions:
                    continue
                print(f"[{context}] {metric.upper()}: {image.name}; ROI: {rois[0].name}", flush=True)
                command = [sys.executable, str(extractor), "-i", str(image), "-r", str(rois[0])]
                if translator is not None:
                    command.extend(("-t", str(translator)))
                try:
                    result = subprocess.run(command, cwd=dwi, check=False)
                except OSError as error:
                    print(f"Error: [{context}] Could not run extractor for '{image.name}': {error}",
                          file=sys.stderr)
                    failed += 1
                    continue
                if result.returncode != 0:
                    print(f"Error: [{context}] Extraction failed for '{image.name}' "
                          f"(exit code {result.returncode}).", file=sys.stderr)
                    failed += 1
                else:
                    succeeded += 1

    print(f"Summary: {succeeded} succeeded, {skipped} skipped, {failed} failed.", flush=True)
    return 1 if failed or not succeeded else 0


if __name__ == "__main__":
    sys.exit(main())
