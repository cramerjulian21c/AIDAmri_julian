"""Run DTIdata_extract.py across proc_folder/subject/session/dwi folders.

Match the ROI name passed with -r immediately before .nii or .nii.gz,
ignoring case (e.g. annosplit does not match AnnoSplit_parental).
For each session, run the extractor sequentially on the selected maps
in DSI_studio. Select FA, MD, RD or AD with -m/--metric. The default
selection, all, processes all four metrics and legacy .fa0 maps.
Forward an optional -t translation table without choosing one automatically.

The extractor creates individual CSV files beside the parameter maps.
Show warnings and errors in the terminal; the extractor writes its own log.
Suppress normal extractor output such as Done messages during batch runs.
List failed and fully or partially skipped sessions in the final summary.
The skipped counter counts skipped entries, which may be whole sessions
or individual missing metrics, rather than unique skipped sessions.
Continue after individual failures. Exit with status 1 if any extraction
or session fails, or no extraction succeeds; otherwise exit with status 0.

Usage:
    python batch_DTIdata_extract.py -i PROC_FOLDER -r annosplit_parental [-t TXT] [-m fa]

Created on 08.04.2019 by Niklas Pallast and Markus Aswendt.
"""

import argparse
from pathlib import Path
import subprocess
import sys
from typing import NamedTuple


METRICS = ("fa", "md", "rd", "ad")


class SessionResult(NamedTuple):
    succeeded: int
    skipped: int
    failed: int
    skip_reason: str = ""


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


def _print_summary(succeeded, skipped, failed, failed_sessions=None, reason="", skipped_sessions=None):
    message = f"Summary: {succeeded} succeeded, {skipped} skipped, {failed} failed."
    print(f"{message} {reason}" if reason else message, flush=True)
    if failed_sessions:
        print("Failed sessions:")
        for session, note in sorted(failed_sessions.items()):
            print(f"  {session}" + (f" ({note})" if note else ""))
    if skipped_sessions:
        print("Skipped sessions:")
        for session, note in sorted(skipped_sessions.items()):
            print(f"  {session} ({note})")
    return 1 if failed or not succeeded else 0


def _run_session(dwi, context, roi_name, metric_selection, translator, extractor):
    succeeded = skipped = failed = 0
    missing_metrics = []
    try:
        rois = _find_rois(dwi, roi_name)
        if not rois:
            print(f"Warning: [{context}] No ROI matching '*{roi_name}.nii[.gz]'. Skipping session.",
                  file=sys.stderr)
            return SessionResult(0, 1, 0, f"fully skipped: ROI matching '*{roi_name}.nii[.gz]' missing")
        if len(rois) > 1:
            print(f"Error: [{context}] Multiple matching ROIs: {', '.join(p.name for p in rois)}. "
                  "Skipping session.", file=sys.stderr)
            return SessionResult(0, 0, 1)
        metrics_folder = dwi / "DSI_studio"
        if not metrics_folder.is_dir():
            print(f"Warning: [{context}] DSI_studio folder missing. Skipping session.", file=sys.stderr)
            return SessionResult(0, 1, 0, "fully skipped: DSI_studio folder missing")
        metrics = _find_metrics(metrics_folder, metric_selection)
    except OSError as error:
        print(f"Error: [{context}] Cannot read session inputs: {error}", file=sys.stderr)
        return SessionResult(0, 0, 1)

    # Different NIfTI encodings of one image would target the same CSV.
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
            missing_metrics.append(metric.upper())
        for image in images:
            if image in collisions:
                continue
            command = [sys.executable, str(extractor), "-i", str(image), "-r", str(rois[0])]
            if translator is not None:
                command.extend(("-t", str(translator)))
            try:
                # The extractor owns its log. Keep warnings/errors visible while
                # suppressing its normal stdout only for batch execution.
                result = subprocess.run(command, cwd=dwi, check=False, stdout=subprocess.DEVNULL, stderr=None)
            except OSError as error:
                print(f"Error: [{context}] Could not run extractor for '{image.name}': {error}", file=sys.stderr)
                failed += 1
                continue
            if result.returncode != 0:
                print(f"Error: [{context}] Extraction failed for '{image.name}' "
                      f"(exit code {result.returncode}).", file=sys.stderr)
                failed += 1
            else:
                succeeded += 1
    skip_reason = ""
    if missing_metrics:
        status = "fully" if succeeded == 0 and failed == 0 else "partially"
        skip_reason = f"{status} skipped: {', '.join(missing_metrics)} missing"
    return SessionResult(succeeded, skipped, failed, skip_reason)


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
        return _print_summary(0, 0, 1, reason=f"Error: Proc folder does not exist or is not a directory: '{root}'.")
    if not roi_name:
        return _print_summary(0, 0, 1, reason="Error: ROI name must not be empty.")
    for path in (extractor, translator):
        if path is not None and not path.is_file():
            return _print_summary(0, 0, 1, reason=f"Error: Required file does not exist: '{path}'.")

    succeeded = skipped = failed = 0
    failed_sessions = {}
    skipped_sessions = {}
    try:
        sessions = findData(root)
    except OSError as error:
        return _print_summary(0, 0, 1, reason=f"Error: Cannot search '{root}': {error}")
    if not sessions:
        return _print_summary(0, 0, 0, reason=f"No subject/session/dwi folders found in '{root}'.")

    for dwi in sessions:
        context = str(dwi.relative_to(root))
        session = f"{dwi.parent.parent.name} / {dwi.parent.name}"
        result = _run_session(dwi, context, roi_name, args.metric, translator, extractor)
        succeeded += result.succeeded
        skipped += result.skipped
        failed += result.failed
        if result.failed:
            failed_sessions.setdefault(session, "")
        if result.skip_reason:
            skipped_sessions[session] = result.skip_reason

    return _print_summary(succeeded, skipped, failed, failed_sessions, skipped_sessions=skipped_sessions)


if __name__ == "__main__":
    sys.exit(main())
