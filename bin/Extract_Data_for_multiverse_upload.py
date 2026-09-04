#!/usr/bin/env python3
"""Collect selected Multiverse results in a BIDS-like upload directory."""

import argparse
import shutil
from pathlib import Path


DEFAULT_PROJECT_ROOT = Path(
    "/home/jcramer/Server_Aswendt/projects/Data_Cooperations/Multiverse/"
    "proc_data_good"
)
OUTPUT_DIRECTORY_NAME = "Extracted_data_for_upload"

# Paths are relative to PROJECT_ROOT. Add or remove patterns here to change the
# files included in an upload. Restricting the search to sub-*/ses-* also means
# that the Report directory is never considered.
FILE_PATTERNS = (
    (
        "sub-*/ses-*/func/**/*"
        "_registered_on_SIGMA_template_corrected_resampled_0p3mm.nii.gz"
    ),
    (
        "sub-*/ses-*/func/**/*"
        "_registered_on_SIGMA_template_temporal_mean_corrected_resampled_0p3mm.nii.gz"
    ),
)

# Each matched directory is copied recursively, including all files and nested
# directories. Its own name is retained, but parent directories below func or
# anat are flattened. Examples:
# DIRECTORY_PATTERNS = (
#     "sub-*/ses-*/func/**/regr",
#     "sub-*/ses-*/anat/**/t2_values_extraction",
# )
DIRECTORY_PATTERNS = ("sub-*/ses-*/func/**/*_bold_EPI_mcf.mat",)


def find_upload_files(project_root):
    """Return unique files matching the hard-coded upload patterns."""
    project_root = Path(project_root).expanduser().resolve()
    if not project_root.is_dir():
        raise NotADirectoryError(
            f"Project root does not exist or is not a directory: {project_root}"
        )

    matches = set()
    pattern_counts = {}
    for pattern in FILE_PATTERNS:
        pattern_matches = {
            path for path in project_root.glob(pattern) if path.is_file()
        }
        matches.update(pattern_matches)
        pattern_counts[("file", pattern)] = len(pattern_matches)

    return sorted(matches), pattern_counts


def find_upload_directories(project_root):
    """Return unique directories matching the hard-coded directory patterns."""
    project_root = Path(project_root).expanduser().resolve()
    if not project_root.is_dir():
        raise NotADirectoryError(
            f"Project root does not exist or is not a directory: {project_root}"
        )

    matches = set()
    pattern_counts = {}
    for pattern in DIRECTORY_PATTERNS:
        pattern_matches = {
            path for path in project_root.glob(pattern) if path.is_dir()
        }
        matches.update(pattern_matches)
        pattern_counts[("directory", pattern)] = len(pattern_matches)

    return sorted(matches), pattern_counts


def upload_relative_path(source_path, project_root):
    """Flatten a source path to sub-*/ses-*/anat|func/<filename>."""
    source_path = Path(source_path)
    project_root = Path(project_root)
    relative_path = source_path.relative_to(project_root)
    path_parts = relative_path.parts

    if (
        len(path_parts) < 4
        or not path_parts[0].startswith("sub-")
        or not path_parts[1].startswith("ses-")
        or path_parts[2] not in {"anat", "func"}
    ):
        raise ValueError(
            "Matched file is not below sub-*/ses-*/anat or sub-*/ses-*/func: "
            f"{source_path}"
        )

    return Path(*path_parts[:3], source_path.name)


def copy_upload_files(project_root, output_root=None, dry_run=False):
    """Copy selected files and directories into BIDS modality directories."""
    project_root = Path(project_root).expanduser().resolve()
    if output_root is None:
        output_root = project_root.parent / OUTPUT_DIRECTORY_NAME
    output_root = Path(output_root).expanduser().resolve()

    source_files, file_pattern_counts = find_upload_files(project_root)
    source_directories, directory_pattern_counts = find_upload_directories(
        project_root
    )
    pattern_counts = file_pattern_counts | directory_pattern_counts

    copy_plan = {}
    sources = [(source_path, "file") for source_path in source_files]
    sources.extend(
        (source_path, "directory") for source_path in source_directories
    )
    for source_path, source_kind in sources:
        relative_path = upload_relative_path(source_path, project_root)
        destination_path = output_root / relative_path
        previous_entry = copy_plan.get(destination_path)
        if previous_entry is not None and previous_entry[0] != source_path:
            raise FileExistsError(
                "Multiple source paths would be copied to the same destination: "
                f"{previous_entry[0]} and {source_path} -> {destination_path}"
            )
        copy_plan[destination_path] = (source_path, source_kind)

    for destination_path, (source_path, source_kind) in copy_plan.items():
        if not dry_run:
            if source_kind == "directory":
                shutil.copytree(
                    source_path,
                    destination_path,
                    dirs_exist_ok=True,
                    copy_function=shutil.copy2,
                )
            else:
                destination_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_path, destination_path)

    return sorted(copy_plan), pattern_counts


def _build_argument_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Copy selected Multiverse output files into a separate BIDS-like "
            "directory for upload."
        )
    )
    parser.add_argument(
        "-r",
        "--project-root",
        type=Path,
        default=DEFAULT_PROJECT_ROOT,
        help=f"Source project root (default: {DEFAULT_PROJECT_ROOT}).",
    )
    parser.add_argument(
        "-o",
        "--output-root",
        type=Path,
        help=(
            "Destination root. By default, Extracted_data_for_upload is "
            "created next to the project root."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List what would be copied without creating or changing files.",
    )
    return parser


def main(argv=None):
    parser = _build_argument_parser()
    args = parser.parse_args(argv)

    try:
        copied_files, pattern_counts = copy_upload_files(
            args.project_root,
            output_root=args.output_root,
            dry_run=args.dry_run,
        )
    except NotADirectoryError as error:
        parser.error(str(error))

    for (path_kind, pattern), count in pattern_counts.items():
        print(f"{count:>4} {path_kind}(s): {pattern}")

    action = "Would copy" if args.dry_run else "Copied"
    output_root = (
        args.output_root.expanduser().resolve()
        if args.output_root is not None
        else args.project_root.expanduser().resolve().parent
        / OUTPUT_DIRECTORY_NAME
    )
    print(f"{action} {len(copied_files)} path(s) to: {output_root}")

    if not copied_files:
        print("No files or directories matched the configured patterns.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
