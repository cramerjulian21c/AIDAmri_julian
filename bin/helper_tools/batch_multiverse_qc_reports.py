"""Create registration and corpus-callosum QC reports for Multiverse output."""

import argparse
import logging
from pathlib import Path

import nibabel as nib
import numpy as np

try:
    from .batch_qc_reports import (
        CC_ATLAS_LABELS,
        _custom_parameter,
        _entry_metadata,
        _plot_registration_overlay,
        _positive_int,
        _write_report,
    )
except ImportError:  # Direct execution from bin/helper_tools.
    from batch_qc_reports import (
        CC_ATLAS_LABELS,
        _custom_parameter,
        _entry_metadata,
        _plot_registration_overlay,
        _positive_int,
        _write_report,
    )


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SIGMA_ATLAS = (
    REPO_ROOT / "lib" / "sigma" / "SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz"
)
CORRECTED_MEAN_PATTERN = (
    "sub-*/ses-*/func/rs-fMRI_niiData/Multiverse_Output/"
    "*SIGMA_template_temporal_mean_corrected.nii.gz"
)


def _validate_overlay_geometry(mean_path, atlas_path):
    """Reject overlays whose voxel arrays or spatial affines do not match."""
    mean_img = nib.load(str(mean_path))
    atlas_img = nib.load(str(atlas_path))

    mean_shape = tuple(mean_img.shape[:3])
    atlas_shape = tuple(atlas_img.shape[:3])
    if mean_img.ndim != 3:
        raise ValueError(
            f"Expected a 3D corrected temporal mean, got shape {mean_img.shape}"
        )
    if atlas_img.ndim != 3:
        raise ValueError(f"Expected a 3D SIGMA atlas, got shape {atlas_img.shape}")
    if mean_shape != atlas_shape:
        raise ValueError(
            "Corrected temporal mean and SIGMA atlas have different shapes: "
            f"{mean_shape} vs {atlas_shape}"
        )
    if not np.allclose(mean_img.affine, atlas_img.affine, rtol=0, atol=1e-4):
        raise ValueError(
            "Corrected temporal mean and SIGMA atlas have different spatial affines"
        )


def _build_multiverse_report(
    project_dir,
    atlas_path,
    report_kind,
    n_slices=10,
    custom_parameters=None,
):
    project_dir = Path(project_dir)
    atlas_path = Path(atlas_path)
    if not atlas_path.is_file():
        raise FileNotFoundError(f"SIGMA atlas does not exist: {atlas_path}")

    is_cc = report_kind == "cc"
    report_label = "CC" if is_cc else "Registration"
    out_dir = project_dir / "Report" / "Multiverse" / report_label
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    mean_files = sorted(project_dir.glob(CORRECTED_MEAN_PATTERN))
    for mean_path in mean_files:
        try:
            _validate_overlay_geometry(mean_path, atlas_path)
            png_path, mean_shape, atlas_shape, zooms = _plot_registration_overlay(
                mean_path,
                atlas_path,
                out_dir,
                project_dir,
                n_slices,
                atlas_labels=CC_ATLAS_LABELS if is_cc else None,
                report_title=(
                    "Multiverse Corpus Callosum Report"
                    if is_cc
                    else "Multiverse Registration Report"
                ),
                filename_suffix=(
                    "multiverse_cc_report"
                    if is_cc
                    else "multiverse_registration_report"
                ),
                png_source_path=mean_path,
            )
            relative_mean, subject, session, modality = _entry_metadata(
                mean_path,
                project_dir,
            )
            info = [
                ("Corrected temporal mean", relative_mean),
                ("SIGMA atlas", atlas_path),
            ]
            if is_cc:
                info.append(
                    (
                        "Atlas labels",
                        ", ".join(str(label) for label in CC_ATLAS_LABELS),
                    )
                )
            info.extend(
                [
                    ("Modality", modality),
                    ("Mean dimensions", mean_shape),
                    ("Atlas dimensions", atlas_shape),
                    ("Voxel size", tuple(round(float(z), 4) for z in zooms)),
                ]
            )
            entries.append(
                {
                    "subject": subject,
                    "session": session,
                    "modality": modality,
                    "report_img_path": png_path.name,
                    "image_alt": f"{mean_path.name} + {atlas_path.name}",
                    "info": info,
                }
            )
        except Exception as exc:
            logging.warning(
                "Could not create Multiverse %s report for %s: %s",
                report_label,
                mean_path,
                exc,
            )

    if not entries:
        return None, 0

    if is_cc:
        title = "Multiverse Corpus Callosum Report: corrected mean + SIGMA atlas"
        report_name = "multiverse_cc_report.html"
    else:
        title = "Multiverse Registration Report: corrected mean + SIGMA atlas"
        report_name = "multiverse_registration_report.html"

    html_path = _write_report(
        entries,
        out_dir,
        title,
        report_name,
        custom_parameters=custom_parameters,
    )
    return html_path, len(entries)


def build_multiverse_registration_qc_report(
    project_dir,
    atlas_path=DEFAULT_SIGMA_ATLAS,
    n_slices=10,
    custom_parameters=None,
):
    """Overlay every corrected temporal mean with the complete SIGMA atlas."""
    return _build_multiverse_report(
        project_dir,
        atlas_path,
        "registration",
        n_slices=n_slices,
        custom_parameters=custom_parameters,
    )


def build_multiverse_cc_qc_report(
    project_dir,
    atlas_path=DEFAULT_SIGMA_ATLAS,
    n_slices=10,
    custom_parameters=None,
):
    """Overlay every corrected temporal mean with SIGMA CC labels 891/892."""
    return _build_multiverse_report(
        project_dir,
        atlas_path,
        "cc",
        n_slices=n_slices,
        custom_parameters=custom_parameters,
    )


def _build_argument_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Create Multiverse registration and corpus-callosum QC reports by "
            "overlaying corrected temporal means with the SIGMA atlas."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="project_dir",
        required=True,
        type=Path,
        help="Processed project directory containing sub-*/ses-* folders.",
    )
    parser.add_argument(
        "--atlas",
        type=Path,
        default=DEFAULT_SIGMA_ATLAS,
        help=f"SIGMA atlas path (default: {DEFAULT_SIGMA_ATLAS}).",
    )
    parser.add_argument(
        "--report",
        choices=("all", "registration", "cc"),
        default="all",
        help="Report to create (default: all).",
    )
    parser.add_argument(
        "--n-slices",
        type=_positive_int,
        default=10,
        help="Number of slices per orientation (default: 10).",
    )
    parser.add_argument(
        "--custom-parameter",
        dest="custom_parameters",
        action="append",
        type=_custom_parameter,
        default=[],
        metavar="NAME=VALUE",
        help="Parameter to list in the HTML report. May be repeated.",
    )
    return parser


def main(argv=None):
    parser = _build_argument_parser()
    args = parser.parse_args(argv)
    project_dir = args.project_dir.expanduser()
    atlas_path = args.atlas.expanduser()

    if not project_dir.is_dir():
        parser.error(
            f"project directory does not exist or is not a directory: {project_dir}"
        )
    if not atlas_path.is_file():
        parser.error(f"SIGMA atlas does not exist: {atlas_path}")

    report_builders = []
    if args.report in ("all", "registration"):
        report_builders.append(
            ("Registration", build_multiverse_registration_qc_report)
        )
    if args.report in ("all", "cc"):
        report_builders.append(("CC", build_multiverse_cc_qc_report))

    for report_label, report_builder in report_builders:
        html_path, count = report_builder(
            project_dir,
            atlas_path=atlas_path,
            n_slices=args.n_slices,
            custom_parameters=args.custom_parameters or None,
        )
        if html_path:
            print(f"{report_label} report written to {html_path} ({count} image(s))")
        else:
            print(f"{report_label} report skipped: no matching files found.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
