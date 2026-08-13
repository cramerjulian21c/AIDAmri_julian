"""Regression tests for AIDAmri manifests and the manifest-based reset.

This module is used for development and quality assurance; it is not required
to process MRI data. It verifies that newly created and modified files are
recorded in the modality-specific manifest and that ``Reset_proc_folder.py``
deletes only registered outputs from stages after the selected target phase.

The tests also cover manifest naming, T2map support, preservation of unmanaged
files, and missing manifests inside or outside the selected modality. Every
test uses an automatically created temporary directory; real project data is
never read or modified.

Run from the repository root::

    python -m unittest discover -s tests -v
"""

import contextlib
import importlib.util
import io
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from common.artifact_manifest import OutputTracker, load_manifest, manifest_filename
from helper_tools.Reset_proc_folder import (
    apply_reset_plan,
    build_reset_plan,
    main as reset_main,
    print_issue_summary,
    print_reset_plan,
    validate_project_manifests,
)


class ArtifactManifestTests(unittest.TestCase):
    def test_multiverse_outputs_are_tracked_in_separate_stage(self):
        module_path = BIN_DIR / "Create_multiverse_output.py"
        spec = importlib.util.spec_from_file_location(
            "create_multiverse_output_test",
            module_path,
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory() as temp_dir:
            subject = Path(temp_dir) / "sub-01" / "ses-1"
            func_root = subject / "func"
            func_output = func_root / "rs-fMRI_niiData"
            anat = subject / "anat"
            func_output.mkdir(parents=True)
            anat.mkdir()
            epi = func_output / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"
            epi.write_text("input", encoding="utf-8")

            registered = func_output / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template.nii.gz"
            )
            temporal_mean = func_output / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template_temporal_mean.nii.gz"
            )

            def fake_transform(**_kwargs):
                registered.write_text("registered", encoding="utf-8")
                return [str(registered)]

            def fake_temporal_mean(_input_path, output_path):
                Path(output_path).write_text("mean", encoding="utf-8")

            with (
                mock.patch.object(
                    module,
                    "apply_inverse_composite_transformation",
                    side_effect=fake_transform,
                ),
                mock.patch.object(
                    module,
                    "compute_temporal_mean",
                    side_effect=fake_temporal_mean,
                ),
            ):
                module.process_subject(str(subject), "template.nii.gz")

            manifest = load_manifest(func_root, "func")
            self.assertEqual(manifest["mode"], "func")
            self.assertEqual(
                manifest["stages"]["multiverse_output"]["created_files"],
                [
                    registered.relative_to(func_root).as_posix(),
                    temporal_mean.relative_to(func_root).as_posix(),
                ],
            )

    def test_reset_to_processing_deletes_only_multiverse_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "func"
            (folder / "brkraw").mkdir(parents=True)

            processing = OutputTracker.start(folder, "func", "processing")
            processing_output = folder / "processed.nii.gz"
            processing_output.write_text("processed", encoding="utf-8")
            processing.finalize()

            multiverse = OutputTracker.start(
                folder,
                "func",
                "multiverse_output",
            )
            multiverse_output = folder / "multiverse.nii.gz"
            multiverse_output.write_text("multiverse", encoding="utf-8")
            multiverse.finalize()

            plans = build_reset_plan(project, "func", "processing")
            self.assertEqual(plans[0]["selected_stages"], ["multiverse_output"])

            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertTrue(processing_output.exists())
            self.assertFalse(multiverse_output.exists())
            manifest = load_manifest(folder, "func")
            self.assertIn("processing", manifest["stages"])
            self.assertNotIn("multiverse_output", manifest["stages"])

    def test_tracker_records_created_and_modified_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "sub-test" / "ses-test" / "anat"
            folder.mkdir(parents=True)
            existing = folder / "input.nii.gz"
            existing.write_text("before", encoding="utf-8")

            tracker = OutputTracker.start(folder, "anat", "preprocessing")
            existing.write_text("after", encoding="utf-8")
            (folder / "output.nii.gz").write_text("new", encoding="utf-8")
            tracker.finalize()

            self.assertEqual(
                manifest_filename(folder, "anat"),
                ".sub-test_ses-test_anat_aidamri_manifest.json",
            )
            manifest_file = folder / manifest_filename(folder, "anat")
            self.assertEqual(
                stat.S_IMODE(manifest_file.stat().st_mode),
                0o644,
            )
            stage = load_manifest(folder, "anat")["stages"]["preprocessing"]
            self.assertEqual(stage["created_files"], ["output.nii.gz"])
            self.assertEqual(stage["modified_files"], ["input.nii.gz"])

    def test_reset_deletes_only_manifested_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "ses-01" / "anat"
            (folder / "brkraw").mkdir(parents=True)
            raw = folder / "input_T2w.nii.gz"
            raw.write_text("raw", encoding="utf-8")

            preprocessing = OutputTracker.start(folder, "anat", "preprocessing")
            preprocessing_output = folder / "preprocessed.nii.gz"
            preprocessing_output.write_text("pre", encoding="utf-8")
            preprocessing.finalize()

            registration = OutputTracker.start(folder, "anat", "registration")
            registration_dir = folder / "registration_outputs"
            registration_dir.mkdir()
            registration_output = registration_dir / "registered.nii.gz"
            registration_output.write_text("reg", encoding="utf-8")
            raw.write_text("changed by registration", encoding="utf-8")
            registration.finalize()

            unknown = registration_dir / "notes.txt"
            unknown.write_text("user data", encoding="utf-8")

            plans = build_reset_plan(project, "anat", "preprocessing")
            self.assertIn(unknown, plans[0]["unknown_files"])
            report = io.StringIO()
            with contextlib.redirect_stdout(report):
                print_reset_plan(plans, "anat", "preprocessing")
            self.assertIn(
                f"Not managed by the manifest; preserved: {unknown}",
                report.getvalue(),
            )
            summary = io.StringIO()
            with contextlib.redirect_stdout(summary):
                print_issue_summary(plans)
            summary_text = summary.getvalue()
            self.assertIn("sub-01 | ses-01 | anat", summary_text)
            self.assertIn(
                f"[INFO] Not managed by the manifest; review manually: {unknown}",
                summary_text,
            )
            self.assertIn(
                f"[WARNING] Modified by registration; original content cannot be "
                f"restored: {raw}",
                summary_text,
            )
            with contextlib.redirect_stdout(io.StringIO()):
                apply_reset_plan(plans)

            self.assertTrue(raw.exists())
            self.assertTrue(preprocessing_output.exists())
            self.assertFalse(registration_output.exists())
            self.assertTrue(unknown.exists())
            self.assertTrue(registration_dir.exists())

    def test_folder_without_manifest_is_skipped(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            folder = project / "sub-01" / "anat"
            (folder / "brkraw").mkdir(parents=True)
            unknown = folder / "unknown.nii.gz"
            unknown.write_text("keep", encoding="utf-8")

            plans = build_reset_plan(project, "anat", "base")

            self.assertTrue(plans[0]["skip"])
            self.assertTrue(unknown.exists())

    def test_t2map_outputs_can_be_tracked(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir) / "sub-test" / "ses-test" / "t2map"
            folder.mkdir(parents=True)
            tracker = OutputTracker.start(folder, "t2map", "processing")
            (folder / "t2values.csv").write_text("value", encoding="utf-8")
            tracker.finalize()

            manifest = load_manifest(folder, "t2map")
            self.assertEqual(manifest["mode"], "t2map")
            self.assertEqual(
                manifest["stages"]["processing"]["created_files"],
                ["t2values.csv"],
            )

    def test_reset_ignores_missing_manifest_in_other_modality(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            anat = project / "sub-01" / "ses-01" / "anat"
            dwi = project / "sub-01" / "ses-01" / "dwi"
            (anat / "brkraw").mkdir(parents=True)
            (dwi / "brkraw").mkdir(parents=True)

            tracker = OutputTracker.start(anat, "anat", "preprocessing")
            anat_output = anat / "preprocessed.nii.gz"
            anat_output.write_text("pre", encoding="utf-8")
            tracker.finalize()

            self.assertEqual(validate_project_manifests(project, "anat"), [])
            dwi_errors = validate_project_manifests(project, "dwi")
            self.assertEqual(len(dwi_errors), 1)
            self.assertIn(
                "/dwi/.sub-01_ses-01_dwi_aidamri_manifest.json",
                dwi_errors[0],
            )

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "anat",
                        "--phase",
                        "base",
                        "--yes",
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertFalse(anat_output.exists())

    def test_reset_aborts_when_selected_modality_manifest_is_missing(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir)
            dwi = project / "sub-01" / "ses-01" / "dwi"
            (dwi / "brkraw").mkdir(parents=True)
            unknown = dwi / "unknown.nii.gz"
            unknown.write_text("keep", encoding="utf-8")

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = reset_main(
                    [
                        "--input",
                        str(project),
                        "--mode",
                        "dwi",
                        "--phase",
                        "base",
                        "--yes",
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertTrue(unknown.exists())


if __name__ == "__main__":
    unittest.main()
