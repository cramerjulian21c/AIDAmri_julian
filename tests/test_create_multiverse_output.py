"""Tests for transforming processed 4D fMRI data into SIGMA space."""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
MODULE_PATH = BIN_DIR / "Create_multiverse_output.py"
sys.path.insert(0, str(BIN_DIR))


def load_module():
    spec = importlib.util.spec_from_file_location(
        "create_multiverse_output_transform_test",
        MODULE_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def save_nifti(path, shape, affine=None):
    if affine is None:
        affine = np.eye(4)
    nib.save(nib.Nifti1Image(np.zeros(shape, dtype=np.float32), affine), path)


class MultiverseTransformTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_inverts_and_applies_complete_composite_transform(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)

            prefix = "sub-01_task-rest_bold_EPIMPSmoothBet"
            reference = func_root / f"{prefix}.nii.gz"
            composite = func_root / f"{prefix}Matrixcomp_rsfMRI.nii.gz"
            template = Path(temp_dir) / "SIGMA_template.nii.gz"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            save_nifti(template, (7, 8, 9))
            save_nifti(atlas, (7, 8, 9))
            save_nifti(input_file, (4, 5, 6, 3))
            composite.write_bytes(b"composite")

            with mock.patch.object(self.module.subprocess, "run") as run:
                outputs = self.module.apply_inverse_composite_transformation(
                    files_list=[str(input_file)],
                    func_folder=str(func_folder),
                    sigma_template_address=str(template),
                )

            inverse = func_folder / f"{prefix}Matrixcomp_rsfMRI_inv.nii.gz"
            output = func_folder / (
                "sub-01_task-rest_bold_EPI_mcf_f_"
                "registered_on_SIGMA_template.nii.gz"
            )
            self.assertEqual(outputs, [str(output)])
            self.assertEqual(
                run.call_args_list,
                [
                    mock.call(
                        [
                            "reg_transform",
                            "-invNrr",
                            str(composite),
                            str(template),
                            str(inverse),
                        ],
                        check=True,
                    ),
                    mock.call(
                        [
                            "reg_resample",
                            "-ref",
                            str(template),
                            "-flo",
                            str(input_file),
                            "-trans",
                            str(inverse),
                            "-res",
                            str(output),
                            "-inter",
                            "3",
                        ],
                        check=True,
                    ),
                ],
            )

    def test_rejects_fmri_on_a_different_spatial_grid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)

            prefix = "sub-01_task-rest_bold_EPIMPSmoothBet"
            reference = func_root / f"{prefix}.nii.gz"
            composite = func_root / f"{prefix}Matrixcomp_rsfMRI.nii.gz"
            template = Path(temp_dir) / "SIGMA_template.nii.gz"
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            save_nifti(template, (7, 8, 9))
            save_nifti(input_file, (4, 5, 7, 3))
            composite.write_bytes(b"composite")

            with (
                mock.patch.object(self.module.subprocess, "run") as run,
                self.assertRaisesRegex(ValueError, "different spatial shapes"),
            ):
                self.module.apply_inverse_composite_transformation(
                    files_list=[str(input_file)],
                    func_folder=str(func_folder),
                    sigma_template_address=str(template),
                )

            run.assert_not_called()

    def test_accepts_minor_ants_affine_normalisation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "fmri.nii.gz"
            reference = Path(temp_dir) / "reference.nii.gz"
            affine = np.diag([-0.4, -0.4, -1.1, 1.0])
            reference_affine = affine.copy()
            reference_affine[1, 2] = 0.0006

            save_nifti(input_file, (64, 64, 18, 3), affine)
            save_nifti(reference, (64, 64, 18), reference_affine)

            with mock.patch("builtins.print") as print_mock:
                self.module.validate_spatial_geometry(
                    str(input_file),
                    str(reference),
                )

            self.assertIn(
                "Accepting a minor fMRI/reference affine difference",
                print_mock.call_args.args[0],
            )

    def test_rejects_material_affine_difference(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir) / "fmri.nii.gz"
            reference = Path(temp_dir) / "reference.nii.gz"
            affine = np.diag([-0.4, -0.4, -1.1, 1.0])
            shifted_affine = affine.copy()
            shifted_affine[0, 3] = 0.2

            save_nifti(input_file, (64, 64, 18, 3), affine)
            save_nifti(reference, (64, 64, 18), shifted_affine)

            with self.assertRaisesRegex(
                ValueError,
                "materially different spatial affines",
            ):
                self.module.validate_spatial_geometry(
                    str(input_file),
                    str(reference),
                )

    def test_rejects_sigma_template_with_wrong_atlas_affine(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            template = Path(temp_dir) / "SIGMA_InVivo_Brain_Template_Masked.nii"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            template_affine = np.eye(4)
            atlas_affine = np.eye(4)
            atlas_affine[1, 3] = 11.045

            save_nifti(template, (128, 128, 218), template_affine)
            save_nifti(atlas, (128, 128, 218), atlas_affine)

            with self.assertRaisesRegex(
                ValueError,
                "does not use the geometry of the SIGMA atlas",
            ):
                self.module.validate_sigma_template_geometry(str(template))

    def test_accepts_sigma_template_matching_atlas_geometry(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            template = Path(temp_dir) / "SIGMA_InVivo_Brain_Template_Masked.nii.gz"
            atlas = Path(temp_dir) / self.module.SIGMA_ATLAS_FILENAME
            affine = np.diag([-0.15, -0.15, -0.15, 1.0])

            save_nifti(template, (128, 128, 218), affine)
            save_nifti(atlas, (128, 128, 218), affine)

            self.module.validate_sigma_template_geometry(str(template))

    def test_rejects_ambiguous_composite_transforms(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            func_root = Path(temp_dir) / "func"
            func_folder = func_root / "rs-fMRI_niiData"
            func_folder.mkdir(parents=True)
            (func_root / "firstMatrixcomp_rsfMRI.nii.gz").write_bytes(b"first")
            (func_root / "secondMatrixcomp_rsfMRI.nii.gz").write_bytes(b"second")

            with self.assertRaisesRegex(RuntimeError, "found 2"):
                self.module.apply_inverse_composite_transformation(
                    files_list=[],
                    func_folder=str(func_folder),
                    sigma_template_address="template.nii.gz",
                )


if __name__ == "__main__":
    unittest.main()
