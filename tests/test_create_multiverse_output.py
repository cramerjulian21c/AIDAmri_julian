"""Tests for transforming processed 4D fMRI data into SIGMA space."""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
MODULE_PATH = BIN_DIR / "Create_multiverse_output.py"


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
            input_file = func_folder / "sub-01_task-rest_bold_EPI_mcf_f.nii.gz"

            save_nifti(reference, (4, 5, 6))
            save_nifti(template, (7, 8, 9))
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
