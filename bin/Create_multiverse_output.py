import os
import glob
import argparse
import subprocess

import nibabel as nib
import numpy as np
from tqdm import tqdm


def find_first_file(search_pattern, description):
    """
    Find the first file matching a glob pattern.

    Raises
    ------
    FileNotFoundError
        If no matching file is found.
    """
    matches = sorted(glob.glob(search_pattern, recursive=True))

    if not matches:
        raise FileNotFoundError(
            f"Could not find {description} using pattern:\n"
            f"{search_pattern}"
        )

    if len(matches) > 1:
        print(
            f"Warning: Found multiple files for {description}. "
            f"Using:\n{matches[0]}"
        )

    return matches[0]


def apply_affine_transformations(
    files_list,
    func_folder,
    anat_folder,
    sigma_template_address,
):
    """
    Transform fMRI files into SIGMA-template space using existing
    functional and anatomical affine transformation matrices.

    The input fMRI files are passed directly to reg_resample.
    No flipping or additional reorientation is performed.
    """
    transformed_files = []

    search_root_func = os.path.dirname(func_folder)

    func_trafo = find_first_file(
        os.path.join(
            search_root_func,
            "**",
            "*transMatrixAff.txt",
        ),
        "functional affine transformation",
    )

    anat_trafo_inv = find_first_file(
        os.path.join(
            anat_folder,
            "**",
            "*MatrixInv.txt",
        ),
        "inverse anatomical transformation",
    )

    func_trafo_inv = os.path.join(
        func_folder,
        "func_trafo_inv.txt",
    )

    merged_inverted = os.path.join(
        func_folder,
        "merged_inverted.txt",
    )

    # Invert the functional affine transformation.
    subprocess.run(
        [
            "reg_transform",
            "-invAff",
            func_trafo,
            func_trafo_inv,
        ],
        check=True,
    )

    # Combine the inverse anatomical and inverse functional transforms.
    subprocess.run(
        [
            "reg_transform",
            "-comp",
            anat_trafo_inv,
            func_trafo_inv,
            merged_inverted,
        ],
        check=True,
    )

    for input_file in tqdm(
        files_list,
        desc="Applying transformations",
        unit="file",
    ):
        file_dir = os.path.dirname(input_file)
        base_filename = os.path.basename(input_file)

        output_filename = base_filename.replace(
            ".nii.gz",
            "_registered_on_SIGMA_template.nii.gz",
        )

        output_file = os.path.join(
            file_dir,
            output_filename,
        )

        # Use the original fMRI file directly as the floating image.
        subprocess.run(
            [
                "reg_resample",
                "-ref",
                sigma_template_address,
                "-flo",
                input_file,
                "-trans",
                merged_inverted,
                "-res",
                output_file,
                "-inter",
                "3",
            ],
            check=True,
        )

        transformed_files.append(output_file)

    return transformed_files


def compute_temporal_mean(four_d_img_path, output_path):
    """
    Calculate the mean across the temporal dimension of a 4D NIfTI file.
    """
    try:
        img = nib.load(four_d_img_path)
        data = img.get_fdata()

        if data.ndim != 4:
            raise ValueError(
                f"Expected a 4D fMRI file, but received data with "
                f"{data.ndim} dimensions and shape {data.shape}."
            )

        mean_data = np.mean(data, axis=3)

        output_header = img.header.copy()
        output_header.set_data_shape(mean_data.shape)

        mean_img = nib.Nifti1Image(
            mean_data,
            img.affine,
            output_header,
        )

        nib.save(mean_img, output_path)
        print(f"Temporal mean saved at: {output_path}")

    except Exception as error:
        print(
            f"Error computing temporal mean for "
            f"{four_d_img_path}: {error}"
        )
        raise


def process_subject(subject_path, template_path):
    """
    Register one subject's fMRI data to the SIGMA template and
    calculate the temporal mean.
    """
    func_folder = os.path.join(
        subject_path,
        "func",
        "rs-fMRI_niiData",
    )

    anat_folder = os.path.join(
        subject_path,
        "anat",
    )

    epi_pattern = os.path.join(
        func_folder,
        "*_task-rest_bold_EPI_mcf_st_f.nii.gz",
    )
    epi_files = sorted(glob.glob(epi_pattern))

    if not epi_files:
        epi_pattern = os.path.join(
            func_folder,
            "*_task-rest_bold_EPI_mcf_f.nii.gz",
        )
        epi_files = sorted(glob.glob(epi_pattern))

    if not epi_files:
        print(
            f"No EPI files found in {func_folder} using pattern:\n"
            f"{epi_pattern}"
        )
        return

    if len(epi_files) > 1:
        print(
            f"Warning: Found multiple EPI files in {func_folder}. "
            f"Using:\n{epi_files[0]}"
        )

    epi_file = epi_files[0]

    registered_files = apply_affine_transformations(
        files_list=[epi_file],
        func_folder=func_folder,
        anat_folder=anat_folder,
        sigma_template_address=template_path,
    )

    if not registered_files:
        print(f"Registration failed for: {subject_path}")
        return

    registered_file = registered_files[0]

    output_mean_path = registered_file.replace(
        ".nii.gz",
        "_temporal_mean.nii.gz",
    )

    compute_temporal_mean(
        registered_file,
        output_mean_path,
    )


def main_batch(root_folder, template_path):
    """
    Process every sub-*/ses-1 directory below the given root folder.
    """
    subject_dirs = sorted(
        glob.glob(
            os.path.join(
                root_folder,
                "sub-*",
                "ses-1",
            )
        )
    )

    if not subject_dirs:
        print(
            f"No subject session directories found under: "
            f"{root_folder}"
        )
        return

    print(f"Found {len(subject_dirs)} subject session directories.")

    for subject_path in subject_dirs:
        print(f"\nProcessing subject: {subject_path}")

        try:
            process_subject(
                subject_path,
                template_path,
            )
        except subprocess.CalledProcessError as error:
            print(
                f"NiftyReg command failed for {subject_path} "
                f"with exit code {error.returncode}."
            )
        except Exception as error:
            print(f"Error processing {subject_path}: {error}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Register 4D fMRI data to a SIGMA template using existing "
            "NiftyReg transformations and compute a temporal mean image."
        )
    )

    parser.add_argument(
        "-r",
        "--root_folder",
        required=True,
        help="Root directory containing sub-*/ses-1 folders.",
    )

    parser.add_argument(
        "-t",
        "--template_path",
        required=True,
        help="Path to the SIGMA template in NIfTI format.",
    )

    args = parser.parse_args()

    main_batch(
        root_folder=os.path.abspath(args.root_folder),
        template_path=os.path.abspath(args.template_path),
    )
