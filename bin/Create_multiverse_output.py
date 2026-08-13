import os
import glob
import argparse
import subprocess

import nibabel as nib
import numpy as np
from tqdm import tqdm

from common.artifact_manifest import start_output_tracking


COMPOSITE_SUFFIX = "Matrixcomp_rsfMRI.nii.gz"
MAX_AFFINE_DISPLACEMENT_VOXELS = 0.1
SIGMA_ATLAS_FILENAME = "SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz"


def find_single_file(search_pattern, description):
    """
    Find exactly one file matching a glob pattern.

    Raises
    ------
    FileNotFoundError
        If no matching file is found.
    RuntimeError
        If more than one matching file is found.
    """
    matches = sorted(glob.glob(search_pattern, recursive=True))

    if not matches:
        raise FileNotFoundError(
            f"Could not find {description} using pattern:\n"
            f"{search_pattern}"
        )

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {description}, found {len(matches)}:\n"
            + "\n".join(matches)
        )

    return matches[0]


def find_registration_reference(composite_transform):
    """Return the fMRI BET image used to create a composite transform."""
    composite_name = os.path.basename(composite_transform)
    if not composite_name.endswith(COMPOSITE_SUFFIX):
        raise ValueError(
            f"Unexpected composite transformation filename: {composite_name}"
        )

    reference_stem = composite_name[: -len(COMPOSITE_SUFFIX)]
    reference_pattern = os.path.join(
        os.path.dirname(composite_transform),
        reference_stem + ".nii*",
    )
    return find_single_file(
        reference_pattern,
        "fMRI registration reference",
    )


def validate_spatial_geometry(input_file, registration_reference):
    """Ensure a 4D fMRI file uses the grid of the registration reference.

    ANTs can orthogonalise a NIfTI affine when qform and sform differ slightly,
    while the FSL processing path retains the original sform. Such differences
    do not imply a different voxel grid. Accept them when the images have the
    same shape and orientation and the maximum corner displacement remains
    below one tenth of the smallest spatial voxel size.
    """
    input_img = nib.load(input_file)
    reference_img = nib.load(registration_reference)

    if input_img.ndim != 4:
        raise ValueError(
            f"Expected a 4D fMRI file, but received {input_img.ndim} dimensions "
            f"for {input_file}."
        )

    input_shape = tuple(input_img.shape[:3])
    reference_shape = tuple(reference_img.shape[:3])
    if input_shape != reference_shape:
        raise ValueError(
            "The fMRI file and registration reference have different spatial "
            f"shapes: {input_shape} vs {reference_shape}.\n"
            f"fMRI: {input_file}\nReference: {registration_reference}"
        )

    input_orientation = nib.aff2axcodes(input_img.affine)
    reference_orientation = nib.aff2axcodes(reference_img.affine)
    if input_orientation != reference_orientation:
        raise ValueError(
            "The fMRI file and registration reference have different spatial "
            f"orientations: {input_orientation} vs {reference_orientation}.\n"
            f"fMRI: {input_file}\nReference: {registration_reference}"
        )

    spatial_corners = np.array(
        [
            [i, j, k, 1.0]
            for i in (0, input_shape[0] - 1)
            for j in (0, input_shape[1] - 1)
            for k in (0, input_shape[2] - 1)
        ]
    ).T
    world_difference = (
        (input_img.affine - reference_img.affine) @ spatial_corners
    )[:3]
    max_displacement_mm = float(
        np.max(np.linalg.norm(world_difference, axis=0))
    )
    min_voxel_size_mm = float(
        min(
            np.min(nib.affines.voxel_sizes(input_img.affine)),
            np.min(nib.affines.voxel_sizes(reference_img.affine)),
        )
    )
    max_allowed_mm = MAX_AFFINE_DISPLACEMENT_VOXELS * min_voxel_size_mm

    if max_displacement_mm > max_allowed_mm:
        raise ValueError(
            "The fMRI file and registration reference have materially different "
            f"spatial affines (maximum corner displacement "
            f"{max_displacement_mm:.6f} mm; allowed {max_allowed_mm:.6f} mm).\n"
            f"fMRI: {input_file}\nReference: {registration_reference}"
        )

    if not np.allclose(input_img.affine, reference_img.affine, rtol=0, atol=1e-4):
        print(
            "Warning: Accepting a minor fMRI/reference affine difference "
            f"(maximum corner displacement {max_displacement_mm:.6f} mm, "
            f"{max_displacement_mm / min_voxel_size_mm:.4f} voxel)."
        )


def validate_sigma_template_geometry(sigma_template_address):
    """Ensure the export grid matches the SIGMA atlas used in registration."""
    sigma_atlas_address = os.path.join(
        os.path.dirname(sigma_template_address),
        SIGMA_ATLAS_FILENAME,
    )
    if not os.path.exists(sigma_atlas_address):
        print(
            "Warning: Could not validate the SIGMA template geometry because "
            f"the matching atlas was not found: {sigma_atlas_address}"
        )
        return

    template_img = nib.load(sigma_template_address)
    atlas_img = nib.load(sigma_atlas_address)
    template_shape = tuple(template_img.shape[:3])
    atlas_shape = tuple(atlas_img.shape[:3])

    if template_shape != atlas_shape or not np.allclose(
        template_img.affine,
        atlas_img.affine,
        rtol=0,
        atol=1e-4,
    ):
        expected_template = os.path.join(
            os.path.dirname(sigma_template_address),
            "SIGMA_InVivo_Brain_Template_Masked.nii.gz",
        )
        raise ValueError(
            "The SIGMA template does not use the geometry of the SIGMA atlas "
            "used for registration. Use the matching compressed template:\n"
            f"{expected_template}\n"
            f"Template: {sigma_template_address}\n"
            f"Atlas: {sigma_atlas_address}"
        )


def apply_inverse_composite_transformation(
    files_list,
    func_folder,
    sigma_template_address,
):
    """
    Transform fMRI files into SIGMA space by inverting the exact composite
    transformation used to resample the SIGMA atlas into fMRI space.

    The input fMRI files are passed directly to reg_resample.
    No flipping or additional reorientation is performed.
    """
    transformed_files = []

    search_root_func = os.path.dirname(func_folder)

    validate_sigma_template_geometry(sigma_template_address)

    composite_transform = find_single_file(
        os.path.join(
            search_root_func,
            f"*{COMPOSITE_SUFFIX}",
        ),
        "fMRI-to-SIGMA composite transformation",
    )

    registration_reference = find_registration_reference(composite_transform)
    inverse_composite = os.path.join(
        func_folder,
        os.path.basename(composite_transform).replace(
            ".nii.gz",
            "_inv.nii.gz",
        ),
    )

    # The composite transform was estimated on the 3D fMRI BET reference.
    # It is valid for the 4D data only if their spatial grids are identical.
    for input_file in files_list:
        validate_spatial_geometry(input_file, registration_reference)

    # Matrixcomp_rsfMRI maps fMRI reference coordinates into SIGMA floating
    # coordinates. Invert the complete (affine + non-linear) transform and
    # discretise the inverse deformation on the SIGMA target grid.
    subprocess.run(
        [
            "reg_transform",
            "-invNrr",
            composite_transform,
            sigma_template_address,
            inverse_composite,
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
                inverse_composite,
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
    func_root = os.path.join(subject_path, "func")
    tracker = start_output_tracking(func_root, "func", "processing")

    try:
        registered_files = apply_inverse_composite_transformation(
            files_list=[epi_file],
            func_folder=func_folder,
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
    finally:
        # Finalize per subject because one batch invocation processes several
        # independent func folders.
        tracker.finalize()


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
