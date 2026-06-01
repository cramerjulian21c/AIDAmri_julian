import os
import subprocess

import nibabel as nib
import nipype.interfaces.fsl as fsl
import numpy as np


FSL_BET_WORLD_SWAPS = [(1, 2)]


def set_default_xyzt_units_if_unknown(target):
    """
       Set NIfTI space/time units to mm/sec only if they are unknown.

       Accepts either:
       - a file path to a NIfTI file, or
       - a nibabel header object
       """
    if isinstance(target, str):
        img = nib.load(target)
        hdr = img.header

        space_unit, time_unit = hdr.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        hdr.set_xyzt_units(space_unit, time_unit)
        nib.save(img, target)
        return target

    if hasattr(target, "get_xyzt_units") and hasattr(target, "set_xyzt_units"):
        space_unit, time_unit = target.get_xyzt_units()
        if not space_unit or space_unit.lower() == "unknown":
            space_unit = "mm"
        if not time_unit or time_unit.lower() == "unknown":
            time_unit = "sec"

        target.set_xyzt_units(space_unit, time_unit)
        return target

    raise TypeError("Expected a NIfTI file path or nibabel header object")



def estimate_center_intensity_based(nifti, percentile=60):
    """
        Estimate BET center (-c) using intensity-weighted center-of-gravity (fslstats -C),
        but excluding low-intensity voxels using a data-adaptive threshold (-l = P{percentile}).
        """
    # 1) Get intensity percentile
    p = subprocess.check_output(
        ["fslstats", nifti, "-P", str(percentile)]
    ).decode().strip()

    # 2) Compute center-of-gravity using only voxels > P{percentile}
    center = subprocess.check_output(
        ["fslstats", nifti, "-l", p, "-C"]
    ).decode().strip().split()

    cx, cy, cz = [float(v) for v in center]
    return [cx, cy, cz], float(p)


def apply_world_ops(mat, swaps=()):
    out = mat.copy()
    # swaps first
    for a, b in swaps:
        out[[a, b], :] = out[[b, a], :]

    return out


def save_header_only_reoriented_copy(src_path, dst_path, swaps=()):
    img = nib.load(src_path)
    data = img.get_fdata(dtype=np.float32)
    aff = img.affine.copy()

    # header-only world-axis operation
    aff[:3, :] = apply_world_ops(aff[:3, :], swaps=swaps)

    hdr = img.header.copy()
    hdr["pixdim"][0] = 1
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    out = nib.Nifti1Image(np.ascontiguousarray(data, np.float32), aff, header=hdr)
    out.set_qform(aff, code=1)
    out.set_sform(aff, code=1)
    nib.save(out, dst_path)

    return dst_path


def _format_result(output_file, mask_file, return_mask):
    '''
    Gives back imageBet.nii.gz for Preprocessing und mask for process_fMRI.py and regress.py
    '''
    if return_mask:
        return output_file, mask_file
    return output_file


def apply_bet(
    input_file,
    frac,
    radius,
    horizontal_gradient=0.0,
    use_bet4animal=False,
    species="mouse",
    center=None,
    return_mask=False,
    output_path=None,
    verbose=False,
):
    """Apply the shared AIDAmri BET implementation."""
    if output_path is None:
        output_path = os.path.dirname(input_file)

    if use_bet4animal:
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2
        species_id = 6 if species == "mouse" else 5
        output_file = os.path.join(
            output_path,
            os.path.basename(input_file).split(".")[0] + "AnimalBet.nii.gz",
        )

        tmp_bet = os.path.join(os.path.dirname(input_file), "bet4animal_tmp_out.nii.gz")
        tmp_mask = tmp_bet.replace(".nii.gz", "_mask.nii.gz")
        mask_file = output_file.replace(".nii.gz", "_mask.nii.gz")
        tmp_std = os.path.join(os.path.dirname(input_file), "bet4animal_reorient2std.nii.gz")

        subprocess.run(["fslreorient2std", input_file, tmp_std], check=True)

        bet_in = tmp_std
        if center is None:
            center, _ = estimate_center_intensity_based(bet_in)
        cx, cy, cz = center

        cmd = [
            "/aida/bin/bet4animal",
            bet_in,
            tmp_bet,
            "-m",
            "-w", str(w_value),
            "-z", str(species_id),
            "-c", str(cx), str(cy), str(cz),
        ]
        subprocess.run(cmd, check=True)

        input_img = nib.load(input_file)
        target_axcodes = nib.aff2axcodes(input_img.affine)

        img = nib.load(tmp_bet)
        data = img.get_fdata(dtype=np.float32)
        aff = img.affine

        ornt_cur = nib.orientations.io_orientation(aff)
        ornt_tgt = nib.orientations.axcodes2ornt(target_axcodes)
        transform = nib.orientations.ornt_transform(ornt_cur, ornt_tgt)

        data_lip = nib.orientations.apply_orientation(data, transform)
        aff_lip = aff @ nib.orientations.inv_ornt_aff(transform, img.shape)

        hdr_final = img.header.copy()
        hdr_final.set_data_dtype(np.float32)
        hdr_final["pixdim"][0] = 1
        hdr_final["pixdim"][4:8] = 1
        set_default_xyzt_units_if_unknown(hdr_final)

        img_final = nib.Nifti1Image(
            np.ascontiguousarray(data_lip, dtype=np.float32),
            aff_lip,
            header=hdr_final,
        )
        img_final.set_qform(aff_lip, code=1)
        img_final.set_sform(aff_lip, code=1)
        nib.save(img_final, output_file)

        if os.path.exists(tmp_mask):
            m_img = nib.load(tmp_mask)
            m_data = m_img.get_fdata(dtype=np.float32)
            m_aff = m_img.affine

            m_ornt_cur = nib.orientations.io_orientation(m_aff)
            m_transform = nib.orientations.ornt_transform(m_ornt_cur, ornt_tgt)

            m_data_lip = nib.orientations.apply_orientation(m_data, m_transform)
            m_aff_lip = m_aff @ nib.orientations.inv_ornt_aff(m_transform, m_img.shape)

            m_bin = (m_data_lip > 0.5).astype(np.uint8)

            m_hdr_lip = m_img.header.copy()
            m_hdr_lip.set_data_dtype(np.uint8)
            m_hdr_lip["pixdim"][0] = 1
            m_hdr_lip["pixdim"][4:8] = 1
            set_default_xyzt_units_if_unknown(m_hdr_lip)

            m_out = nib.Nifti1Image(
                np.ascontiguousarray(m_bin, dtype=np.uint8),
                m_aff_lip,
                header=m_hdr_lip,
            )
            m_out.set_qform(m_aff_lip, code=1)
            m_out.set_sform(m_aff_lip, code=1)

            nib.save(m_out, mask_file)
        else:
            print("Warning: BET mask not found:", tmp_mask)

        for tmp_file in [tmp_std, tmp_bet, tmp_mask]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

        print(f"Brain extraction completed, output saved to {output_file}")
        return _format_result(output_file, mask_file, return_mask)

    data = nib.load(input_file)
    img_temp = data.get_fdata()
    if verbose:
        print("Image dimensions before scaling:", data.header.get_zooms())

    scale = np.eye(4)
    scale[0, 0] = 10
    scale[1, 1] = 10
    scale[2, 2] = 10
    if verbose:
        print("Image dimensions after scaling:", (data.affine * scale)[:3, :3])

    scaled_affine = data.affine @ scale
    scaled_img = nib.Nifti1Image(img_temp, scaled_affine)
    hdr_in = scaled_img.header
    set_default_xyzt_units_if_unknown(hdr_in)

    fsl_path = os.path.join(os.path.dirname(input_file), "fslScaleTemp.nii.gz")
    nib.save(scaled_img, fsl_path)
    if verbose:
        print("Saved scaled image to:", fsl_path)
        print("Image dimensions:", scaled_img.header.get_zooms())

    bet_input_tmp = os.path.join(os.path.dirname(input_file), "fslScaleTemp_LPIhdr.nii.gz")
    save_header_only_reoriented_copy(
        fsl_path,
        bet_input_tmp,
        swaps=FSL_BET_WORLD_SWAPS,
    )

    output_file = os.path.join(
        output_path,
        os.path.basename(input_file).split(".")[0] + "Bet.nii.gz",
    )
    bet_output_tmp = os.path.join(
        os.path.dirname(input_file),
        os.path.basename(input_file).split(".")[0] + "Bet_LPIhdr_tmp.nii.gz",
    )

    if center is not None:
        center_int = [int(round(c)) for c in center]
        print(f"Using user-defined center (rounded): {center_int}")
        my_bet = fsl.BET(
            in_file=bet_input_tmp,
            out_file=bet_output_tmp,
            frac=frac,
            radius=radius,
            vertical_gradient=horizontal_gradient,
            center=center_int,
            mask=True,
        )
    else:
        print("Using robust center estimation (-R)")
        my_bet = fsl.BET(
            in_file=bet_input_tmp,
            out_file=bet_output_tmp,
            frac=frac,
            radius=radius,
            vertical_gradient=horizontal_gradient,
            robust=True,
            mask=True,
        )
    if verbose:
        print(my_bet.cmdline)
    my_bet.run()

    save_header_only_reoriented_copy(
        bet_output_tmp,
        output_file,
        swaps=FSL_BET_WORLD_SWAPS,
    )

    mask_tmp_file = bet_output_tmp.replace(".nii.gz", "_mask.nii.gz")
    mask_file = output_file.replace(".nii.gz", "_mask.nii.gz")

    if os.path.exists(mask_tmp_file):
        save_header_only_reoriented_copy(
            mask_tmp_file,
            mask_file,
            swaps=FSL_BET_WORLD_SWAPS,
        )

    data_out = nib.load(output_file)
    img_out = data_out.get_fdata(dtype=np.float32)

    inv_scale = np.eye(4)
    inv_scale[0, 0] = 0.1
    inv_scale[1, 1] = 0.1
    inv_scale[2, 2] = 0.1

    unscaled_affine = data_out.affine @ inv_scale
    unscaled_img = nib.Nifti1Image(img_out, unscaled_affine)
    unscaled_img.set_qform(unscaled_affine, code=1)
    unscaled_img.set_sform(unscaled_affine, code=1)
    hdr_out = unscaled_img.header
    set_default_xyzt_units_if_unknown(hdr_out)
    if verbose:
        print("Image dimensions after unscaling:", unscaled_img.header.get_zooms())
    nib.save(unscaled_img, output_file)

    if os.path.exists(mask_file):
        mask_data = nib.load(mask_file)
        bet_ref = nib.load(output_file)
        mask_img = (mask_data.get_fdata() > 0.5).astype(np.uint8)

        final_mask = nib.Nifti1Image(mask_img, bet_ref.affine)
        final_mask.set_qform(bet_ref.affine, code=1)
        final_mask.set_sform(bet_ref.affine, code=1)

        hdr_mask = final_mask.header
        hdr_mask.set_data_dtype(np.uint8)
        set_default_xyzt_units_if_unknown(hdr_mask)

        nib.save(final_mask, mask_file)

    for tmp_file in [fsl_path, bet_input_tmp, bet_output_tmp, mask_tmp_file]:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

    print(f"Brain extraction completed, output saved to {output_file}")
    return _format_result(output_file, mask_file, return_mask)


def skip_bet_function(input_file, return_mask=False):
    print("Skipping BET")
    print("Creating BET-compatible outputs for pipeline compatibility")

    output_bet = os.path.join(
        os.path.dirname(input_file),
        os.path.basename(input_file).split(".")[0] + "Bet.nii.gz",
    )

    src = nib.load(input_file)
    data = src.get_fdata(dtype=np.float32)
    aff = src.affine.copy()

    hdr = src.header.copy()
    hdr.set_data_dtype(np.float32)
    hdr["pixdim"][0] = 1
    hdr["pixdim"][4:8] = 1
    set_default_xyzt_units_if_unknown(hdr)

    final_img = nib.Nifti1Image(
        np.ascontiguousarray(data, dtype=np.float32),
        aff,
        header=hdr,
    )
    final_img.set_qform(aff, code=1)
    final_img.set_sform(aff, code=1)
    nib.save(final_img, output_bet)

    print(f"BET skipped -> created compatibility image: {output_bet}")

    mask_file = output_bet.replace(".nii.gz", "_mask.nii.gz")
    mask = (data > 0).astype(np.uint8)

    mask_hdr = hdr.copy()
    mask_hdr.set_data_dtype(np.uint8)
    mask_hdr["pixdim"][0] = 1
    mask_hdr["pixdim"][4:8] = 1
    set_default_xyzt_units_if_unknown(mask_hdr)

    mask_img = nib.Nifti1Image(
        np.ascontiguousarray(mask, dtype=np.uint8),
        aff,
        header=mask_hdr,
    )
    mask_img.set_qform(aff, code=1)
    mask_img.set_sform(aff, code=1)
    nib.save(mask_img, mask_file)

    print(f"BET mask created at {mask_file}")
    return _format_result(output_bet, mask_file, return_mask)
