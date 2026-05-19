"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""



import nipype.interfaces.fsl as fsl
import os,sys
import nibabel as nib
import numpy as np
import nipype.interfaces.ants as ants
import subprocess
import shutil

FATAL_LIP_HEADER_EXIT_CODE = 86


def create_brkraw_backup(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    shutil.copyfile(input_file, dst_path)

    data = nib.load(input_file)
    raw_img = data.dataobj.get_unscaled().astype(np.float32)

    hdr = data.header.copy()
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    raw_nii = nib.Nifti1Image(raw_img, data.affine, header=hdr)
    raw_nii.set_qform(data.affine, code=1)
    raw_nii.set_sform(data.affine, code=1)
    nib.save(raw_nii, input_file)


def header_check(input_file):
    img = nib.load(input_file)
    axcodes = nib.aff2axcodes(img.affine)

    if axcodes != ("L", "I", "P"):
        print(
            f"Fatal header check failure: expected LIP orientation, found {axcodes} in {input_file}",
            file=sys.stderr,
        )
        sys.exit(FATAL_LIP_HEADER_EXIT_CODE)

    data = img.get_fdata(dtype=np.float32)

    out = nib.Nifti1Image(data, img.affine, header=img.header.copy())
    out.set_qform(img.affine, code=1)
    out.set_sform(img.affine, code=1)

    hdr = out.header
    hdr.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdr)

    nib.save(out, input_file)
    return input_file


def set_default_xyzt_units_if_unknown(target):
    if isinstance(target, str):
        img = nib.load(target)
        hdr = img.header
        set_default_xyzt_units_if_unknown(hdr)
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


def set_xform_codes_to_one(input_file):
    img = nib.load(input_file)
    img.set_qform(img.affine, code=1)
    img.set_sform(img.affine, code=1)
    nib.save(img, input_file)
    return input_file


def estimate_center_intensity_based(nifti, percentile=60):
    """
    Estimate BET center using intensity-weighted center-of-gravity, excluding
    low-intensity voxels with a data-adaptive percentile threshold.
    """
    p = subprocess.check_output(
        ["fslstats", nifti, "-P", str(percentile)]
    ).decode().strip()

    center = subprocess.check_output(
        ["fslstats", nifti, "-l", p, "-C"]
    ).decode().strip().split()

    cx, cy, cz = [float(v) for v in center]
    return [cx, cy, cz], float(p)


FSL_BET_WORLD_SWAPS = [(1, 2)]


def apply_world_ops(mat, swaps=()):
    out = mat.copy()

    for a, b in swaps:
        out[[a, b], :] = out[[b, a], :]

    return out


def save_header_only_reoriented_copy(src_path, dst_path, swaps=()):
    img = nib.load(src_path)
    data = img.get_fdata(dtype=np.float32)
    aff = img.affine.copy()

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


def applyBET(input_file,frac,radius,horizontal_gradient,
             use_bet4animal=False, species='mouse', center=None):
    """Apply BET"""
    if use_bet4animal:
        # Use BET for animal brains
        print("Using BET for animal brains")
        print("Note: bet4animal requires that the AC-PC line of brain is parallel to Y-axis")
        w_value = 2 #smooth the surface (lissencephalic weighting)
        species_id = 6 if species == 'mouse' else 5
        output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AnimalBet.nii.gz')

        tmp_bet = os.path.join(os.path.dirname(input_file), "bet4animal_tmp_out.nii.gz")
        tmp_mask = tmp_bet.replace(".nii.gz", "_mask.nii.gz")
        final_mask = output_file.replace(".nii.gz", "_mask.nii.gz")

        # ----- fslreorient2std -----
        tmp_std = os.path.join(os.path.dirname(input_file), "bet4animal_reorient2std.nii.gz")

        cmd = ["fslreorient2std", input_file, tmp_std]
        subprocess.run(cmd, check=True)

        # OPTIONAL: sanity prints
        # print("tmp_hdr axcodes:", nib.aff2axcodes(nib.load(tmp_hdr).affine))
        # print("tmp_std axcodes:", nib.aff2axcodes(nib.load(tmp_std).affine))

        bet_in = tmp_std

        #print("Header-only reorientation saved:", tmp_hdr)
        #print("New axcodes:", nib.aff2axcodes(aff))
        if center is None:
            center, p = estimate_center_intensity_based(bet_in)
        cx, cy, cz = center

        cmd = [
            "/aida/bin/bet4animal",
            bet_in,
            tmp_bet,
            "-m", #mask
            "-w", str(w_value),
            "-z", str(species_id),
            "-c", str(cx), str(cy), str(cz),
        ]
        subprocess.run(cmd, check=True)

        # ===== AFTER bet4animal =====
        #Nifti has to be reoriented to match the expected orientation and geometry of the AIDAmri pipeline (similar to real BET output) so that downstream steps remain compatible

        # ---------- (1) Reorient to original ----------
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

        # ---------- (2) Set affine ----------
        hdr_final = img.header.copy()
        hdr_final.set_data_dtype(np.float32)
        hdr_final["pixdim"][0] = 1
        hdr_final["pixdim"][4:8] = 1
        set_default_xyzt_units_if_unknown(hdr_final)

        img_final = nib.Nifti1Image(
            np.ascontiguousarray(data_lip, dtype=np.float32),
            aff_lip,
            header=hdr_final
        )
        # To match FSL BET output
        img_final.set_qform(aff_lip, code=1)
        img_final.set_sform(aff_lip, code=1)

        nib.save(img_final, output_file)

        #print("Final orientation:", nib.aff2axcodes(aff_final))
        #print("Final offset:", aff_final[:3, 3])

        # ===== APPLY SAME POST-PROCESSING TO BET MASK =====
        if os.path.exists(tmp_mask):
            # (1) Reorient to LIP
            m_img = nib.load(tmp_mask)
            m_data = m_img.get_fdata(dtype=np.float32)
            m_aff = m_img.affine

            m_ornt_cur = nib.orientations.io_orientation(m_aff)
            m_transform = nib.orientations.ornt_transform(m_ornt_cur, ornt_tgt)

            m_data_lip = nib.orientations.apply_orientation(m_data, m_transform)
            m_aff_lip = m_aff @ nib.orientations.inv_ornt_aff(m_transform, m_img.shape)

            # Binarize mask + uint8
            m_bin = (m_data_lip > 0.5).astype(np.uint8)

            m_hdr_lip = m_img.header.copy()
            m_hdr_lip.set_data_dtype(np.uint8)
            m_hdr_lip["pixdim"][0] = 1
            m_hdr_lip["pixdim"][4:8] = 1
            set_default_xyzt_units_if_unknown(m_hdr_lip)

            m_out = nib.Nifti1Image(
                np.ascontiguousarray(m_bin, dtype=np.uint8),
                m_aff_lip,
                header=m_hdr_lip
            )
            # To match FSL BET output
            m_out.set_qform(m_aff_lip, code=1)
            m_out.set_sform(m_aff_lip, code=1)

            nib.save(m_out, final_mask)
        else:
            print("Warning: BET mask not found:", tmp_mask)

        # remove temp
        try:
            for tmp_file in [tmp_std, tmp_bet, tmp_mask]:
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
        except Exception:
            pass
    #FSL BET (human modified version)
    else:
        data = nib.load(input_file)
        imgTemp = data.get_fdata()
        # create 4x4 scaling matrix and scale by 10 to match human like brain size
        scale = np.eye(4)
        scale[0, 0] = 10
        scale[1, 1] = 10
        scale[2, 2] = 10

        #Create new Nifti image with scaled affine
        scaled_affine = data.affine @ scale
        scaledNiiData = nib.Nifti1Image(imgTemp, scaled_affine)
        hdrIn = scaledNiiData.header
        set_default_xyzt_units_if_unknown(hdrIn)

        fslPath = os.path.join(os.path.dirname(input_file), 'fslScaleTemp.nii.gz')
        nib.save(scaledNiiData, fslPath)

        # temporary BET input with header-only LIP -> LPI world swap
        #This insures that the vertical gradient works in horizontally (so snout to cerebellum). This is needed bc this is a issue in FSL BET used for mice
        #Any questions regarding this ask Julian and hope he still knows
        bet_input_tmp = os.path.join(os.path.dirname(input_file), 'fslScaleTemp_LPIhdr.nii.gz')
        save_header_only_reoriented_copy(
            fslPath,
            bet_input_tmp,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # final output path
        output_file = os.path.join(
            os.path.dirname(input_file),
            os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz'
        )

        # temporary BET output in LPI-header space
        bet_output_tmp = os.path.join(
            os.path.dirname(input_file),
            os.path.basename(input_file).split('.')[0] + 'Bet_LPIhdr_tmp.nii.gz'
        )
        if center is not None:
            # nipype BET requires ints
            center_int = [int(round(c)) for c in center]
            print(f"Using user-defined center (rounded): {center_int}")
            myBet = fsl.BET(
                in_file=bet_input_tmp,
                out_file=bet_output_tmp,
                frac=frac,
                radius=radius,
                vertical_gradient=horizontal_gradient,
                center=center_int,
                mask=True
            )
        else:
            print("Using robust center estimation (-R)")
            myBet = fsl.BET(
                in_file=bet_input_tmp,
                out_file=bet_output_tmp,
                frac=frac,
                radius=radius,
                vertical_gradient=horizontal_gradient,
                robust=True,
                mask=True
            )
        myBet.run()

        # backswap BET image: LPI -> LIP (same swap again, because swap is its own inverse)
        save_header_only_reoriented_copy(
            bet_output_tmp,
            output_file,
            swaps=FSL_BET_WORLD_SWAPS,
        )

        # backswap BET mask: LPI -> LIP
        mask_tmp_file = bet_output_tmp.replace('.nii.gz', '_mask.nii.gz')
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')

        if os.path.exists(mask_tmp_file):
            save_header_only_reoriented_copy(
                mask_tmp_file,
                mask_file,
                swaps=FSL_BET_WORLD_SWAPS,
            )

        # unscale result data by factor 10ˆ(-1)
        dataOut = nib.load(output_file)
        imgOut = dataOut.get_fdata(dtype=np.float32)
        #rescale nifti
        inv_scale = np.eye(4)
        inv_scale[0, 0] = 0.1
        inv_scale[1, 1] = 0.1
        inv_scale[2, 2] = 0.1

        unscaled_affine = dataOut.affine @ inv_scale
        unscaledNiiData = nib.Nifti1Image(imgOut, unscaled_affine)
        unscaledNiiData.set_qform(unscaled_affine, code=1)
        unscaledNiiData.set_sform(unscaled_affine, code=1)
        hdrOut = unscaledNiiData.header
        set_default_xyzt_units_if_unknown(hdrOut)
        nib.save(unscaledNiiData, output_file)

        # also unscale BET mask
        mask_file = output_file.replace('.nii.gz', '_mask.nii.gz')
        if os.path.exists(mask_file):
            mask_data = nib.load(mask_file)
            bet_ref = nib.load(output_file)
            #make binary mask and apply affine of BET NIFTI
            mask_img = (mask_data.get_fdata() > 0.5).astype(np.uint8)

            finalMask = nib.Nifti1Image(mask_img, bet_ref.affine)
            finalMask.set_qform(bet_ref.affine, code=1)
            finalMask.set_sform(bet_ref.affine, code=1)

            hdrMask = finalMask.header
            hdrMask.set_data_dtype(np.uint8)
            set_default_xyzt_units_if_unknown(hdrMask)

            nib.save(finalMask, mask_file)
        # delete temporary files
        for tmp_file in [fslPath, bet_input_tmp, bet_output_tmp, mask_tmp_file]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
    print(f"Brain extraction completed, output saved to {output_file}")
    return output_file

def n4biasfieldcorr(input_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0] + 'AntsBias.nii.gz')
    myAnts = ants.N4BiasFieldCorrection(
        input_image=input_file,
        output_image=output_file,
        shrink_factor=4,
        bspline_fitting_distance=20,
        bspline_order=3,
        n_iterations=[1000, 0],
        dimension=3,
    )
    myAnts.run()

    img = nib.load(output_file)
    hdr = img.header
    hdr["pixdim"][4:8] = 1
    nib.save(img, output_file)

    print("Biasfield correction completed")
    return output_file


def skip_bet_function(input_file):
    """
    Create BET-compatible outputs when BET is skipped.
    """
    print("Skipping BET")
    print("Creating BET-compatible outputs for pipeline compatibility")

    outputBET = os.path.join(
        os.path.dirname(input_file),
        os.path.basename(input_file).split('.')[0] + 'Bet.nii.gz'
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
        header=hdr
    )
    final_img.set_qform(aff, code=1)
    final_img.set_sform(aff, code=1)
    nib.save(final_img, outputBET)

    print(f"BET skipped -> created compatibility image: {outputBET}")

    bet_mask_path = outputBET.replace('.nii.gz', '_mask.nii.gz')
    mask = (data > 0).astype(np.uint8)

    mask_hdr = hdr.copy()
    mask_hdr.set_data_dtype(np.uint8)
    mask_hdr["pixdim"][0] = 1
    mask_hdr["pixdim"][4:8] = 1
    set_default_xyzt_units_if_unknown(mask_hdr)

    mask_img = nib.Nifti1Image(
        np.ascontiguousarray(mask, dtype=np.uint8),
        aff,
        header=mask_hdr
    )
    mask_img.set_qform(aff, code=1)
    mask_img.set_sform(aff, code=1)
    nib.save(mask_img, bet_mask_path)

    print(f"BET mask created at {bet_mask_path}")

    return outputBET

def smoothIMG(input_file,outputPath):
    """
    Smoothes image via FSL. Only input and output has do be specified. Parameters are fixed to box shape and to the kernel size of 0.1 voxel.
    """
    source_base = os.path.basename(input_file).split('.')[0]
    data = nib.load(input_file)
    vol = data.get_fdata(dtype=np.float32)
    if vol.ndim == 4:
        ImgSmooth = np.min(vol, axis=3).astype(np.float32)
    elif vol.ndim == 3:
        ImgSmooth = vol.astype(np.float32)
    else:
        raise ValueError(f"Unsupported image dimensionality: {vol.ndim}")

    unscaledNiiData = nib.Nifti1Image(ImgSmooth, data.affine)
    unscaledNiiData.set_qform(data.affine, code=1)
    unscaledNiiData.set_sform(data.affine, code=1)
    hdrOut = unscaledNiiData.header
    hdrOut.set_data_dtype(np.float32)
    set_default_xyzt_units_if_unknown(hdrOut)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0] + 'DN.nii.gz')
    # hdrOut['sform_code'] = 1
    nib.save(unscaledNiiData, output_file)
    input_file = output_file
    #output_file =  os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0] + 'Smooth.nii.gz')
    output_file = os.path.join(outputPath, source_base + 'Smooth.nii.gz')
    myGauss =  fsl.SpatialFilter(in_file=input_file,out_file=output_file,operation='median',kernel_shape='box',kernel_size=0.1)
    myGauss.run()
    set_xform_codes_to_one(output_file)
    set_default_xyzt_units_if_unknown(output_file)
    print("Smoothing completed")
    return output_file

def thresh(input_file,outputPath):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]+ 'Thres.nii.gz')
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Thres.nii.gz')
    myThres = fsl.Threshold(in_file=input_file,out_file=output_file,thresh=20)#,direction='above')
    myThres.run()
    print("Thresholding completed")
    return output_file

def cropToSmall(input_file,outputPath):
    #output_file = os.path.join(os.path.dirname(input_file),os.path.basename(input_file).split('.')[0]  + 'Crop.nii.gz')
    output_file = os.path.join(outputPath, os.path.basename(input_file).split('.')[0] + 'Crop.nii.gz')
    myCrop = fsl.ExtractROI(in_file=input_file,roi_file=output_file,x_min=40,x_size=130,y_min=50,y_size=110,z_min=0,z_size=12)
    myCrop.run()
    print("Cropping done")
    return  output_file


if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='Preprocessing of rsfMRI Data')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the RAW data of rsfMRI NIfTI file', required=True)

    parser.add_argument('-f', '--frac',
                        help='Fractional intensity threshold - default=0.3, smaller values give larger brain outline estimates',
                        nargs='?', type=float, default=0.15)
    parser.add_argument('-r', '--radius', help='Head radius (mm not voxels) - default=45', nargs='?', type=int ,default=45)
    parser.add_argument(
        '-g',
        '--horizontal_gradient',
        '--vertical_gradient',
        dest='horizontal_gradient',
        help='Horizontal gradient in fractional intensity threshold - default=0.0. Not for bet4animals! Higher positive values make the BET stricter posterior and less stricter anterior (snout)',
        nargs='?',
        type=float,
        default=0.0,
    )
    parser.add_argument(
        '-c', '--center',
        help='Brain center in voxel coordinates: x y z',
        nargs=3,
        type=float,
        default=None
    )
    parser.add_argument(
        '--use_bet4animal',
        help='Use BET for animal brains. If not set it uses FSL BET.',
        action='store_true'
    )
    parser.add_argument(
        '--bet_skip',
        help='Skip BET during fMRI preprocessing (still creates *Bet.nii.gz and *_mask.nii.gz for pipeline compatibility). '
             'If not set it uses FSL BET (modified human version)',
        action='store_true'
    )
    parser.add_argument(
        '-b',
        '--bias_method',
        help='Biasfield correction method - default=None, other options are "ants" or "none"',
        choices=["none", "ants"],
        type=str.lower,
        default=None,
    )
    args = parser.parse_args()

    # set parameters
    inputFile = None
    if args.input is not None and args.input is not None:
        inputFile = args.input

    if not os.path.exists(inputFile):
        sys.exit(f"Error: input file does not exist: {inputFile}")

    frac = args.frac
    radius = args.radius
    horizontal_gradient = args.horizontal_gradient
    bias_method = args.bias_method
    outputPath = os.path.dirname(inputFile)

    print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    create_brkraw_backup(inputFile)
    header_check(inputFile)

    outputSmooth = smoothIMG(input_file=inputFile,outputPath=outputPath)

    if bias_method is None or bias_method == "none":
        print("No bias field correction applied")
        outputBiasCorr = outputSmooth
    elif bias_method == "ants":
        print("Starting Biasfieldcorrection with ANTS:")
        try:
            outputBiasCorr = n4biasfieldcorr(input_file=outputSmooth)
            set_xform_codes_to_one(outputBiasCorr)
            set_default_xyzt_units_if_unknown(outputBiasCorr)
            print("Biasfield correction was successful")
        except Exception as e:
            print(f'Error in bias field correction\nError message: {str(e)}')
            raise

    if args.bet_skip:
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(outputBiasCorr)
    else:
        # get rid of your skull
        outputBET = applyBET(
            input_file=outputBiasCorr,
            frac=frac,
            radius=radius,
            horizontal_gradient=horizontal_gradient,
            use_bet4animal=args.use_bet4animal,
            center=args.center,
        )

    print("Preprocessing completed")
