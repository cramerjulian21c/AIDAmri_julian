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
import shutil
#makes sure to import bet.py
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from common.bet import applyBET, skip_bet_function
from common.script_logging import setup_script_logging

FATAL_LIP_HEADER_EXIT_CODE = 86
SUPERIOR_NOISE_DEPTH = 2


def create_brkraw_backup(input_file):

    brkraw_dir = os.path.join(os.path.dirname(input_file), "brkraw")
    if os.path.exists(brkraw_dir):
        return 

    os.mkdir(brkraw_dir)
    dst_path = os.path.join(brkraw_dir, os.path.basename(input_file))

    shutil.copyfile(input_file, dst_path)

    data = nib.load(input_file)
    # Preserve nibabel scaling (scl_slope/scl_inter). Using get_unscaled()
    # would write raw int16 scanner values into the working file.
    raw_img = data.get_fdata(dtype=np.float32)

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


def suppress_superior_noise(input_file, edge_depth=SUPERIOR_NOISE_DEPTH):
    """
    Remove a thin foreground layer along the anatomical superior BET edge.

    Selected voxels are set to zero in the existing BET image. Two additional
    QC images are written: ``SuperiorNoiseMask`` marks removed voxels with 1,
    while ``SuperiorNoiseWeight`` marks retained voxels with 1.
    """
    if edge_depth < 1:
        raise ValueError("Superior noise depth must be at least one voxel.")

    # This function must run on the 3D reference image returned by applyBET,
    # not on the original 4D fMRI time series.
    image = nib.load(input_file)
    data = image.get_fdata(dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(
            "Superior-noise suppression requires a 3D BET image, got shape "
            f"{data.shape} for {input_file}."
        )

    # BET normally represents background as zero. Positive finite voxels are
    # therefore treated as the current foreground from which to remove a layer.
    foreground = np.isfinite(data) & (data > 0)
    if not np.any(foreground):
        raise ValueError(f"BET image contains no positive foreground voxels: {input_file}")

    # Determine the anatomical superior-inferior axis
    axis_codes = nib.aff2axcodes(image.affine)
    try:
        superior_axis = next(
            axis for axis, code in enumerate(axis_codes)
            if code in ("S", "I")
        )
    except StopIteration as error:
        raise ValueError(
            f"Could not determine the superior-inferior axis of {input_file} "
            f"(orientation: {axis_codes})."
        ) from error

    # Temporarily move the superior-inferior dimension to axis 0 and orient it
    # so index 0 always corresponds to the anatomical superior side.
    superior_first = np.moveaxis(foreground, superior_axis, 0)
    # Nifti has to be in LIP orientation!!

    # The cumulative foreground count is calculated independently for every
    # column through the brain. Selecting counts 1..edge_depth makes the mask
    # follow the curved superior BET surface rather than a flat array slice.
    noise_superior_first = superior_first & (
        np.cumsum(superior_first, axis=0, dtype=np.int32) <= edge_depth
    )

    noise_mask = np.moveaxis(noise_superior_first, 0, superior_axis)

    # Mask conventions:
    #   noise_mask: 1 = remove, 0 = retain
    #   keep_mask:  0 = remove, 1 = retain
    keep_mask = np.ones(data.shape, dtype=np.float32)
    keep_mask[noise_mask] = 0.0

    prefix = os.path.basename(input_file).split('.')[0]
    output_dir = os.path.dirname(input_file)
    noise_mask_path = os.path.join(output_dir, prefix + '_SuperiorNoiseMask.nii.gz')
    noise_weight_path = os.path.join(output_dir, prefix + '_SuperiorNoiseWeight.nii.gz')

    # Save both conventions so the removed area can be inspected directly and
    # the inverse mask can be reused without another inversion operation.
    weight_header = image.header.copy()
    weight_header.set_data_dtype(np.float32)
    weight_image = nib.Nifti1Image(keep_mask, image.affine, weight_header)
    weight_image.set_qform(image.affine, code=1)
    weight_image.set_sform(image.affine, code=1)
    nib.save(weight_image, noise_weight_path)

    mask_header = image.header.copy()
    mask_header.set_data_dtype(np.uint8)
    mask_image = nib.Nifti1Image(
        noise_mask.astype(np.uint8),
        image.affine,
        mask_header,
    )
    mask_image.set_qform(image.affine, code=1)
    mask_image.set_sform(image.affine, code=1)
    nib.save(mask_image, noise_mask_path)

    # Apply the keep mask directly to the BET result. Geometry, dimensions and
    # filename remain unchanged, so downstream code still finds *Bet.nii.gz.
    cleaned_data = np.ascontiguousarray(data * keep_mask, dtype=np.float32)
    cleaned_header = image.header.copy()
    cleaned_header.set_data_dtype(np.float32)
    cleaned_image = nib.Nifti1Image(cleaned_data, image.affine, cleaned_header)
    cleaned_image.set_qform(image.affine, code=1)
    cleaned_image.set_sform(image.affine, code=1)
    nib.save(cleaned_image, input_file)

    # process_fMRI.py later applies the companion BET mask to the 4D time
    # series. Update it as well so preprocessing, registration and processing
    # all use the same cleaned brain extent.
    bet_mask_path = input_file.replace('.nii.gz', '_mask.nii.gz')
    if os.path.exists(bet_mask_path):
        bet_mask_image = nib.load(bet_mask_path)
        bet_mask_data = bet_mask_image.get_fdata() > 0
        if bet_mask_data.shape != keep_mask.shape:
            raise ValueError(
                f"BET image and BET mask shapes differ: {data.shape} vs "
                f"{bet_mask_data.shape}."
            )
        #Keep voxel if its inside BET mask and not inside noise mask
        cleaned_bet_mask = bet_mask_data & ~noise_mask
        bet_mask_header = bet_mask_image.header.copy()
        bet_mask_header.set_data_dtype(np.uint8)
        cleaned_bet_mask_image = nib.Nifti1Image(
            cleaned_bet_mask.astype(np.uint8),
            bet_mask_image.affine,
            bet_mask_header,
        )
        cleaned_bet_mask_image.set_qform(bet_mask_image.affine, code=1)
        cleaned_bet_mask_image.set_sform(bet_mask_image.affine, code=1)
        nib.save(cleaned_bet_mask_image, bet_mask_path)

    print(
        "Superior-noise suppression completed after BET: "
        f"orientation={axis_codes}, depth={edge_depth}, "
        f"removed_voxels={int(noise_mask.sum())}, BET={input_file}, "
        f"noise_mask={noise_mask_path}, keep_mask={noise_weight_path}"
    )
    return input_file


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


def smoothIMG(input_file, outputPath, skip_smoothing=False):
    """
    Prepare a 3D reference image and optionally apply FSL's median spatial filter.
    For 4D inputs, a voxel-wise median projection across the 4th dimension is
    written as *MP.nii.gz. For 3D inputs, the MP image is just a
    float32/header-normalized copy.
    """
    source_base = os.path.basename(input_file).split('.')[0]
    data = nib.load(input_file)
    vol = data.get_fdata(dtype=np.float32)
    if vol.ndim == 4:
        ImgSmooth = np.median(vol, axis=3).astype(np.float32)
        source_base = source_base + 'MP'
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
                               os.path.basename(input_file).split('.')[0] + 'MP.nii.gz')
    # hdrOut['sform_code'] = 1
    nib.save(unscaledNiiData, output_file)
    input_file = output_file

    if skip_smoothing:
        print("Spatial smoothing skipped")
        return input_file

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
                        help='Fractional intensity threshold - default=0.1, smaller values give larger brain outline estimates',
                        nargs='?', type=float, default=0.1)
    parser.add_argument('-r', '--radius', help='Head radius (mm not voxels) - default=60', nargs='?', type=int ,default=60)
    parser.add_argument(
        '-g',
        '--horizontal-gradient',
        help='Horizontal gradient in fractional intensity threshold - default=0.13. Not for bet4animals! Higher positive values make the BET stricter posterior and less stricter anterior (snout)',
        nargs='?',
        type=float,
        default=0.13,
    )
    parser.add_argument(
        '-c', '--center',
        help='Brain center in voxel coordinates: x y z',
        nargs=3,
        type=float,
        default=None
    )
    parser.add_argument(
        '--bet',
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help='Brain extraction method for fMRI preprocessing: skip, bet or bet4animal. Default: bet'
    )
    parser.add_argument(
        '-b',
        '--bias-method',
        help='Biasfield correction method - default=None, other options are "ants" or "skip"',
        choices=["skip", "ants"],
        type=str.lower,
        default=None,
    )
    parser.add_argument(
        '--skip-smoothing',
        action='store_true',
        help='Skip the FSL spatial median smoothing step; still creates the 3D median reference image'
    )
    parser.add_argument(
        '--suppress-superior-noise',
        action='store_true',
        help='remove the superior foreground edge after BET',
    )
    parser.add_argument(
        '--superior-noise-depth',
        type=int,
        default=SUPERIOR_NOISE_DEPTH,
        help='number of superior foreground voxels removed per image column after BET (default: %(default)s)',
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
    setup_script_logging(outputPath, "preprocess.log")

    if args.bet == "bet":
        print(f"Frac: {frac} Radius: {radius} Gradient {horizontal_gradient}")

    create_brkraw_backup(inputFile)
    header_check(inputFile)

    outputSmooth = smoothIMG(input_file=inputFile, outputPath=outputPath, skip_smoothing=args.skip_smoothing)

    if bias_method is None or bias_method == "skip":
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

    if args.bet == "skip":
        print("Skipping brain extraction.")
        outputBET = skip_bet_function(outputBiasCorr)
    else:
        # get rid of your skull
        outputBET = applyBET(
            input_file=outputBiasCorr,
            frac=frac,
            radius=radius,
            horizontal_gradient=horizontal_gradient,
            use_bet4animal=args.bet == "bet4animal",
            center=args.center,
        )

    # Apply the optional cleanup only after BET has established the foreground
    # and created the companion BET mask used by downstream fMRI processing.
    if args.suppress_superior_noise:
        outputBET = suppress_superior_noise(
            outputBET,
            edge_depth=args.superior_noise_depth,
        )

    print("Preprocessing completed")
