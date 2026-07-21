"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""


import sys, os
import nipype.interfaces.fsl as fsl
import nibabel as nii
import numpy as np
import glob
import shutil
import regress
import getSingleRegTable
import cv2
import create_seed_rois
import fsl_mean_ts
from pathlib import Path 
import json

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


def copyAtlasOfData(path,post,labels):
    fileALL = glob.glob(path + '/*' + post + '.nii.gz')
    if fileALL.__len__()>1:
        sys.exit("Error: '%s' has no related Atlas File." % (path,))
    else:
        fileALL = fileALL[0]

    print("Copy Atlas Data and generate seed ROIs")
    #pathfMRI = os.path.join(os.path.dirname(path),'fMRI')
    outputRois = create_seed_rois.startSeedPoint(in_atlas=os.path.join(path, os.path.basename(fileALL)),in_labels=labels)
    return outputRois

def imgScaleResize(img):
    newImg = np.zeros([128,128,20,355])
    for i in range(img.shape[3]):
        for j in range(img.shape[2]):
            newImg[:, :, j, i] = cv2.resize(
                img[:, :, j, i],
                (newImg.shape[1], newImg.shape[0]),
                interpolation=cv2.INTER_NEAREST
            )

    return newImg

def scaleBy10(input_path,inv):
    img = nii.load(input_path)
    imgTemp = img.get_fdata(dtype=np.float32)
    aff = img.affine.copy()

    factor = 0.1 if inv else 10.0
    aff[:3, :3] *= factor

    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    out_img = nii.Nifti1Image(imgTemp, aff, header=hdr)
    out_img.header.set_xyzt_units('mm')
    out_img.set_qform(aff, code=1)
    out_img.set_sform(aff, code=1)

    if inv is False:
        fslPath = os.path.join(
            os.path.dirname(input_path),
            os.path.basename(input_path).split('.')[0] + "_fslScaleTemp.nii.gz"
        )
        nii.save(out_img, fslPath)
        return fslPath
    elif inv is True:
        nii.save(out_img, input_path)
        return input_path
    else:
        sys.exit("Error: inv - parameter should be a boolean.")

def findSlicesData(path,pre):
    regMR_list = []
    fileALL = glob.iglob(path+'/'+pre+'*.nii.gz',recursive=True)
    for filename in fileALL:
        regMR_list.append(filename)
    regMR_list.sort()
    return regMR_list

def copyEPIToProcessingFolder(file_name, proc_Path):
    output_file = os.path.join(proc_Path, os.path.basename(file_name))
    shutil.copyfile(file_name, output_file)
    return output_file


def _strip_nii_gz(path):
    name = os.path.basename(path)
    if name.endswith(".nii.gz"):
        return name[:-7]
    return os.path.splitext(name)[0]


def _is_bet_file(path):
    name = os.path.basename(path)
    return name.endswith("Bet.nii.gz") and not name.endswith("_mask.nii.gz")


def findExistingBET(raw_file, bet_file=None):
    origin_path = os.path.dirname(raw_file)

    if bet_file is not None:
        if not os.path.exists(bet_file):
            sys.exit(f"Error: BET file does not exist: {bet_file}")
        if not _is_bet_file(bet_file):
            sys.exit(f"Error: BET file must end with Bet.nii.gz: {bet_file}")
        return bet_file

    if _is_bet_file(raw_file):
        return raw_file

    candidates = sorted(
        file_name for file_name in glob.glob(os.path.join(origin_path, "*Bet.nii.gz"))
        if _is_bet_file(file_name)
    )

    if len(candidates) == 0:
        sys.exit(
            "Error: No existing BET file found in '%s'. "
            "Run preProcessing_fMRI.py first or pass --bet-file." % (origin_path,)
        )

    raw_stem = _strip_nii_gz(raw_file)
    preferred = [
        candidate for candidate in candidates
        if _strip_nii_gz(candidate).startswith(raw_stem)
    ]

    if len(preferred) == 1:
        return preferred[0]

    if len(preferred) > 1:
        candidate_list = "\n  ".join(preferred)
        sys.exit(
            "Error: Multiple matching BET files found for '%s':\n  %s\n"
            "Pass the intended file with --bet-file." % (raw_file, candidate_list)
        )

    if len(candidates) == 1:
        return candidates[0]

    candidate_list = "\n  ".join(candidates)
    sys.exit(
        "Error: Multiple BET files found in '%s':\n  %s\n"
        "Pass the intended file with --bet-file." % (origin_path, candidate_list)
    )


def resolveRawEPIFromBETInput(bet_file):
    origin_path = os.path.dirname(bet_file)
    candidates = sorted(glob.glob(os.path.join(origin_path, "*EPI.nii.gz")))

    if len(candidates) == 0:
        sys.exit(
            "Error: Input is a BET file, but no raw *EPI.nii.gz file was found in '%s'."
            % (origin_path,)
        )

    bet_stem = _strip_nii_gz(bet_file)
    preferred = [
        candidate for candidate in candidates
        if bet_stem.startswith(_strip_nii_gz(candidate))
    ]

    if len(preferred) == 1:
        return preferred[0]

    if len(preferred) > 1:
        candidate_list = "\n  ".join(preferred)
        sys.exit(
            "Error: Multiple raw EPI files match BET file '%s':\n  %s"
            % (bet_file, candidate_list)
        )

    if len(candidates) == 1:
        return candidates[0]

    candidate_list = "\n  ".join(candidates)
    sys.exit(
        "Error: Input is a BET file, but the raw EPI file is ambiguous in '%s':\n  %s"
        % (origin_path, candidate_list)
    )


def copyExistingBETToProcessingFolder(bet_file, proc_Path):
    mask_file = bet_file.replace(".nii.gz", "_mask.nii.gz")
    if not os.path.exists(mask_file):
        sys.exit(
            "Error: Existing BET mask is missing:\n  %s\n"
            "The processing step reuses the BET mask from preprocessing and no longer creates a new one."
            % (mask_file,)
        )

    copied_bet = os.path.join(proc_Path, os.path.basename(bet_file))
    copied_mask = os.path.join(proc_Path, os.path.basename(mask_file))
    shutil.copyfile(bet_file, copied_bet)
    shutil.copyfile(mask_file, copied_mask)
    print("Copied existing BET to %s" % (copied_bet,))
    print("Copied existing BET mask to %s" % (copied_mask,))
    return copied_bet, copied_mask


def getEPIMean(file_name,proc_Path):
    output_file = os.path.join(proc_Path, os.path.basename(file_name).split('.')[0]) + 'mean.nii.gz'
    myMean = fsl.MeanImage(in_file=file_name, out_file=output_file)
    print(myMean.cmdline)
    myMean.run()
    return output_file

def applyMask(input_file,mask_file):
    output_file = os.path.join(os.path.dirname(input_file), os.path.basename(input_file).split('.')[0]) + 'BET.nii.gz'
    myMaskapply = fsl.ApplyMask(in_file=input_file, out_file=output_file, mask_file=mask_file)
    print(myMaskapply.cmdline)
    myMaskapply.run()
    return output_file

def fsl_SeparateSliceMoCo(input_file,par_folder):
    # scale Nifti data by factor 10
    dataName = os.path.basename(input_file).split('.')[0]

    aidamri_dir = os.getcwd()
    temp_dir = os.path.join(os.path.dirname(input_file), "temp")
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    fslPath = scaleBy10(input_file, inv=False)
    os.chdir(temp_dir)
    mySplit= fsl.Split(in_file=fslPath,dimension='z',out_base_name = dataName)
    print(mySplit.cmdline)
    mySplit.run()
    os.remove(fslPath)


    # sparate ref and src volume in slices
    sliceFiles = findSlicesData(os.getcwd(),dataName)

    #start to correct motions slice by slice
    for i in  range(len(sliceFiles)):
        slc = sliceFiles[i]
        # ref = refFiles[i]
        # take epi as ref
        output_file = os.path.join(par_folder,os.path.basename(slc))
        myMCFLIRT = fsl.preprocess.MCFLIRT(in_file=slc,out_file=output_file,save_plots=True,terminal_output='none')
        myMCFLIRT.run()
        os.remove(slc)
        # os.remove(ref)

    # merge slices to a single volume

    mcf_sliceFiles = findSlicesData(par_folder,dataName)
    output_file = os.path.join(os.path.dirname(input_file),
                               os.path.basename(input_file).split('.')[0]) + '_mcf.nii.gz'
    myMerge = fsl.Merge(in_files=mcf_sliceFiles,dimension='z',merged_file =output_file)
    print(myMerge.cmdline)
    myMerge.run()

    for slc in mcf_sliceFiles: os.remove(slc)

    # unscale result data by factor 10ˆ(-1)
    output_file = scaleBy10(output_file, inv=True)
    
    #os.remove(temp_dir)
    os.chdir(aidamri_dir)

    return output_file

def copyRawPhysioData(file_name,i32_Path):
    img_name = Path(file_name).name
    json_name = img_name.replace(".nii.gz", ".json")
    
    json_file = os.path.join(os.path.dirname(file_name), json_name)
    sub_name = (Path(file_name).name.split("_")[0]).split("-")[1]
    studyName = (Path(os.path.dirname(os.path.dirname(file_name))).name).split("-")[1]
    physioPath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(file_name)))),'Physio')
    scanid = None
    
    relatedPhysioData = []
    if not os.path.exists(json_file):
        print("Error: '%s' has no metadata JSON file for physio lookup." % (file_name,))
        return []

    with open(json_file, 'r') as infile:
        content = json.load(infile)

    if "ScanID" not in content:
        print("Error: '%s' has no ScanID in metadata JSON for physio lookup." % (json_file,))
        return []

    scanid = str(content["ScanID"]) + ".I32"

    if not os.path.exists(physioPath):
        print("Error: '%s' is not an existing Physio directory." % (physioPath,))
        return []

    conditions = [sub_name , studyName]
    for file in glob.iglob(os.path.join(physioPath, "**", "*" + scanid), recursive=True):
        filename = os.path.basename(file)
        if all(condition in filename for condition in conditions):
            relatedPhysioData.append(file)

    if len(relatedPhysioData)>1:
        sys.exit("Warning: '%s' has no unique physio data for scan %s." % (physioPath, scanid,))
    if len(relatedPhysioData) == 0:
        print("Error: '%s' has no related physio data for scan %s." % (physioPath, scanid,))
        return []

    physioFile_name = relatedPhysioData[0]
    print('Copy related physio data %s to rawMoData' % (physioFile_name,))
    shutil.copyfile(physioFile_name, os.path.join(i32_Path,os.path.basename(physioFile_name)))

    return physioFile_name

def create_txt_file(file, data):
    import time
    with open(file, "w") as outfile:
        for data_point in data:
            outfile.write(''.join([str(data_point), '\n']))

    time.sleep(1.0)

def delete_txt_file(file):
    os.remove(file)

def startProcess(
    Rawfile_name,
    bet_file=None,
):
    # generate folder for images
    origin_Path = os.path.dirname(Rawfile_name)
    proc_Path = os.path.join(origin_Path, 'rs-fMRI_niiData')
    if os.path.exists(proc_Path):
        shutil.rmtree(proc_Path)
    os.mkdir(proc_Path)

    # generate folder for motion correction files
    par_Path = os.path.join(origin_Path, 'rs-fMRI_mcf')
    if os.path.exists(par_Path):
        shutil.rmtree(par_Path)
    os.mkdir(par_Path)

    # generate folder for motion correction files
    subFile= os.path.basename(Rawfile_name).split('.')[0]
    subFile = '%s_mcf.mat' % subFile
    par_Path = os.path.join(par_Path,subFile)
    if os.path.exists(par_Path):
        shutil.rmtree(par_Path)
    os.mkdir(par_Path)

    # generate folder for physio data
    i32_Path = os.path.join(origin_Path, 'rawMonData')
    if  os.path.exists(i32_Path):
        shutil.rmtree(i32_Path)
    os.mkdir(i32_Path)

    # copy raw EPI without changing voxel array orientation
    file_name = copyEPIToProcessingFolder(Rawfile_name,proc_Path)

    # calculate EPIMean
    file_nameEPI = getEPIMean(file_name,proc_Path)

    # reuse the existing BET from fMRI preprocessing; keep originals in func unchanged
    existing_bet_file = findExistingBET(Rawfile_name, bet_file=bet_file)
    _, mask_file = copyExistingBETToProcessingFolder(existing_bet_file, proc_Path)

    #apply Mask on original dataset
    applyMask(file_name,mask_file)

    # apply motion correction on original dataset with EPImean as reference
    mcfFile_name=fsl_SeparateSliceMoCo(file_name,par_Path)

    # apply mean on motion corrected data
    getEPIMean(mcfFile_name, proc_Path)

    # copy physio data to rawMonData-Folder
    relatedPhysioFolder = copyRawPhysioData(Rawfile_name,i32_Path)

    # get Regression Values
    if len(relatedPhysioFolder) != 0:
        getSingleRegTable.getRegrTable(os.path.dirname(Rawfile_name),relatedPhysioFolder,par_Path)
    else:
        print("Error: Processing not possible, because either there is no folder called Physio or the related physio data for the scan is missing there.")

    return mcfFile_name, mask_file

if __name__ == "__main__":

    TR = 1.0 # all samples have TR = 1s
    FWHM = 8.0 #bc of rat

    import argparse
    parser = argparse.ArgumentParser(description='Process fMRI data')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--input', help='Path to the BET NIFTI', required=True)

    parser.add_argument('-t', '--TR', default=TR, type=float, help='Current TR value')
    parser.add_argument('-f', '--FWHM', default=FWHM, type=float, help='Full width at half maximum')
    parser.add_argument('-stc', '--slicetimecorrection', default="False", type=str, help='choose to perform slice time correction or not')
    parser.add_argument('--bet-file', default=None, help='Existing func/*Bet.nii.gz file to reuse when auto-detection is ambiguous')

    args = parser.parse_args()

    TR = args.TR
    FWHM = args.FWHM

    if args.slicetimecorrection == "True":
        stc = True
    else:
        stc = False

    sigma_labels = os.path.join(REPO_ROOT, 'lib', 'sigma', 'SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt')
    labels = sigma_labels
    labelNames = sigma_labels
    labels2000 = sigma_labels
    labelNames2000 = sigma_labels
    input_file = None
    if args.input is not None and args.input is not None:
        input_file = args.input
    if not os.path.exists(input_file):
        sys.exit(f"Error: input file does not exist: {input_file}")

    bet_file = args.bet_file
    if _is_bet_file(input_file):
        if bet_file is not None and os.path.abspath(bet_file) != os.path.abspath(input_file):
            sys.exit("Error: Use either -i with a BET file or --bet-file, not both.")
        bet_file = input_file
        input_file = resolveRawEPIFromBETInput(bet_file)
        print("Using raw EPI file for processing: %s" % (input_file,))

    mcfFile_name, mask_file = startProcess(
        input_file,
        bet_file=bet_file
    )

    
    # if stc is activated find parameters
    if stc: 
        print("Starting Regression with slice time correction:")
        # find meta data json file
        meta_data_file_name = Path(input_file).name.replace(".nii.gz", ".json")
        meta_data_file = os.path.join(Path(input_file).parent, meta_data_file_name)

        with open(meta_data_file, "r") as infile:
            meta_data = json.load(infile)


        slice_timing = meta_data["SliceTiming"]
        slice_order = [index + 1 for index in sorted(range(len(slice_timing)), key=lambda index: slice_timing[index])]
        n_slices = len(slice_timing)
        costum_timings = [timing / TR for timing in slice_timing]

        # create costum timings txt file
        costum_timings_path = os.path.join(Path(meta_data_file).parent, "tcostum.txt")
        create_txt_file(costum_timings_path, costum_timings)
        
        # create slice order txt file
        slice_order_path = os.path.join(Path(meta_data_file).parent, "slice_order.txt")
        create_txt_file(slice_order_path, slice_order)
        
        rgr_file, srgr_file, sfrgr_file = regress.startRegression(
            mcfFile_name,
            FWHM=FWHM,
            TR=TR,
            stc=stc,
            slice_order=slice_order_path,
            costum_timings=costum_timings_path,
            mask_file=mask_file
        )

        # delete temp txt files
        delete_txt_file(costum_timings_path)
        delete_txt_file(slice_order_path)

    else:
        print("Starting Regression without slice time correction:")
        rgr_file, srgr_file, sfrgr_file = regress.startRegression(
            mcfFile_name,
            FWHM=FWHM,
            TR=TR,
            stc=stc,
            mask_file=mask_file
        )
        print(f"sfrgr_file {sfrgr_file}")

    

    atlasPath = os.path.dirname(input_file)
    roisPath = copyAtlasOfData(atlasPath,'Anno_parental',labels)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames, 'MasksTCs.')

    roisPath = copyAtlasOfData(atlasPath, 'AnnoSplit_parental', labels2000)

    fslMeantsFile = fsl_mean_ts.start_fsl_mean_ts(sfrgr_file, roisPath, labelNames2000, 'MasksTCsSplit.')
