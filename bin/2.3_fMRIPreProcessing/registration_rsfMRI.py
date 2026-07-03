"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""

import sys,os
import glob
import shutil as sh
import subprocess
import shlex
import logging


LOGGER = logging.getLogger(__name__)
DISABLE_LOG_ENV = "AIDAMRI_DISABLE_SCRIPT_LOG"


def setup_logging(outfile):
    handlers = [logging.StreamHandler()]
    if os.environ.get(DISABLE_LOG_ENV) != "1":
        handlers.append(logging.FileHandler(os.path.join(outfile, "registration.log"), mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=handlers,
    )


def require_single_match(matches, description):
    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one %s, found %d: %s"
            % (description, len(matches), matches)
        )
    return matches[0]


def regSIG2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitAnno, splitAnno_rsfMRI, anno_rsfMRI,
                  bsplineMatrix, dref, outfile, use_atlas_mask=False):
    outputT2w = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_T2w.nii.gz')
    outputAff = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + 'transMatrixAff.txt')

    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*T2w.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI T2 image"), outputT2w)
    else:
        registrationT2 = T2data
        if use_atlas_mask:
            prefix = os.path.basename(inputVolume).split('.')[0]
            atlasMask = os.path.join(outfile, prefix + '_T2AtlasMask.nii.gz')
            maskedT2 = os.path.join(outfile, prefix + '_T2AtlasMasked.nii.gz')

            # Create a slightly dilated mask from the well-registered atlas
            # annotation and apply it to the imperfect T2 BET.
            for command_args in [
                ["fslmaths", brain_anno, "-bin", "-dilM", atlasMask], #binarize atlas and dilate it by one voxel layer for tolerance
                ["fslmaths", T2data, "-mas", atlasMask, maskedT2], #apply mask to the T2 BET
            ]:
                result = subprocess.run(
                    command_args,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                LOGGER.info(
                    "Output of %s:\n%s",
                    subprocess.list2cmdline(command_args),
                    result.stdout,
                )
            registrationT2 = maskedT2

        command = f"reg_aladin -ref {inputVolume} -flo {registrationT2} -res {outputT2w} -rigOnly -aff {outputAff}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
        #  resample Annotation
        outputAnno = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno.nii.gz')

        command = f"reg_resample -ref {inputVolume} -flo {brain_anno} -cpp {outputAff} -inter 0 -res {outputAnno}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise

    # resample split annotation
    outputAnnoSplit = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI split annotation"), outputAnnoSplit)
    else:
        outputAnnoSplit_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_T2w.nii.gz')
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_T2}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise

        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_T2} -trans {outputAff} -inter 0 -res {outputAnnoSplit}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
        os.remove(outputAnnoSplit_T2)

    # resample split parental annotation
    outputAnnoSplit_rsfMRI = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*AnnoSplit_parental.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI parental split annotation"), outputAnnoSplit_rsfMRI)
    else:
        outputAnnoSplit_rsfMRI_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_AnnoSplit_parental_T2w.nii.gz')
        command = f"reg_resample -ref {brain_anno} -flo {splitAnno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnnoSplit_rsfMRI_T2}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
        
        command = f"reg_resample -ref {inputVolume} -flo {outputAnnoSplit_rsfMRI_T2} -trans {outputAff} -inter 0 -res {outputAnnoSplit_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
        os.remove(outputAnnoSplit_rsfMRI_T2)

    # resample parental annotation
    outputAnno_rsfMRI = os.path.join(outfile,
                                          os.path.basename(inputVolume).split('.')[0] + '_Anno_parental.nii.gz')
    if dref:
        pathT2 = glob.glob(os.path.dirname(outfile) + '*/dwi/*Anno_parental.nii.gz', recursive=False)
        sh.copy(require_single_match(pathT2, "DTI parental annotation"), outputAnno_rsfMRI)
    else: 
        outputAnno_rsfMRI_T2 = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Anno_parental_T2w.nii.gz')
        command = f"reg_resample -ref {brain_anno} -flo {anno_rsfMRI} -trans {bsplineMatrix} -inter 0 -res {outputAnno_rsfMRI_T2}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise

        command = f"reg_resample -ref {inputVolume} -flo {outputAnno_rsfMRI_T2} -trans {outputAff} -inter 0 -res {outputAnno_rsfMRI}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
        os.remove(outputAnno_rsfMRI_T2)
        
        # resample in-house developed tempalate
        outputTemplate = os.path.join(outfile, os.path.basename(inputVolume).split('.')[0] + '_Template.nii.gz')
        
        command = f"reg_resample -ref {inputVolume} -flo {brain_template} -trans {outputAff} -res {outputTemplate}"
        command_args = shlex.split(command)
        try:
            result = subprocess.run(command_args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,text=True)
            LOGGER.info("Output of %s:\n%s", command, result.stdout)
        except Exception as e:
            LOGGER.error("Error while executing the command: %s\nErrorcode: %s", command_args, e)
            raise
    
    return outputAnnoSplit

def find_RefStroke(refStrokePath,inputVolume):
    search_patterns = [
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9] + '*', 'anat', 'IncidenceData', 'IncidenceData_Lesion_mask.nii.gz'),
        os.path.join(refStrokePath, os.path.basename(inputVolume)[0:9] + '*', 'anat', '*IncidenceData_mask.nii.gz'),
    ]
    path = []
    for pattern in search_patterns:
        path.extend(glob.glob(pattern, recursive=False))
    return path

def find_RefAff(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/anat/*MatrixAff.txt', recursive=False)
    return path

def find_RefTemplate(inputVolume):
    path =  glob.glob(os.path.dirname(os.path.dirname(inputVolume))+'/anat/*TemplateAff.nii.gz', recursive=False)
    return path


def find_relatedData(pathBase):
    pathT2 =  glob.glob(pathBase+'*/anat/*Bet.nii.gz', recursive=False)
    pathStroke_mask = glob.glob(pathBase + '*/anat/*Stroke_mask.nii.gz', recursive=False)
    pathAnno = glob.glob(pathBase + '*/anat/*Anno.nii.gz', recursive=False)
    pathSIGMA = glob.glob(pathBase + '*/anat/*SIGMA.nii.gz', recursive=False)
    bsplineMatrix =  glob.glob(pathBase + '*/anat/*MatrixBspline.nii', recursive=False)
    return pathT2,pathStroke_mask,pathAnno,pathSIGMA,bsplineMatrix


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Registration of SIGMA Brain Atlas to rsfMRI')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputVolume', help='Path to rsfMRI data after preprocessing', required=True)
    parser.add_argument('-d', '--dtiasRef', action='store_true', help='use DTI as reference if data quality is low')
    parser.add_argument('--atlas-mask-t2', action='store_true',
                        help='mask the T2 BET with the registered atlas annotation before registration')
    parser.add_argument('-r', '--referenceDay', help='Reference Stroke mask', nargs='?', type=str,
                        default=None)
    parser.add_argument('-s', '--splitAnno', help='Split annotations atlas', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('-f', '--splitAnno_rsfMRI', help='Split annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')
    parser.add_argument('-a', '--anno_rsfMRI', help='Annotations atlas for rsfMRI', nargs='?', type=str,
                        default=os.path.abspath(os.path.join(os.getcwd(), os.pardir,os.pardir))+'/lib/sigma/SIGMA_InVivo_Anatomical_Brain_Atlas.nii.gz')



    args = parser.parse_args()



    stroke_mask = None
    inputVolume = None
    refStrokePath = None
    splitAnno = None
    splitAnno_rsfMRI = None
    anno_rsfMRI = None

    if args.inputVolume is not None:
        inputVolume = args.inputVolume
    if not os.path.exists(inputVolume):
        sys.exit("Error: '%s' is not an existing directory." % (inputVolume,))

    outfile = os.path.join(os.path.dirname(inputVolume))
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    setup_logging(outfile)


    # find related  data
    pathT2, pathStroke_mask, pathAnno, pathTemplate, bsplineMatrix = find_relatedData(os.path.dirname(outfile))
    if len(pathT2) == 0:
        T2data = []
        sys.exit("Error: %s' has no reference T2 template." % (os.path.basename(inputVolume),))
    else:
        T2data = require_single_match(pathT2, "reference T2 template")

    if len(pathStroke_mask) == 0:
        pathStroke_mask = []
        LOGGER.info("Notice: '%s' has no defined reference (stroke) mask - will proceed without.", os.path.basename(inputVolume))
    else:
        stroke_mask = require_single_match(pathStroke_mask, "stroke mask")

    if len(pathAnno) == 0:
        pathAnno = []
        sys.exit("Error: %s' has no reference annotations." % (os.path.basename(inputVolume),))
    else:
        brain_anno = require_single_match(pathAnno, "reference annotation")

    if len(pathTemplate) == 0:
        pathTemplate = []
        sys.exit("Error: %s' has no reference template." % (os.path.basename(inputVolume),))
    else:
        brain_template = require_single_match(pathTemplate, "reference template")

    if len(bsplineMatrix) == 0:
        bsplineMatrix = []
        sys.exit("Error: %s' has no bspline Matrix." % (os.path.basename(inputVolume),))
    else:
        bsplineMatrix = require_single_match(bsplineMatrix, "B-spline matrix")


    # find reference stroke mask
    refStroke_mask = None
    if args.referenceDay is not None:
        referenceDay = args.referenceDay
        refStrokePath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(outfile))), referenceDay)

        if not os.path.exists(refStrokePath):
            sys.exit("Error: '%s' is not an existing directory." % (refStrokePath,))
        refStroke_mask = find_RefStroke(refStrokePath, inputVolume)
        if len(refStroke_mask) == 0:
            refStroke_mask = []
            LOGGER.info("Notice: '%s' has no defined reference (stroke) mask - will proceed without.", os.path.basename(inputVolume))
        else:
            refStroke_mask = require_single_match(refStroke_mask, "reference stroke mask")

    if args.splitAnno is not None:
        splitAnno = args.splitAnno
    if not os.path.exists(splitAnno):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno,))

    if args.splitAnno_rsfMRI is not None:
        splitAnno_rsfMRI = args.splitAnno_rsfMRI
    if not os.path.exists(splitAnno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (splitAnno_rsfMRI,))

    if args.anno_rsfMRI is not None:
        anno_rsfMRI = args.anno_rsfMRI
    if not os.path.exists(anno_rsfMRI):
        sys.exit("Error: '%s' is not an existing directory." % (anno_rsfMRI,))

    output = regSIG2rsfMRI(inputVolume, T2data, brain_template, brain_anno, splitAnno, splitAnno_rsfMRI,
                           anno_rsfMRI, bsplineMatrix, args.dtiasRef, outfile, args.atlas_mask_t2)
    sys.stdout = sys.__stdout__

    current_dir = os.path.dirname(inputVolume)
    search_string = os.path.join(current_dir, "*EPI.nii.gz")
    currentFile = glob.glob(search_string)

    search_string = os.path.join(current_dir, "*.nii*")
    created_imgs = glob.glob(search_string, recursive=True)

    os.chdir(os.path.dirname(os.getcwd()))
    for idx, img in enumerate(created_imgs):
        if img == None:
            continue
        #os.system('python adjust_orientation.py -i '+ str(img) + ' -t ' + currentFile[0])

    LOGGER.info("Registration done")
