# Atlas-based Imaging Data Analysis Pipeline for Functional and Structural MRI Data

**AIDAmri v3.0**  
Aref Kalantari, Leon Scharwächter, Niklas Pallast, Michael Diedenhofen, Victor Vera Frazão, Marc Schneider, Markus Aswendt  
**Status:** May 2026  
Department of Neurology, University Hospital Frankfurt, Germany

## Contents

- [Introduction](#introduction)
  - [Atlas-based analysis](#atlas-based-analysis)
  - [Modular structure](#modular-structure)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Docker usage](#docker-usage)
  - [General overview](#general-overview)
  - [Creating image](#creating-image)
  - [Running container and mounting data](#running-container-and-mounting-data)
  - [Built-in installation legacy](#built-in-installation-legacy)
- [Functions](#functions)
- [Batch processing](#batch-processing)
- [Processing single files step-by-step](#processing-single-files-step-by-step)
  - [Convert raw data](#convert-raw-data)
  - [Processing of T2w and T2mapping data](#processing-of-t2w-and-t2mapping-data)
  - [Processing of ROI stroke mask data](#processing-of-roi-stroke-mask-data)
  - [Processing of DTI data](#processing-of-dti-data)
  - [Processing of fMRI data](#processing-of-fmri-data)
  - [Peri-infarct ROI analysis](#peri-infarct-roi-analysis)

## Introduction

The Atlas-based Processing Pipeline for functional and structural MRI data (**AIDAmri**) was developed for automated processing of mouse brain MRI. AIDAmri works with T2-weighted MRI (`anat`), diffusion weighted MRI or diffusion tensor imaging (`dwi`), resting-state functional MRI (`func`) and T2 map values (`t2map`).

### Atlas-based analysis

The Allen Mouse Brain Reference Atlas (ARA, CCF v3) is registered on each of these MRI data sets and is used to analyse regions of interest. Furthermore, the regions of the ARA are used as seed points for the connectivity and activity matrices. User-defined ROIs and masks can be generated separately and used for analysis, for example stroke lesion masks and peri-infarct regions.

AIDAmri comes with different atlas and template versions, which are necessary for the registration to work:

1. `annotation_50_changed_anno`: original ARA CCF v3 labels, with regions whose grey values are greater than 100,000 changed to new values starting with 2000.
2. `anno_volume_2000_rsfMRI`: atlas from item 1 with reduced number of atlas regions by selective region fusion, 96 regions in total, split between hemispheres, with the right side +2000.
3. Same as item 2 but not split.
4. Same as item 1 but split.
5. Original ARA template with 50 µm isotropic resolution.
6. Custom-made MRI template.

For the complete list of atlas labels, see:

```text
annoVolume+2000_rsfMRI.nii.txt
```

![Figure 1: Atlases included in the /lib folder.](images/figure-1-atlases.png)

*Figure 1: Atlases included in the `/lib` folder.*

### Modular structure

The AIDAmri processing pipeline consists of several modular units designed to process structural and functional MRI data. It is possible to apply an ROI, such as a stroke lesion mask, during the processing steps.

#### Data conversion

Converts raw MRI data into the BIDS structure and NIfTI format, including all necessary header information.

#### Pre-processing

- **Image re-orientation:** aligns images to a standard orientation.
- **Bias-field correction:** corrects inhomogeneities in the MRI data using the MICO algorithm.
- **Brain extraction:** removes non-brain tissues from the images.
- **Region-of-interest segmentation:** allows the user to define specific regions for detailed analysis, such as stroke lesions.

#### Registration

- Utilizes affine and non-linear transformations to align MRI data with the Allen Brain Reference Atlas (ARA).
- Applies transformations to ensure that all imaging data are mapped accurately to a common coordinate system for further analysis.

#### DTI and rs-fMRI processing

##### DTI

- **Motion correction:** motion artifacts are corrected on a slice-by-slice basis using FSL’s MCFLIRT tool, tailored to handle rapid breathing-induced movements typical in small animals.
- **Brain extraction:** non-brain tissues are removed using a binary mask generated during brain extraction.
- **Data reconstruction:** diffusion data are reconstructed using the electrostatically optimized protocol of Jones30 with 30 gradient directions.
- **Tractography:** performed using deterministic streamline propagation, starting from random voxel positions, with fiber tracking parameters optimized for true and false fiber generation.
- **Connectivity matrix creation:** connectivity matrices are generated, representing the strength of connections between ARA-defined brain regions.

##### rs-fMRI

- **Physiological recording correction:** respiratory artifacts are identified and corrected based on recorded breathing signals.
- **Motion correction:** similar to DTI, rs-fMRI data undergo slice-wise motion correction.
- **Smoothing and filtering:** spatial smoothing and high-pass filtering, with cut-off at 0.01 Hz, are applied to reduce noise and enhance the signal.
- **Time-series extraction:** time-series data for each ARA-defined region are calculated by averaging voxel intensities over time within each region.
- **Functional connectivity analysis:** correlation of BOLD signals between brain regions is computed to assess functional connectivity.
- **Output:** produces connectivity matrices that represent structural and functional connections within the brain, which can be used for further analysis, such as graph theory applications.

## Installation

AIDAmri is distributed as a Docker image. If you prefer the built-in installation instead, see the [legacy installation](#built-in-installation-legacy) section. We highly recommend the usage of the Docker image.

### Prerequisites

The following are required to launch your AIDAmri instance:

- Docker engine
  - Getting started tutorial for Docker
- Windows only: Bash subsystem, for example Git Bash
- At least 10.6 GB free disk memory

We advise you to get comfortable with shell or command-line interface usage.

Download or clone the repository:

```text
git clone https://github.com/aswendtlab/AIDAmri.git
```

### Docker usage

This guide introduces the usage of Docker-based containers of the AIDAmri tools for Unix/Linux-based systems, including Linux and macOS. The commands shown are written for such systems, and you may copy the commands into your shell including backslashes, as they indicate line breaks.

Windows users may use a subsystem like Git Bash to use the software. You may also use the Docker Desktop application to get access to the Docker image.

### General overview

The AIDAmri pipeline is containerized and structured as depicted in Figure 2. The `Dockerfile` located in the repository provides the installation routine for every required dependency. The `docker build` command constructs the image, meaning the installed software on your system. The `docker run` command creates a runnable instance of this image, called a container.

The container provides a command-line interface. It can be accessed by directly attaching to an interactive interface that lets you input AIDAmri commands within the isolated file system of the container, or via the `docker exec` command from your host shell.

![Figure 2: Docker architecture draft.](images/figure-2-docker-architecture.png)

*Figure 2: Docker architecture draft. The blue boxes within the Dockerfile box depict the main content layers. The boxes preceded with `$` are command-line codes.*

The container follows a basic Ubuntu 18.04 system, meaning that the root directory is called `/`. To share a volume between the host system and the container, bind mounts are commonly used. See [Running container and mounting data](#running-container-and-mounting-data) for further information.

When referring to the mounted volume while in the container, use the path given at mounting, for example `/<MOUNTED DIRECTORY>`. This path will likely be different on your host system.

![Figure 3: General structure of container file structure.](images/figure-3-container-filesystem.png)

*Figure 3: General structure of container file structure. The paths are constructed from left to right, for example the absolute path of the `bin` folder in `aida` would be `/aida/bin`. Keep in mind that `/` is its own directory. The `<MOUNTED DIRECTORY>` and `<DATA>` directories are the same but can be named differently, depending on how it was named when mounted.*

### Creating image

Before you can build the Docker image, you need to open a terminal. In the terminal, change into the AIDAmri repository folder that you previously cloned from GitHub:

```text
cd PATH/TO/AIDAmri
```

Replace `PATH/TO/AIDAmri` with the actual path to your local AIDAmri folder.

Check the folder contents with:

```text
ls
```

A file named `Dockerfile`, as well as `fslinstaller_mod.py`, a `bin/` folder and a `lib/` folder should be located in this directory. Then launch the Docker daemon to build the image:

```text
docker build -t aidamri:latest -f Dockerfile .
```

The created image is a template for running containers and instantiating the pipeline. Be aware that the **period** at the end is part of the command and refers to the corresponding directory.

The `-t` flag sets the name and tag of the image. In this example, it is called `aidamri` and tagged `latest`. You may change the name and tag, but remember to change them accordingly in later steps that invoke the image. The `-f` flag refers to the `Dockerfile` in the current directory.

You only need to build the Docker image once during the initial installation.  
After updating the GitHub repository, for example with `git pull`, you should rebuild the image so that the new code is included.

Docker usually uses its build cache during this process. This means that unchanged parts of the image are reused, and only the layers affected by the update are rebuilt. Therefore, rebuilding the image after a code update is usually faster than the initial build.

Keep in mind that existing containers are not updated automatically. To use the updated image, stop and remove the old container, then start a new container from the rebuilt image.
### Running container and mounting data

Before starting the container, check where your MRI data is stored on your computer. You will need the **absolute path** to this folder.
To start an AIDAmri container and make your data available inside it, run:

```text
docker run -dit \
  --name aidamri_container \
  --mount type=bind,source=PATH/TO/DATA,target=/aida/DATA \
  aidamri:latest
```

The command performs the following function:
`docker run` starts a new container from the AIDAmri Docker image.
`-dit` starts the container in the background while keeping it interactive. This means the container keeps running, and you can enter it later.

`--name aidamri_container` gives the container a name. The container name is independent from the image name. In this example, the image is called `aidamri:latest`, while the container is called `aidamri_container`. If you wish to use more than one running container instance, for example to process multiple datasets simultaneously, each container needs a different name.

`--mount type=bind,source=PATH/TO/DATA,target=/aida/DATA` connects a folder from your computer to a folder inside the container. This is called a bind mount. It allows AIDAmri to access and process your MRI data without copying it into the container.

`source=PATH/TO/DATA` is the data folder where your data that you want to process is stored on your computer. Always use an absolute path here! 
`target=/aida/DATA` is the location where the same folder will appear inside the container. So you do not need to use the absolute path inside the container, but can refer to the mounted data with `/aida/DATA`.

`aidamri:latest` is the Docker image that is used to create the container.

After running the command, Docker will print a container ID. You can check whether the container is running with:

```text
docker aidamri_container ls
```

To enter the running container, use:

```text
docker attach aidamri_container
```

Windows users can leave a running container without stopping it by pressing `CTRL+P` and `CTRL+Q` consecutively. Alternatively, typing `exit` will stop the container.

Use the following to re-run the container:

```text
docker start aidamri_container
```

Type `stop` instead of `start` to stop an already running container.

After successfully attaching to the container, your terminal should look similar to this:

```text
root@<SOME NUMBERS AND CHARACTERS>:/aida#
```

The number-and-letter combination is the first part of the container ID.  
You are now inside the running AIDAmri container.

From here, you can use the AIDAmri commands described in the usage sections.  
As a first test, change into the `bin/` folder and open the help page of the batch processing script:

```text
cd bin/
python batchProg.py -h
```

If the help page is displayed, the container is working and AIDAmri can be used.

Alternatively, you can use `docker exec` to run a command inside the container without entering the container shell.

For example, to open the help page of `batchProg.py`, run:

```text
docker exec -w /aida/bin aidamri-container \
  python batchProg.py -h
```

Here, `-w /aida/bin` sets the working directory inside the container to `/aida/bin`.

## Functions

List of functions and script groups:

- `conv2Nifti_auto.py`: batch conversion of raw Bruker ParaVision data to BIDS-like NIfTI output folders.
- `PV2NIfTiConverter/`: ParaVision-to-NIfTI conversion scripts used by `conv2Nifti_auto.py`, including DTI sidecar generation for `.bval` and `.bvec` files.
- `batchProc.py`: batch processing entry point for converted datasets. It runs the selected preprocessing, registration and processing steps for the requested data types and sessions.
- `2.1_T2PreProcessing/`: T2w preprocessing, including reorientation, smoothing, bias-field correction, brain extraction and atlas registration.
- `2.2_DTIPreProcessing/`: DTI preprocessing, including b0 averaging, motion correction, smoothing, bias-field correction, brain extraction and atlas registration.
- `2.3_fMRIPreProcessing/`: rs-fMRI preprocessing, including slice-time correction, motion correction, smoothing and atlas registration.
- `3.1_T2Processing/`: T2w analysis scripts for stroke mask statistics, incidence maps, incidence sizes and SNR calculations.
- `3.2_DTIConnectivity/`: DTI reconstruction and whole-brain fiber tracking with DSI Studio, including connectivity matrix generation and matrix plotting.
- `3.2.1_DTIdata_extract/`: extraction of region-wise DTI measures such as FA, AD, RD and MD from registered atlas regions, including iterative batch helpers.
- `3.3_fMRIActivity/`: rs-fMRI activity and connectivity analysis, including seed ROI creation, regression tables, mean time-series extraction, correlation matrices and matrix plotting.
- `4.1_T2mapPreProcessing/`: T2 map preprocessing, atlas registration and extraction of region-wise T2 map values.
- `5.1_ROI_analysis/`: ROI-based analyses for user-defined regions, for example peri-infarct regions around stroke lesions, including mask dilation, transform application, seed ROI creation and ROI inspection.
- `helper_tools/`: additional utilities for data preparation and quality control, including naming cleanup, batch reorientation, fieldmap JSON updates, stroke mask distribution, source-data plotting and atlas region size summaries.
- `adjustbvecRep.py`: helper for adjusting repeated b-vector files in DTI datasets.

All program examples are listed only with the mandatory input parameters. For more details or help, call:

```text
python <command> -h
```

The command-line examples use the identifier `testData<No.>.nii.gz` and can be identically applied to other data. The test data are freely available at DOI `10.12751/g-node.70e11f`.

After a successful download, you can choose either to process single files manually or automate the processing for the whole dataset. In both cases, processing includes file conversion from the raw Bruker format into NIfTI format, several preprocessing steps and registration with the Allen Brain Reference Atlas. The functions in `/bin` are named according to the MRI sequence to be processed, such as `T2`, `DTI` and `fMRI`.

> [!IMPORTANT]
> Process your T2 data first so that preprocessing of DTI, fMRI and T2 map data works correctly.

## Batch processing

AIDAmri provides functions for data conversion and batch processing. Complete processing requires two scripts:

1. `conv2Nifti_auto.py` Converts all files to NIfTI format and stores them in the new folder (proc folder).
2. `batchProc.py` applies preprocessing steps and registration with the atlas.

> [!NOTE]
> If multiple reconstructions exist, conversion will only use the first folder correctly. The test data set is already converted into NIfTI format (see nifti folder), so only the second script needs to be applied.

In general, raw Bruker data must be in the following structure for the first script to work:

```text
projectfolder/days/subjects/
```

To convert the whole project folder into NIfTI format, open the terminal and change the directory to the `/bin` folder of the AIDAmri installation:

```text
cd <path to AIDAmri>/bin
```

Start Bruker to NIfTI conversion:

```text
python conv2Nifti_auto.py -i /path/to/raw_dataset -o /path/to/output
```

This script automatically finds all raw Bruker datasets saved within the input path. You can specify the output directory via the `-o` flag. Without the `-o` flag, the output is saved next to the input directory:

```text
/path/to/proc_data
```


After successful Bruker to NIfTI conversion, the second script can be applied to the new project folder `proc_data`. The data need to be ordered in BIDS format like the output of `conv2Nifti_auto.py`:

```text
projectfolder/sub-/ses-/datatype
```

> [!WARNING]
> Batch processing may slow down the system depending on the CPU load-out. Use the `-cpu` or `-e cpu` flag to specify CPU usage. Run `python batchProg.py -h` for more information.

Example:

```text
python batchProc.py -i /path/to/proc_data \
  -t anat dwi func t2map -s Baseline P3 P12 --t2-bias-method mico
```

This script runs every necessary script for preprocessing, registration and processing steps. You can specify which data types (`-t`) e.g anatomical or dwi and sessions (`-s`) to compute. You can also specify which bias method should be used on the data.

You do not need to specify data types and sessions or any other argument except the input argument. If no `-t` and no `-s` flag are given, every data type and session of every subject will be processed.

> [!IMPORTANT]
> The scripts executed by `batchProc.py` are related to each other. Therefore, `anat` processing always needs to be performed before `dwi`.

Depending on the size of your project, this process may take a while. After finishing, the project folder is ready for network graph analysis, for example using AIDAconnect.

Please have a look into the help page of `batchProc.py` for more information on all the options which you can select:

```text
python batchProc.py -h
```

## Processing single files step-by-step
Every script has arguments that can be specified when calling the script. For more information on the arguments, every script has a help page which can be accessed by calling the script with the `-h` flag. The following sections provide examples of how to process single files step-by-step.

### Convert raw data

Convert Bruker raw data to NIfTI files by specifying the folder containing all raw folders of each scan.

> [!NOTE]
> If multiple reconstructions exist, conversion will only use the first folder correctly. A file with exactly the same name is created in the given input folder. It contains all sorted NIfTI files. The raw data should have the same orientation as the example dataset.

```text
python pv_conv2Nifti.py -i /aida/DATA/raw_data_folder
```

Move the newly generated file to a new project folder if you want to separate raw Bruker files from processed NIfTI files. We recommend the following folder structure, especially if you want to use AIDAconnect for graph analysis:

```text
projectfolder/days/groups/subjects/data/
```

### Processing of anatomical data
> [!WARNING]
> Before you process any data please visually inspect the NIFTIs to check for correct orientation and quality. For more information please see the "Data Format and Orientation Requirements" section in the [README](README.md).

Apply bias field correction and brain extraction to the anatomical data set. The automatically attached endings of the processed filenames indicate which steps have been performed.


```text
python preProcessing_T2.py -i /aida/DATA/PATH/TO/anat/testData.5.1.nii.gz
```
You can specify the bias correction method with the `-b` flag.  
By default, AIDAmri uses `MICO`, which is well suited for small animal MRI data. Alternatively, you can use `ANTS`. The bias-corrected output file is saved with the ending `...Bias.nii.gz`.
For brain extraction, several parameters can be adjusted. For example, the `-f` flag defines the fractional intensity threshold used by FSL's BET method. You can also specify parameters such as the brain radius or the horizontal gradient if you want to make the extraction stricter in the anterior or posterior direction.
The default brain extraction method is FSL's `BET`, which was originally developed for human brain MRI but is used here in a modified way. As an alternative, AIDAmri also supports `bet4animal`, which is specifically designed for small animal brains, such as mouse or rat brains. `bet4animal` is often easier to use because parameters such as the fractional intensity threshold do not need to be selected manually.
However, in practice, neither method is always superior. In some cases, the default `BET` method gives better results, while in other cases `bet4animal` works better. Therefore, we recommend testing both methods and visually checking the resulting brain extraction.
The brain-extracted output file is saved with the ending `...BET.nii.gz`.

### Registration of anatomical data
The next step includes registration of the Allen Brain Reference Atlas with the brain-extracted T2 dataset.  

There is an option to segment an additional region of interest, such as the stroke lesion. You can segment the region using the brain-extracted dataset as reference, ending with `...BET.nii.gz`. We recommend conducting this step with ITK-SNAP. The saved file should end with: `...Stroke_mask.nii.gz`

Run registration:

```text
python registration_T2.py -i /aida/DATA/PATH/TO/anat/testDataBiasBet.nii
```

For an input file called `<input>.nii.gz`, `registration_T2.py` creates the following output files in the same `anat` folder:

```text
<input>_TemplateAff.nii.gz
<input>_Template.nii.gz
<input>_TemplateAllen.nii.gz
<input>MatrixAff.txt
<input>MatrixInv.txt
<input>MatrixBspline.nii
<input>_Anno.nii.gz
<input>_AnnoSplit.nii.gz
<input>_Anno_parental.nii.gz
<input>_AnnoSplit_parental.nii.gz
```

- `<input>_TemplateAff.nii.gz`: MRI template (NP_template_sc0.nii.gz) after affine registration to the T2 image.
- `<input>_Template.nii.gz`: MRI template (NP_template_sc0.nii.gz) after non-linear registration to the T2 image.
- `<input>_TemplateAllen.nii.gz`: Allen Brain Reference Template registered to the T2 image.
- `<input>MatrixAff.txt`: affine transformation matrix from template space to T2 space.
- `<input>MatrixInv.txt`: inverse affine transformation matrix from T2 space to Allen template space.
- `<input>MatrixBspline.nii`: non-linear B-spline transformation field.
- `<input>_Anno.nii.gz`: registered detailed Allen atlas annotation in T2 space.
- `<input>_AnnoSplit.nii.gz`: registered detailed Allen atlas annotation in T2 space with separated left and right hemispheres.
- `<input>_Anno_parental.nii.gz`: registered parental atlas annotation in T2 space with larger brain regions.
- `<input>_AnnoSplit_parental.nii.gz`: registered parental atlas annotation in T2 space with separated left and right hemispheres.

Check the registration result visually, for example by overlaying the brain-extracted image with the registered atlas annotation file.
If the registration result is not satisfactory, try to improve the brain extraction first. For example, you can adjust the brain extraction parameters or manually correct the generated brain mask using tools such as ImageJ. After improving the mask, run the registration again.

The script also creates an `IncidenceData` subfolder:

```text
IncidenceData/<input>_IncidenceData.nii.gz
IncidenceData/<input>_IncidenceData_Lesion_mask.nii.gz
```

- `<input>_IncidenceData.nii.gz`: T2 image registered into Allen template space, used for incidence-map processing.
- `<input>_IncidenceData_Lesion_mask.nii.gz`: stroke or lesion mask registered into Allen template space. This file is only created if a matching `*Stroke_mask.nii.gz` file exists in the same folder.

### Processing of anatomical data
If you previously defined a region of interest, such as a stroke lesion, you can calculate the ROI size and determine which atlas regions overlap with it. In this step, the segmented ROI file, for example `...Stroke_mask.nii.gz`, is overlaid with the registered atlas annotation.

Two atlas variants can be evaluated:

- `getIncidenceSize.py` uses the regular left/right-separated ARA atlas `ARA_annotationR+2000.nii.gz` and the subject-space annotation file `*_AnnoSplit.nii.gz`.
- `getIncidenceSize_par.py` uses the parental atlas `annoVolume+2000_rsfMRI.nii.gz` and the subject-space annotation file `*_AnnoSplit_parental.nii.gz`.

Both scripts expect the corresponding `.../anat` folder as input. The folder must contain exactly one stroke mask, one BET image, one matching annotation file and one `IncidenceData_Lesion_mask.nii.gz` file.

```text
python getIncidenceSize.py -i .../testData/anat
python getIncidenceSize_par.py -i .../testData/anat
```

The non-parental affected-region results from `getIncidenceSize.py` are stored in `.../anat/affected_Regions`:

```text
*affectedRegions.csv
*affectedRegions.nii.gz
*labelCount.mat
```

The labelled non-parental incidence lesion mask is stored in `.../anat/IncidenceData`:

```text
*IncidenceData_Anno_lesion_mask.nii.gz
```

The parental affected-region results from `getIncidenceSize_par.py` are stored in `.../anat/affected_Regions`:

```text
*affectedRegions_Parental.csv
*affectedRegions_Parental.nii.gz
*labelCount_par.mat
```

The labelled parental incidence lesion mask is stored in `.../anat/IncidenceData`:

```text
*IncidenceData_Anno_parental_lesion_mask.nii.gz
```

### Processing of ROI stroke mask data

From masks drawn on the T2-weighted images, it is possible to determine both the incidence map and the size of affected regions. For example, if a `day1` folder contains multiple `Mouse 1` to `Mouse 15` folders and the processed T2 data are in those folders, the command would be:

```text
python getIncidenceMap.py -i .../day1 -s "Mouse*"
```

It is also possible to determine the region size as voxels and volume in mm³:

```text
python getRegionSize_par.py -i .../T2w
```

### Processing of DTI data

The DTI processing procedure includes dimension reduction, bias correction, threshold application and subsequent brain extraction. The endings on the filenames indicate which steps have been performed.

```text
python preProcessing_DTI.py -i .../DTI/testData.7.1.nii.gz
```

The next step includes registration of the Allen Brain Reference Atlas with the brain-extracted DTI dataset. For processing a reference stroke mask, two options are available:

1. Registration of a reference mask that is related to another dataset or day, for example to always use the same mask. Append `-r <filename of ref>`.
2. Otherwise, the algorithm automatically uses the corresponding reference mask from the respective subject folder. If no mask is defined, the registration proceeds without a mask.

```text
python registration_DTI.py -i .../DTI/testDataSmoothMicoBet.nii.gz
```

Connectivity is finally calculated using DSI Studio. All connectivity matrices are based on the reference atlas.

```text
python dsi_main.py -i .../DTI/testData.7.1.nii.gz
```

The connectivity matrices of the parental ARA, the original ARA and the related ROI are stored in the folder `.../DTI/connectivity` as `.txt` and `.mat` files. DSI Studio differentiates between matrices that count how many fibers pass through and end in each region.

The adjacency matrices can be visualized using the related plot function:

```text
python plotDTI_mat.py -i .../testData/DTI/connectivity/testData*.connectivity.mat
```

The folder `.../DTI/DSI_studio` also contains diffusion value maps, for example FA maps, registered with the atlas. This data can be extracted and saved as `.txt` with the region name and corresponding FA, RD, MD and AD values using:

```text
python DTIdata_extract.py image_file roi_file
```

Use this command in the `3.2.1 DTIdata_extract` folder. To iteratively process all subjects, use the `iterativeRun.py` function.

### Processing of fMRI data

The fMRI processing is roughly comparable to preprocessing of DTI datasets.

> [!IMPORTANT]
> Brain extraction should be of good quality and must be manually checked or corrected by adapting the given parameters.

```text
preProcessing_fMRI.py -i .../fMRI/testData.6.1.nii.gz
```

The next step includes registration of the Allen Brain Reference Atlas with the brain-extracted fMRI dataset. The result is a variety of files. An impression of the registration can be obtained by superimposing the brain-extracted file with the annotations of the Allen Brain, ending with `...Anno.nii.gz`.

```text
python registration_fMRI.py -i .../testData/fMRI/testSmoothBet.nii
```

If physiological data are not available, the step will be conducted without the included regression. All activity matrices are based on the reference atlas.

```text
python process_fMRI -i .../fMRI/testData.6.1.nii.gz
```

The activity matrices of the parental Atlas and original Atlas are stored in the folder `.../fMRI/regr` as `.txt` and `.mat` files with the prefixes `MasksTCs.` and `MasksTCsSplit.`.

The related adjacency matrices can be visualized using the related plot function:

```text
python plotfMRI_mat.py -i .../testData/fMRI/regr/MasksTCsSplit*.mat
```

### Peri-infarct ROI analysis

You can create custom peri-infarct masks to further analyze stroke-related regions. Go to the folder:

```text
bin/4.1 ROI analysis
```

Open `proc_tools.py` with an editor that can open Python files. Adjust all directories, paths and further specifications as described in the script.

To decide which regions to include in the peri-infarct region, modify:

```text
cortex_labels_1.txt
cortex_labels_2.txt
```

For the full list of atlas labels, see:

```text
../lib/annoVolume+2000_rsfMRI.nii.txt
```

Proceed with the scripts in order from 1 to 4.

The first script creates peri-infarct masks for all time points:

```text
python 01_dilate_mask_process.py
```

The second script aligns the peri-infarct masks in the rs-fMRI and DTI space:

```text
python 02_apply_xfm_process.py
```

The result of the third script depends on the imaging type:

- For rs-fMRI, a MATLAB file is created which contains two text files:
  1. For each region, one column with the averaged rs-fMRI time series.
  2. The atlas label names.
- For DTI, a modified atlas labels file is created which includes individually shaped peri-infarct brain regions. These newly generated regions replace the original regions in the file.

```text
python 03_create_seed_rois_process.py
```

The fourth script is not mandatory, but is a helper tool to compare the number of voxels included in the peri-infarct region for each subject.

```text
python 04_examine_rois.py
```

> [!WARNING]
> The scripts for peri-infarct ROI analysis are provided for analysis of time point 7 only, for example 7 days post-stroke. For other time points, manual modifications are necessary.
