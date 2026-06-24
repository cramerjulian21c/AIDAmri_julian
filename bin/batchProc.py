"""
Created on 18/11/2020

@author: Marc Schneider
AG Neuroimaging and Neuroengineering of Experimental Stroke
Department of Neurology, University Hospital Cologne

This script runs every needed script for all (pre-)processing and registration
steps. The data needs to be ordered like after Bruker2NIfTI conversion:
project_folder/days/groups/subjects/.
For the script to work, it needs to be placed within the /bin folder of AIDAmri.

Example:
python batchProc.py -i /Volumes/Desktop/MRI/proc_data -t anat dwi func t2map
"""

import os
import fnmatch
from pathlib import Path
import concurrent.futures
import subprocess
from tqdm import tqdm
import multiprocessing
import logging
import shlex
import time

FATAL_LIP_HEADER_EXIT_CODE = 86


def findData(projectPath, sessions, data_types):
    if not data_types:
        data_types = ["anat", "dwi", "func", "t2map"]
    # This function screens all existing paths. Within these paths, this function collects all subject
    # folders, which are all folders that are not named 'Physio'.
    full_path_list = sorted(os.listdir(projectPath))
    all_wanted_paths, anat_files, dwi_files, func_files, t2map_files = [], [], [], [], []

    # collect ses paths
    for path in full_path_list:
        if path.startswith("sub-"):
            sub_root = os.path.join(projectPath, path)
            wanted_paths = sorted(os.listdir(sub_root))
            wanted_paths = [
                os.path.join(sub_root, wp)
                for wp in wanted_paths
                if wp.startswith("ses-")
            ]
            all_wanted_paths.extend(wanted_paths)

    # filter sessions (exact match on path components)
    if sessions:
        wanted = {f"ses-{s}" for s in sessions}
        matching_paths = []
        for p in all_wanted_paths:
            parts = os.path.normpath(p).split(os.sep)
            if any(part in wanted for part in parts):
                matching_paths.append(p)
        all_wanted_paths = matching_paths

    # collect datatype folders
    for path in all_wanted_paths:
        for sub_dir in sorted(os.listdir(path)):
            if sub_dir == "anat" and "anat" in data_types:
                anat_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "dwi" and "dwi" in data_types:
                dwi_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "func" and "func" in data_types:
                func_files.append(os.path.join(path, sub_dir))

            elif sub_dir == "t2map" and "t2map" in data_types:
                t2map_files.append(os.path.join(path, sub_dir))

    return {"anat": anat_files, "dwi": dwi_files, "func": func_files, "t2map": t2map_files}

def _get_arg_after(flags, argv):
    for f in flags:
        if f in argv:
            i = argv.index(f)
            if i + 1 < len(argv):
                return argv[i + 1]
    return None

def _log_base_from_input(input_path: str) -> str:
    # If input is a dir -> log in that dir
    # If input is a file -> log in its parent dir
    return input_path if os.path.isdir(input_path) else os.path.dirname(input_path)

def _quote(value) -> str:
    '''
    Safely quote a value for shell command usage, handling spaces and special characters.
    '''
    return shlex.quote(str(value))

def run_subprocess(command, datatype, step, anat_process=False):
    timeout = 5400 #timeout (sec) for subprocess
    command_args = shlex.split(command)

    inp = _get_arg_after(["-i", "--input", "--input-file"], command_args)
    if inp is None:
        inp = next((a for a in reversed(command_args)
                    if a.endswith(".nii") or a.endswith(".nii.gz")), command_args[-1])

    base = _log_base_from_input(inp)

    # default location
    log_file = os.path.join(base, f"{step}.log")

    # special case: anat/process wants different filenames
    if datatype == "anat" and step == "process":
        log_name = f"{step}.log" if anat_process else f"{step}_par.log"
        log_file = os.path.join(base, log_name)

    #Determine sub / ses
    normalized_path = os.path.normpath(inp)
    directories = normalized_path.split(os.path.sep)
    sub = next((d for d in directories if d.startswith("sub-")), "sub-UNKNOWN")
    ses = next((d for d in directories if d.startswith("ses-")), "ses-UNKNOWN")

    try:
        logging.info(f"Running command: {command}.\nCheck {log_file} for further information.")
        with open(log_file, 'w') as outfile:
            time.sleep(2) # make sure logging file is created before starting the subprocess
            child_env = os.environ.copy()
            # dsi_main.py can create its own process.log during interactive
            # runs. Disable that side log here because batchProc.py already
            # captures stdout/stderr into the step-specific batch log.
            if any(arg.endswith("dsi_main.py") for arg in command_args):
                child_env["AIDAMRI_DISABLE_PROCESS_LOG"] = "1"
            child_env["AIDAMRI_DISABLE_SPINNER"] = "1"
            result = subprocess.run(
                command_args,
                stdout=outfile,
                stderr=outfile,
                text=True,
                timeout=timeout,
                env=child_env,
            )
            if result.returncode != 0:
                if (
                    result.returncode == FATAL_LIP_HEADER_EXIT_CODE
                    and datatype in {"anat", "dwi"}
                    and step == "preprocess"
                ):
                    raise RuntimeError(
                        f"Fatal header check failure in {inp}. Expected LIP orientation."
                    )
                return sub,ses,datatype,step
            else:
                return 0
    except subprocess.TimeoutExpired:
        logging.error(f'Timeout expired for command: {command_args}')
        return sub,ses,datatype,step
    except Exception as e:
        logging.error(f'Error while executing the command: {command_args} Errorcode: {str(e)}')
        raise
    

def executeScripts(currentPath_wData, dataFormat, step, cfg, stc=False):
    # For every datatype (T2w, fMRI, DTI), go in all days/group/subjects folders
    # and execute the respective (pre-)processing/registration-scripts.
    # If a certain file does not exist, a note will be created in the errorList.
    # cwd should contain the path of the /bin folder (the user needs to navigate to the /bin folder before executing this script)
    #KEEP IN MIND DUE TO PARALLEL COMPUTING NO ERRORS IN THIS FUNCTION WILL BE PRINTED OUT => GREY ZONE
    errorList = [];
    message = '';
    cwd = str(Path(__file__).resolve().parent)
    currentPath_wData = Path(currentPath_wData)
    # currentPath_wData = projectfolder/sub/ses/dataFormat (e.g. anat, func, dwi)
    if os.path.isdir(currentPath_wData):
        if dataFormat == 'anat':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.1_T2PreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*T2w.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_T2.py -i {_quote(currentFile[0])}'

                    # Bias field correction for T2: none | mico | ants
                    if cfg.get("t2_bias_method"):
                        command += f' -b {cfg["t2_bias_method"]}'

                    command += f' --bet {cfg["t2_bet"]}'

                    # BET-Parameter
                    if cfg.get("t2_frac") is not None:
                        command += f' -f {cfg["t2_frac"]}'
                    if cfg.get("t2_radius") is not None:
                        command += f' -r {cfg["t2_radius"]}'
                    if cfg.get("t2_gradient") is not None:
                        command += f' -g {cfg["t2_gradient"]}'
                    if cfg.get("t2_center") is not None:
                        cx, cy, cz = cfg["t2_center"]
                        command += f' -c {cx} {cy} {cz}'

                    result = run_subprocess(command, dataFormat, step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *T2w.nii.gz in {str(currentPath_wData)}'
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)

            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.1_T2PreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile) > 0:
                    r1 = run_subprocess(f'python registration_T2.py -i {_quote(currentFile[0])}', dataFormat, step)
                    if r1 != 0:
                        errorList.append(r1)
                    r2 = run_subprocess(f'python t2_value_extraction.py -i {_quote(currentFile[0])}', dataFormat, step)
                    if r2 != 0:
                        errorList.append(r2)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}'
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)

            elif step == "process":
                has_stroke_mask = any(currentPath_wData.glob("**/*Stroke_mask.nii.gz"))
                if not has_stroke_mask:
                    message = f"No stroke mask found for {str(currentPath_wData)}, proceeding without mask."
                    logging.info(message)  #write in log-file
                    #print(message, flush=True)
                    return 0
                os.chdir(os.path.join(cwd, '3.1_T2Processing'))
                if cfg.get("t2_incidence_script") == "detailed":
                    command = f'python getIncidenceSize.py -i {_quote(currentPath_wData)}'
                    result = run_subprocess(command, dataFormat, step, anat_process=True)
                else:
                    command = f'python getIncidenceSize_par.py -i {_quote(currentPath_wData)}'
                    result = run_subprocess(command, dataFormat, step)

                if result != 0:
                    errorList.append(result)

                os.chdir(cwd)


        elif dataFormat == 'func':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.3_fMRIPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_fMRI.py -i {_quote(currentFile[0])}'
                    if cfg.get("func_bias_method") is not None:
                        command += f' -b {cfg["func_bias_method"]}'
                    command += f' --bet {cfg["func_bet"]}'
                    if cfg.get("func_frac") is not None:
                        command += f' -f {cfg["func_frac"]}'
                    if cfg.get("func_radius") is not None:
                        command += f' -r {cfg["func_radius"]}'
                    if cfg.get("func_gradient") is not None:
                        command += f' -g {cfg["func_gradient"]}'
                    if cfg.get("func_center") is not None:
                        cx, cy, cz = cfg["func_center"]
                        command += f' -c {cx} {cy} {cz}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *EPI.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.3_fMRIPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_rsfMRI.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = sorted(currentPath_wData.glob("*EPI.nii.gz"))
                if len(currentFile)>0:
                    os.chdir(os.path.join(cwd, '3.3_fMRIActivity'))
                    command = f'python process_fMRI.py -i {_quote(currentFile[0])} -stc {stc} --bet {cfg["func_bet"]}'
                    if cfg.get("func_frac") is not None:
                        command += f' --bet-frac {cfg["func_frac"]}'
                    if cfg.get("func_radius") is not None:
                        command += f' --bet-radius {cfg["func_radius"]}'
                    if cfg.get("func_gradient") is not None:
                        command += f' --bet-gradient {cfg["func_gradient"]}'
                    if cfg.get("func_center") is not None:
                        cx, cy, cz = cfg["func_center"]
                        command += f' -ctr {cx} {cy} {cz}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                    os.chdir(cwd)
                else:
                    message = f'Could not find *EPI.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
        elif dataFormat == 't2map':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*MEMS.nii.gz"))
                if len(currentFile)>0:
                    command = f'python preProcessing_T2MAP.py -i {_quote(currentFile[0])}'
                    if cfg.get("t2map_bias_method"):
                        command += f' -b {cfg["t2map_bias_method"]}'
                    command += f' --bet {cfg["t2map_bet"]}'
                    if cfg.get("t2map_frac") is not None:
                        command += f' -f {cfg["t2map_frac"]}'
                    if cfg.get("t2map_radius") is not None:
                        command += f' -r {cfg["t2map_radius"]}'
                    if cfg.get("t2map_gradient") is not None:
                        command += f' -g {cfg["t2map_gradient"]}'
                    if cfg.get("t2map_center") is not None:
                        cx, cy, cz = cfg["t2map_center"]
                        command += f' -c {cx} {cy} {cz}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *MEMS.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = sorted(
                    currentPath_wData.glob("*Bet.nii.gz"),
                    key=lambda path: path.stat().st_mtime,
                    reverse=True,
                )
                if len(currentFile)>0:
                    command = f'python registration_T2MAP.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {str(currentPath_wData)}';
                    print(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                os.chdir(os.path.join(cwd, '4.1_T2mapPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*T2w_MAP.nii.gz"))
                if len(currentFile)>0:
                    command = f'python t2map_data_extract.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *T2w_MAP.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
        elif dataFormat == 'dwi':
            if step == "preprocess":
                os.chdir(os.path.join(cwd, '2.2_DTIPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*dwi.nii.gz"))
                if len(currentFile) > 0:
                    command = f'python preProcessing_DTI.py -i {_quote(currentFile[0])}'

                    # DWI BET parameter (only append if set, otherwise use script defaults)
                    if cfg.get("dwi_frac") is not None:
                        command += f' -f {cfg["dwi_frac"]}'
                    if cfg.get("dwi_radius") is not None:
                        command += f' -r {cfg["dwi_radius"]}'
                    if cfg.get("dwi_gradient") is not None:
                        command += f' -g {cfg["dwi_gradient"]}'

                    # Bias field
                    # dwi_bias_method with choices ["mico",“ants”], default=None
                    if cfg.get("dwi_bias_method") is not None:
                        command += f' -b {cfg["dwi_bias_method"]}'

                    # Denoiser
                    if cfg.get("dwi_denoiser"):
                        command += f' --denoiser {cfg["dwi_denoiser"]}'

                    command += f' --bet {cfg["dwi_bet"]}'

                    if cfg.get("dwi_average_b0"):
                        command += ' --average-b0'

                    if cfg.get("dwi_skip_min_projection"):
                        command += ' --skip-min-projection'

                    result = run_subprocess(command, dataFormat, step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *dwi.nii.gz in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "registration":
                os.chdir(os.path.join(cwd, '2.2_DTIPreProcessing'))
                currentFile = sorted(currentPath_wData.glob("*Bet.nii.gz"))
                if len(currentFile)>0:
                    command = f'python registration_DTI.py -i {_quote(currentFile[0])}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find *Bet.nii.gz in {currentPath_wData}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
            elif step == "process":
                currentFile = sorted(currentPath_wData.glob("*dwi.nii.gz"))
                if cfg.get("dwi_denoiser") == "patch2self":
                    currentFile = sorted(currentPath_wData.glob("*Patch2SelfDenoised.nii.gz"))
                # Appends optional (fa0, nii_gz) flags to DTI main process if passed
                if len(currentFile)>0:
                    # Pull values from cfg (with defaults)
                    track_param = cfg.get("dsi_track_param", "default")
                    if isinstance(track_param, (list, tuple)):
                        track_param_args = ' '.join(_quote(item) for item in track_param)
                    else:
                        track_param_args = _quote(track_param)
                    recon_method = cfg.get("dsi_recon_method", "dti")
                    vivo = cfg.get("dsi_vivo", "in_vivo")
                    make_iso = cfg.get("dsi_make_isotropic", "0")
                    b_table = cfg.get("dsi_b_table", "auto")
                    optional = cfg.get("dsi_optional")
                    thread_count = cfg.get("num_processes", 1)
                    legacy = bool(cfg.get("dsi_legacy", False))
                    skip_motion_correction = bool(cfg.get("dsi_skip_motion_correction", False))

                    cli_str = (
                        f'dsi_main.py -i {_quote(currentFile[0])} '
                        f'-b {_quote(b_table)} '
                        f'-t {track_param_args} -r {_quote(recon_method)} '
                        f'-v {_quote(vivo)} -m {_quote(make_iso)} '
                        f'--thread-count {thread_count}'
                    )
                    if legacy:
                        cli_str += ' -l'
                    if skip_motion_correction:
                        cli_str += ' --skip-motion-correction'
                    if optional:
                        cli_str += ' -o ' + ' '.join(_quote(item) for item in optional)

                    os.chdir(cwd + '/3.2_DTIConnectivity')
                    command = f'python {cli_str}'
                    result = run_subprocess(command,dataFormat,step)
                    if result != 0:
                        errorList.append(result)
                else:
                    message = f'Could not find DWI input for DSI processing in {str(currentPath_wData)}';
                    logging.error(message)
                    errorList.append(message)
                os.chdir(cwd)
        else:
            message = 'The data folders'' names do not match anat, dwi, func or t2map';
            logging.error(message);
            errorList.append(message)
    else:
        message = 'The folder ' + dataFormat + ' does not exist in ' + str(currentPath_wData)
        logging.error(message)
        errorList.append(message)
    
    if errorList:
        return errorList
    else:
        return 0
    
 
def find(pattern, path):
    # This function finds all files with a specified fragment within
    # the given path
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def create_qc_reports(project_path, steps):
    requested_steps = set(steps)

    try:
        from helper_tools.batch_qc_reports import (
            build_bet_qc_report,
            build_registration_qc_report,
        )
    except Exception as exc:
        logging.warning("Could not import batch report tools: %s", exc)
        print(f"Report generation skipped: {exc}")
        return

    if "preprocess" in requested_steps:
        try:
            html_path, count = build_bet_qc_report(project_path, n_slices=10)
            if html_path:
                print(f"BET report written to {html_path} ({count} image(s))")
                logging.info("BET report written to %s (%s images)", html_path, count)
            else:
                print("BET report skipped: no BET files found.")
                logging.info("BET report skipped: no BET files found.")
        except Exception as exc:
            logging.warning("BET report generation failed: %s", exc)
            print(f"BET report generation failed: {exc}")

    if "registration" in requested_steps:
        try:
            html_path, count = build_registration_qc_report(project_path, n_slices=7)
            if html_path:
                print(f"Registration report written to {html_path} ({count} image(s))")
                logging.info("Registration report written to %s (%s images)", html_path, count)
            else:
                print("Registration report skipped: no BET/AnnoSplit_parental pairs found.")
                logging.info("Registration report skipped: no BET/AnnoSplit_parental pairs found.")
        except Exception as exc:
            logging.warning("Registration report generation failed: %s", exc)
            print(f"Registration report generation failed: {exc}")


def format_step_label(step):
    labels = {
        "preprocess": "Preprocessing",
        "registration": "Registration",
        "process": "Processing",
    }
    return labels.get(step, step.capitalize())


TQDM_BAR_FORMAT = (
    "{desc}: {percentage:3.0f}%|{bar}| "
    "{n_fmt}/{total_fmt} done | elapsed {elapsed} | remaining {remaining}"
)

if __name__ == "__main__":
    import argparse

    def parse_cpu_percent(value):
        cpu_count = multiprocessing.cpu_count()
        value = str(value).strip()

        if value.endswith("%"):
            value = value[:-1]

        try:
            percent = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "--cpu-percent must be a percentage from 1 to 100, e.g. 50 or 50%"
            )
        if percent <= 0 or percent > 100:
            raise argparse.ArgumentTypeError(
                "--cpu-percent must be greater than 0 and at most 100"
            )
        return max(1, int(cpu_count * percent / 100 + 0.5))

    def parse_cpu_cores(value):
        value = str(value).strip().lower()
        if value in {"min", "half", "max"}:
            return value
        try:
            cores = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "--cpu-cores must be min, half, max or a positive integer"
            )
        if cores < 1:
            raise argparse.ArgumentTypeError("--cpu-cores must be at least 1")
        cpu_count = multiprocessing.cpu_count()
        if cores > cpu_count:
            raise argparse.ArgumentTypeError(
                f"--cpu-cores must not exceed the available CPU cores ({cpu_count})"
            )
        return cores

    parser = argparse.ArgumentParser(
        description=(
            "Batch processing of all data (AIDAmri). "
            "Runs preprocessing, registration and processing steps for T2, DWI, fMRI and T2map.\n\n"
            "Example:\n"
            "python batchProc.py -i /path/to/proc_data -t anat dwi "
            "--t2-frac 0.1 --t2-bias-method mico "
            "--dwi-denoiser patch2self"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ============================================================
    # REQUIRED
    # ============================================================
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the parent project folder (e.g. proc_data)"
    )

    # ============================================================
    # GLOBAL / BATCH CONTROL
    # ============================================================
    batch = parser.add_argument_group("batch control")
    batch.add_argument(
        "-s", "--sessions",
        nargs="+",
        help="Process only selected sessions (e.g. Baseline P7 P14)"
    )
    batch.add_argument(
        "-t", "--data-types",
        nargs="+",
        choices=["anat", "dwi", "func", "t2map"],
        help="Data types to process (anat, dwi, func, t2map). Default: all"
    )
    batch.add_argument(
        "-d", "--debug-steps",
        dest="debug_steps",
        nargs="+",
        choices=["preprocess", "registration", "process"],
        help="Processing steps to run (preprocess registration process). Default: all"
    )
    batch.add_argument(
        "--slice-time-correction",
        action="store_true",
        help="Enable slice time correction for fMRI"
    )
    # ============================================================
    # CPU / PARALLELIZATION
    # ============================================================
    cpu = parser.add_argument_group("cpu / parallelization")
    cpu.add_argument(
        "-c", "--cpu-cores",
        default="half",
        type=parse_cpu_cores,
        help="CPU usage preset (min, half, max) or explicit number of parallel processes"
    )
    cpu.add_argument(
        "-p", "--cpu-percent",
        dest="cpu_percent",
        type=parse_cpu_percent,
        help="CPU percentage for parallel processes, e.g. 50 or 50%%"
    )

    # ============================================================
    # T2 PREPROCESSING (preProcessing_T2.py)
    # ============================================================
    t2 = parser.add_argument_group("T2 preprocessing (preProcessing_T2.py)")
    t2.add_argument(
        "--t2-bias-method",
        choices=["none", "mico", "ants"],
        type=str.lower,
        default = "mico",
        help="Bias field correction method for T2 (none, mico or ants). Default: mico"
    )
    t2.add_argument(
        "--t2-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help="Brain extraction method for T2: skip, bet or bet4animal. Default: bet"
    )

    t2.add_argument(
        "--t2-frac",
        type=float,
        help="BET fractional intensity threshold"
    )
    t2.add_argument(
        "--t2-radius",
        type=int,
        help="BET head radius in mm"
    )
    t2.add_argument(
        "--t2-gradient",
        type=float,
        help="BET horizontal gradient"
    )
    t2.add_argument(
        "--t2-center",
        nargs=3,
        type=int,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates"
    )
    t2.add_argument(
        "--t2-incidence-script",
        choices=["par", "detailed"],
        default="par",
        help="T2 incidence script for anat process: par (parental Atlas) runs getIncidenceSize_par.py, detailed runs getIncidenceSize.py. Default: par"
    )

    # ============================================================
    # DWI PREPROCESSING (preProcessing_DTI.py)
    # ============================================================
    dwi = parser.add_argument_group("DWI preprocessing (preProcessing_DTI.py)")
    dwi.add_argument(
        "--dwi-denoiser",
        choices=["patch2self"],
        type=str.lower,
        help="DWI denoising method"
    )
    dwi.add_argument(
        "--dwi-average-b0",
        action="store_true",
        help="Average b0 volumes before DWI processing"
    )
    dwi.add_argument(
        "--dwi-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help="Brain extraction method for DWI: skip, bet or bet4animal. Default: bet"
    )
    dwi.add_argument(
        "--dwi-skip-min-projection",
        action="store_true",
        help="Skip minimum intensity projection step"
    )
    dwi.add_argument(
        "--dwi-frac",
        type=float,
        help="BET fractional intensity threshold for DWI"
    )
    dwi.add_argument(
        "--dwi-radius",
        type=int,
        help="BET head radius (mm) for DWI"
    )
    dwi.add_argument(
        "--dwi-gradient",
        type=float,
        help="BET horizontal gradient for DWI"
    )
    dwi.add_argument(
        "--dwi-bias-method",
        choices=["none", "mico", "ants"],
        type=str.lower,
        default=None,
        help="Bias field correction for DWI: none, MICO or ANTs (default: None)"
    )

    # ============================================================
    # fMRI PREPROCESSING (preProcessing_fMRI.py)
    # ============================================================
    func = parser.add_argument_group("fMRI preprocessing (preProcessing_fMRI.py)")
    func.add_argument(
        "--func-bias-method",
        choices=["none", "ants"],
        type=str.lower,
        default=None,
        help="Bias field correction for fMRI: none or ANTs (default: None)"
    )
    func.add_argument(
        "--func-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help="Brain extraction method for fMRI preprocess/process: skip, bet or bet4animal. Default: bet"
    )
    func.add_argument(
        "--func-frac",
        type=float,
        help="BET fractional intensity threshold for fMRI"
    )
    func.add_argument(
        "--func-radius",
        type=int,
        help="BET head radius (mm) for fMRI"
    )
    func.add_argument(
        "--func-gradient",
        type=float,
        help="BET horizontal gradient for fMRI"
    )
    func.add_argument(
        "--func-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates for fMRI"
    )

    # ============================================================
    # T2MAP PREPROCESSING (preProcessing_T2MAP.py)
    # ============================================================
    t2map = parser.add_argument_group("T2map preprocessing (preProcessing_T2MAP.py)")
    t2map.add_argument(
        "--t2map-bet",
        choices=["skip", "bet", "bet4animal"],
        type=str.lower,
        default="bet",
        help="Brain extraction method for T2map: skip, bet or bet4animal. Default: bet"
    )
    t2map.add_argument(
        "--t2map-bias-method",
        choices=["none", "mico"],
        type=str.lower,
        default="mico",
        help='Biasfield correction method for T2map: none or mico. Default: mico'
    )
    t2map.add_argument(
        "--t2map-frac",
        type=float,
        help="BET fractional intensity threshold for T2map"
    )
    t2map.add_argument(
        "--t2map-radius",
        type=int,
        help="BET head radius (mm) for T2map"
    )
    t2map.add_argument(
        "--t2map-gradient",
        type=float,
        help="BET horizontal gradient for T2map"
    )
    t2map.add_argument(
        "--t2map-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="BET center in voxel coordinates for T2map"
    )

    # ============================================================
    # DSI STUDIO / TRACTOGRAPHY (dsi_main.py)
    # ============================================================
    dsi = parser.add_argument_group("DSI Studio / tractography (dsi_main.py)")
    dsi.add_argument(
        "--dsi-b-table",
        default="auto",
        help='Diffusion gradient source: "auto" or explicit b-table path'
    )
    dsi.add_argument(
        "--dsi-recon-method",
        default="dti",
        type=str.lower,
        choices=["dti", "gqi"],
        help="DSI reconstruction method"
    )
    dsi.add_argument(
        "--dsi-vivo",
        default="in_vivo",
        type=str.lower,
        choices=["in_vivo", "ex_vivo"],
        help="In vivo or ex vivo data (controls sampling length)"
    )
    dsi.add_argument(
        "--dsi-make-isotropic",
        default='0',
        help="Voxel size (mm) for isotropic resampling (0 = off, auto = header)"
    )
    dsi.add_argument(
        "--dsi-track-param",
        nargs="+",
        default="default",
        help="Tracking parameter preset or 8 custom values"
    )
    dsi.add_argument(
        "--dsi-skip-motion-correction",
        dest="dsi_skip_motion_correction",
        action="store_true",
        help="Skip slice-wise motion correction"
    )
    dsi.add_argument(
        "--dsi-legacy",
        action="store_true",
        help="Enable legacy .fib.gz / .src.gz support"
    )
    dsi.add_argument(
        "--dsi-optional",
        nargs="*",
        choices=["fa0", "nii_gz"],
        help="Optional dsi_main compatibility outputs (fa0, nii_gz)"
    )

    args = parser.parse_args()

    pathToData = args.input
    sessions = args.sessions
    
    #configurate the logging module
    log_file_path = os.path.join(pathToData, "batchproc_log.txt")
    logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', force =True)

    stc = args.slice_time_correction

    if args.data_types is None:
        data_types = ["anat", "dwi", "func", "t2map"]
    else:
        data_types = args.data_types

    if args.debug_steps is None:
        steps = ["preprocess","registration","process"]
    else:
        steps = args.debug_steps
    
    print('Entered information:')
    print(pathToData)
    print('data_types %s' % data_types)
    print('Slice time correction [%s]' % stc)
    print('Steps %s' % steps)
    print()

    all_files = findData(pathToData, sessions, data_types)

    num_processes = 1

    if isinstance(args.cpu_cores, int):
        num_processes = args.cpu_cores
    elif args.cpu_cores == "min":
        num_processes = 1
    elif args.cpu_cores == "half":
        num_processes = int(multiprocessing.cpu_count() / 2)
    elif args.cpu_cores == "max":
        num_processes = multiprocessing.cpu_count()

    print(args)
    
    if args.cpu_percent is not None:
        num_processes = args.cpu_percent
    
    print(f"Running with {num_processes} CPUs for the parallelization!")
    logging.info(f"Using {num_processes} CPUs for the parallelization")
    logging.info(f"Processing following datasets:\n{all_files}")
    # turns argparse.Namespace into a dict
    cfg = vars(args)
    cfg["num_processes"] = num_processes
    logging.info(
        "DSI settings: b_table=%s recon=%s vivo=%s make_isotropic=%s track_param=%s skip_motion_correction=%s legacy=%s optional=%s",
        args.dsi_b_table,
        args.dsi_recon_method,
        args.dsi_vivo,
        args.dsi_make_isotropic,
        args.dsi_track_param,
        args.dsi_skip_motion_correction,
        args.dsi_legacy,
        args.dsi_optional)

    for key, value in all_files.items():
        if value:
            error_list_all = []
            print()
            print(f"Entered {key} data: \n{value}")
            print()
            print(f"\nStarting {key} pipeline \33[5m...\33[0m")
            print()
            for step in steps: 
                error_list_step = []
                step_label = format_step_label(step)
                step_display = f"{step_label} {key} data"
                progress_bar = tqdm(
                    total=len(value),
                    desc=step_display,
                    unit="dataset",
                    bar_format=TQDM_BAR_FORMAT,
                )
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
                    futures = [executor.submit(executeScripts, path, key, step, cfg, stc) for path in value]

                    # --- collect errors robustly ---
                    flat_errors_step = []

                    for future in concurrent.futures.as_completed(futures):
                        progress_bar.update(1)

                        res = future.result()

                        # normalize result into a flat list
                        if res == 0 or res is None:
                            continue

                        if isinstance(res, list):
                            flat_errors_step.extend(res)
                        else:
                            flat_errors_step.append(res)

                    concurrent.futures.wait(futures)
                    progress_bar.close()

                    # keep a per-step and per-datatype summary
                    if not flat_errors_step:
                        print(f"{step_display}  \033[0;30;42m Complete \33[0m")
                    else:
                        print(f"{step_display}  \033[0;30;41m Incomplete \33[0m")
                        error_list_all.extend(flat_errors_step)

                    logging.info(f"{key} {step} processing completed")

            logging.error(f"Following errors were occurring: {error_list_all}")
            logging.info(f"{key} pipeline completed")

            if not error_list_all:
                print(f"\n{key} pipeline \033[0;30;42m COMPLETED \33[0m")
            else:
                print(f"\n{key} pipeline \033[0;30;41m INCOMPLETE \33[0m")
                print()
                for err in error_list_all:
                    if isinstance(err, tuple) and len(err) == 4:
                        sub, ses, dtype, stepname = err
                        print(
                            f"Error in sub: {sub} in session: {ses} in datatype: {dtype} and step: {stepname}. Check log.")
                    else:
                        # strings or unexpected types
                        print(f"Error: {err}")

    create_qc_reports(pathToData, steps)

 
