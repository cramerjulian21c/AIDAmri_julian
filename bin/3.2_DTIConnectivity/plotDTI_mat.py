"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne



"""

import matplotlib.pyplot as plt
import os, sys
import numpy as np
import scipy.io as sio

np.seterr(divide='ignore', invalid='ignore')
import seaborn as sns

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
DEFAULT_SIGMA_LABEL_FILE = os.path.join(
    REPO_ROOT,
    "lib",
    "sigma",
    "SIGMA_InVivo_Anatomical_Brain_Atlas_Labels.txt",
)


def intersect_mtlb(a, b):
    a1, ia = np.unique(a, return_index=True)
    b1, ib = np.unique(b, return_index=True)
    aux = np.concatenate((a1, b1))
    aux.sort()
    c = aux[:-1][aux[1:] == aux[:-1]]
    return ia[np.isin(a1, c)]


def parse_label_line(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None

    parts = stripped.split()
    try:
        label_id = int(float(parts[0]))
    except (IndexError, ValueError):
        return None

    if '"' in stripped:
        label_name = stripped.split('"', 1)[1].rsplit('"', 1)[0]
    elif "\t" in stripped:
        tab_parts = stripped.split("\t")
        label_name = tab_parts[1].strip() if len(tab_parts) > 1 else ""
    else:
        label_name = " ".join(parts[1:])

    return label_id, label_name


def getRefLabels(label_file=DEFAULT_SIGMA_LABEL_FILE):
    ref_labels = []
    with open(label_file, "r") as label_handle:
        for line in label_handle:
            parsed_label = parse_label_line(line)
            if parsed_label is None:
                continue

            label_id, label_name = parsed_label
            if label_id > 0:
                ref_labels.append(label_name)

    if len(ref_labels) == 0:
        sys.exit("Error: No labels could be read from '%s'." % (label_file,))

    return np.array(ref_labels)


def matrixMaker(inputPath):
    # Read pass and end
    if "pass" in inputPath:
        matData = sio.loadmat(inputPath)
        connectivityPass = matData['connectivity']
        matData = sio.loadmat(inputPath.replace('.pass.', '.end.'))
        connectivityEnd = matData['connectivity']
        connectivity = connectivityEnd + connectivityPass

    elif "end" in inputPath:
        matData = sio.loadmat(inputPath)
        connectivityEnd = matData['connectivity']
        matData = sio.loadmat(inputPath.replace('.end.', '.pass.'))
        connectivityPass = matData['connectivity']
        connectivity = connectivityEnd + connectivityPass

    else:
        sys.exit("Error: %s path do not conatain path or end data." % (inputPath,))

    labels = matData['name']
    tempLabels = ""
    labels = tempLabels.join([chr(a) for a in labels[0]]).split('\n')

    # Get reference Labels
    refLabels = getRefLabels()

    # Intersection between ref and cur labels
    ia = intersect_mtlb(refLabels, labels)
    missingLabels = np.setdiff1d(np.arange(1, len(refLabels)), ia)

    # Adapt labels to pyplot
    labels = [s.replace('_', ' ') for s in labels]

    zeroVec = np.zeros([len(refLabels), len(refLabels)])
    zeroVec[np.ix_(np.sort(ia), np.sort(ia))] = connectivity

    connectivityFilled = zeroVec

    fig, ax = plt.subplots()

    sns.heatmap(connectivityFilled)
    ax.axis('tight')

    # Set labels
    ax.set(xticks=np.arange(len(labels)), xticklabels=labels,
           yticks=np.arange(len(labels)), yticklabels=labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    ax.set_title("DTI connectivity between SIGMA regions")
    plt.show()

    return connectivity


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Visualize mat file of DTI ')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--inputMat', help='file name: DTI mat-File')
    args = parser.parse_args()

    inputPath = None
    if args.inputMat is not None and args.inputMat is not None:
        inputPath = args.inputMat
    if not os.path.exists(inputPath):
        sys.exit("Error: %s path is not an existing directory." % (args.inputPath,))

    inputPath = args.inputMat

    # generate Matrix
    matrixMaker(inputPath)
