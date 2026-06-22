'''
Created on 08.04.2019
Updated: 26.09.2020
@author: Niklas Pallast and Markus Aswendt

process all DTI data
'''


import glob
import os
import numpy as np

def findData(path):


    regAtlas_list = []
    patterns = [
        os.path.join(path, '**', 'dwi', '*_AnnoSplit_parental.nii.gz'),
        os.path.join(path, '**', 'DTI', '*_AnnoSplit_parental.nii.gz'),
    ]
    for pattern in patterns:
        for filename in glob.iglob(pattern, recursive=True):
            regAtlas_list.append(filename)



    return sorted(set(regAtlas_list))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Find all related DTI data')
    parser.add_argument('-p','--pathData', help='Path to study')

    args = parser.parse_args()

    pathData = args.pathData

    listAtlas = findData(pathData)
    print(listAtlas)
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.md.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.fa0.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.rd.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
    for i in range(np.size(listAtlas)):
        print(listAtlas[i])
        curPath = os.path.dirname(listAtlas[i])
        dti = glob.glob(curPath + '/*.ad.nii.gz')
        if dti:
            print('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
            os.system('python DTIdata_extract.py ' + dti[0] + ' ' + listAtlas[i])
