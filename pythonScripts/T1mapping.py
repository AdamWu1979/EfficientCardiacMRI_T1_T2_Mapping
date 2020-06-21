import sys
sys.path.append("C:\\Fani\\sourceCode\\fittingCodePy\\pythonScripts\\pyImage")

import pyImageDy

import dicom
import os
from os import listdir
from os.path import isfile, isdir, join
import numpy as np
from matplotlib import pyplot, cm
from matplotlib.widgets import Slider

mainDir = "C:\\XXX\\"
subjDirs = ["01","02","03","04","05","06","07"]
dirM = "PRE-T1"
ext = "IMA"
outputMainDir = "C:\\XXX\\Results\\MLeigen\\"
namesOutput = ["A","B","TI","dataT1star","dataXsquare","minIndexImg","statusImg","numIterImg"]

saveFlag = False

if(saveFlag):
    if not os.path.exists(outputMainDir):
        os.makedirs(outputMainDir)


for subj in subjDirs:
    pathC = join(mainDir,subj,dirM)
    print(pathC)
    outputDir = join(outputMainDir,subj)

    if(saveFlag):
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        
    outputDir = join(outputDir,dirM)

    if(saveFlag):
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)    
    print("sort out files")
    onlyfiles = [f for f in listdir(pathC) if (isfile(join(pathC, f)) and f.endswith(ext))]

    headerInfo = list()
    seriesNo = list()
    invTime = list()
    DicomDataS = list()
    for cc, f in enumerate(onlyfiles):
        fileN = join(pathC,f)
        print(fileN)
        headerInfo.append( dicom.read_file(fileN) )
        seriesNo.append(float(headerInfo[cc].SeriesNumber))
        invTime.append(float(headerInfo[cc].InversionTime))

    
    slices = list(set(seriesNo))
    indS = sorted(range(len(slices)), key=lambda k: slices[k])
    slices = [slices[i] for i in indS]
    for ss in slices:
        indS = [i for i, x in enumerate(seriesNo) if x == ss]
        invTimeS = [invTime[i] for i in indS ]
        indSS = sorted(range(len(invTimeS)), key=lambda k: invTimeS[k])
        indF = [indS[i] for i in indSS]
        invTimeF = [invTime[i] for i in indF]

        DicomDataS = list()
        ImgDimsAll = list()
        #intantiate the c++ class for fitting
        myImagePack = pyImageDy.Images()
        for ii, vv in enumerate(indF):
            dataL = (headerInfo[vv].pixel_array).tolist()
            DicomDataS.append( [float(item) for sublist in dataL for item in sublist] )
            dims = [ int( ((headerInfo[vv].pixel_array).shape )[i] ) for i,xx in enumerate(( headerInfo[vv].pixel_array).shape)]
            ImgDimsAll.append( [len((headerInfo[vv].pixel_array).shape)] + dims  ) 
            #initialise the c++ class
            myImagePack.initImagePy( DicomDataS[ii], ImgDimsAll[ii], invTimeF[ii] )
        
        #apply fitting function
        print("Fitting...")
        myImagePack.fitImagesPy()
        print("...Done")

        #save files
        for ii, name in enumerate(namesOutput):
            outImg = myImagePack.getImagePy(ii)
            outImg = np.asarray(outImg)
            finTest = np.isfinite(outImg)
            ind = np.where(finTest!=True)
            outImg[ind] = 0
            outImg = outImg.reshape(ImgDimsAll[0][1:3])
            outputFileName = join(outputDir,str(abs(int(ss)))+"_"+name)
            print(outputFileName)
            if(saveFlag):
                np.save(outputFileName,outImg,allow_pickle=False)
            

            
    


